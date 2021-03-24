%% Read EGM-84 coefficients.

fid = fopen('c:\Temp\egm84.egm.cof', 'r');
name = fread(fid, 8, 'uint8=>char')'
num_n = fread(fid, 1, 'uint32')
num_m = fread(fid, 1, 'uint32')
coeffs = fread(fid, 'double');
fclose(fid);

% Last coefficient is 8xFF, ie NaN.
assert(isnan(coeffs(end)), 'Final value is not NaN');
coeffs = coeffs(1:end-1);

% Coefficents should be n+1,m+1.
assert((num_n+1)*(num_m+1) == length(coeffs));
%coeffs = reshape(coeffs, [num_n+1, num_m+1]);

% Create lambda for easy access to C coefficents, C(n,m):
C = @(n,m) coeffs( m*(num_n+1) - m*(m+1)/2 + n + 1 );

% Check against a few known coefficients from PDF.
ASSERT_C = @(n,m,v) ...
    assert(C(n,m) == v, 'Coefficient C(%d,%d) is incorrect.', n, m);

ASSERT_C( 2, 0, -0.48416685E-03)
ASSERT_C( 3, 0, +0.95706390E-06)
ASSERT_C(18, 0, +0.10196218E-07)
ASSERT_C( 2, 2, +0.24395796E-05)
ASSERT_C(16,12, +0.18888013E-07)
ASSERT_C(18,18, +0.11121796E-08)

% Create a lambda for access to S coefficients, S(n,m):
S = @(n,m) coeffs( 16471 + (m-1)*(num_n+1) - (m-1)*m/2 - m + n + 1 );

% Check against a few known coefficients from PDF.
ASSERT_S = @(n,m,v) ...
    assert(S(n,m) == v, 'Coefficient S(%d,%d) is incorrect.', n, m);

ASSERT_S( 1, 1,  0.0);
ASSERT_S( 2, 2, -0.13979548E-05)
ASSERT_S( 3, 3, +0.14152388E-05)
ASSERT_S( 8, 8, +0.12210258E-06)
ASSERT_S( 6, 1, +0.32780040E-07)
ASSERT_S(16,12, +0.46949615E-08)
ASSERT_S(18,18, -0.94806182E-08)


%% Compute WGS-84 Geoid Height (pg 220).

% Distance from the earth's geocentric center.
r = 6371008.7714;

% Setup constants.
GM = 3986005E+08;
a = 6378137;

% Normalized Legendre polynomial. (4-3) 
%  - Implemented as a recursive funtion from PDF 44310926 pg 14..17
% (see bottom of script for function)

% Gamma, theoretical gravity.
gamma_e = 978032.67714;
k = 0.00193185138639;
e2 = 0.00669437999013;
gamma = @(lat) gamma_e * (1 + k*sin(lat)^2) / sqrt(1 - e2*sin(lat)^2);

Nxx = @(n,m,lat,lon) ...
    (a/r)^n * ( C(n,m)*cos(m*lon) + S(n,m)*sin(m*lon) ) * Pnm(n,m,sin(lat));
Nx  = @(n,lat,lon) ...
    sum(arrayfun(@(m) Nxx(n,m,lat,lon), 0:n));  % m = 0..n
N   = @(lat,lon) ...
    GM / (r*gamma(lat)) ...
    * sum(arrayfun(@(n) Nx(n,lat,lon), 2:18)); % n = 2..18

U   = @(lat,lon) ...
    GM / r ...
    * sum(arrayfun(@(n) Nx(n,lat,lon), 2:18)); % n = 2..18



%% How I determined the indexing into data array.
% First row is m=0, n=0,1,2,3...180    (181 values)    1..181
% Second row is m=1, n=1,2,3,4,5,6,... (180 values)  182..361
% Third row is m=2, n=2,3,4,5,...      (179 values)  362..540
% Fourth row is m=3, n=3,4,5,6,...     (178 values)  541..718
% Offset is Gauss's sum on m.
% Offset = 181*m - m*(m+1)/2 + n + 1
Cindex = @(m,n) m*(num_n+1) - m*(m+1)/2 + n + 1;

% Now find the S values in the data array:
% S(3,3) is at index 16831

% C(180,180) is at 16471

% S(1,1) is at index 16472 and is 0
% S(2,1) is at index 16473 and is 0.
% S(3,1) is at index 16474
% S(4,1) is at index 16475

% S(1,2) is at index 16651? NO its a value. last value - S(180,1)
% S(2,2) is at index 16652
% S(3,2) is at index 16653
% S(4,2) is at index 16654

% First row is m=1 n=1,2,3,...   (180 values)  16472..16651
% Second row is m=2, n=2,3,4,... (179 values)  16652..16830
% Third row is m=3, n=3,4,5,...  (178 values)  16831..
Sindex = @(n,m) 16471 + (m-1)*(num_n+1) - (m-1)*(m)/2 - m + (n-1) + 1 + 1;



function x  = Pnm(n,m,sin_lat)
    v = sin_lat;
    if (n == 0 && m == 0)                       % (II.18)
        x = 1;  
    elseif (n == 1 && m == 1)                   % (II.19)
        x = sqrt(1-v^2) * sqrt(3) * 1.0;
    elseif (n == m)                             % (II.20)
        x = sqrt(1-v^2) * sqrt((2*m+1)/(2*m)) * Pnm(m-1,m-1,v);
    elseif ((n - m) == 1)                       % (II.21)
        x = v * sqrt(2*m+1) * Pnm(m,m,v);
    else                                        % (II.22)
        e1 = sqrt( (2*n+1)/(2*n-1) * (n-m)/(n+m) );
        e2 = sqrt( (2*n+1)/(2*n-3) * (n-m)/(n+m) * (n-m-1)/(n+m-1) );        
        x = 1/(n-m)*( ...
              (2*n-1)*v*e1*Pnm(n-1,m,v)  ...
            - (n+m-1)*e2*Pnm(n-2,m,v)    ...
            );
    end
end


