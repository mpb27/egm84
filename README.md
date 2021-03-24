# EGM84 Gravitational Model

Coefficients are from:

https://geographiclib.sourceforge.io/html/gravity.html

Direct link to EGM84 n=m=180: [Download ZIP](/references/egm84.zip)

EGM84 n=m=18 coefficients (see page 194): [Download PDF](/references/Supplement_to_Department_of_Defense_Worl.pdf)

GeographicLib coefficient file format:

```c
struct {
    char name[8];
    int n; // or m?
    int m; // or n?
    double values[]; // all values
}
```
