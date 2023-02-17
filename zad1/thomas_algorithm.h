// Natalia Kiełbasa - metody numeryczne 
#include <iostream>			// f - wyniki

void thomas_algorithm(double* a, double* b, double* c, double* d, int size) {

    size--;
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < size; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[size] = (d[size] - a[size]*d[size-1]) / (b[size] - a[size]*c[size-1]);

    for (int i = size; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}