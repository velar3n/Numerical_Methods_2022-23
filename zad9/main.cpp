#include <iostream>

int main() {

    float x[8] = {-0.875, -0.625, -0.375, -0.125, 0.125, 0.375, 0.625, 0.875};
    float a[8], h[7], xa[7], xl[8], xu[8], xz[8], b[8], c[8], d[8];
    int n = 7, m = 6;

    for (int i=0; i <= n; i++) {
        a[i] = 1 / (1 + (5*x[i]*x[i]));
    }

    for (int i=0; i <= m; i++) {
        h[i] = x[i+1] - x[i];
    }

    for (int i=1; i <= m; i++){
        xa[i] = 3*(a[i+1]*h[i-1] - a[i]*(x[i+1] - x[i-1]) + a[i-1]*h[i]) / (h[i]*h[i-1]);
    }

    xl[0] = 1;
    xu[0] = 0;
    xz[0] = 0;

    for (int i=1; i <= m; i++) {
        xl[i] = 2 * (x[i+1] - x[i-1]) - h[i-1]*xu[i-1];
        xu[i] = h[i] / xl[i];
        xz[i] = (xa[i] - h[i-1]*xz[i-1]) / xl[i];
    }

    xl[8] = 1;
    xz[8] = 0;
    c[8] = 0;

    for (int i=0; i<=m; i++) {
        int j = m - i;
        c[j] = xz[j] - xu[j]*c[j+1];
        b[j] = (a[j+1] - a[j]) / h[j] - h[j]*(c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    std::cout<<"\n\nWspółczynniki wielomianów:\n";
    std::cout<<"     a(i)        b(i)          c(i)          d(i)\n";
    for(int i=0; i <= m; i++){
        printf("%11.8f  %11.8f  %11.8f  %11.8f \n",a[i], b[i], c[i], d[i]);
    }
    std::cout<<std::endl;

    for (int i=0; i < n; i++) {
        std::cout<<"Wielomian dla przedziału ("<<x[i]<<")-("<<x[i+1]<<"):   "<<a[i]<<" + "<<b[i]<<"(x - "<<x[i]<<") + "<<c[i]<<"(x - "<<x[i]<<")^2 + "<<d[i]<<"(x - "<<x[i]<<")^3"<<std::endl;
    }
    std::cout<<std::endl<<std::endl;


    return 0;
}