/*

kwadratura - wzór na całkowanie przybliżone, uzyskany przez całkowanie odpowiednich wielomianów interpolacyjnych
	z całkowania interolacyjnego Lagrange'a --> kwadratury Newtona-Cotesa (zamknięte jeśli krańce przedaiłu całkowania to węzły interpolacji)
	metoda trapezów - interpolacja wielomianu niskiego stopnia, jest wzorek do tego 

Całkowanie po przedziałach nieskończonych - trzeba uważać by nie wyliczać rozbieżnej całki. Funkcja podcałkowa, musi w nieskończoności zmierzać
	dostatecznie szybko do 0, aby całka istniała (wzór 9 z zadania). Potem całkę I1 obliczamy metodą Romberga / trapezów numerycznie, a I ogon 
	szacujemy analitycznie

Metoda Romberga - wykorzystuje błąd metody trapezów - że zawiera on wyłącznie parzyste potęgi średnicy podziału

	Trzeba znaleźć A, I1, całkę obliczoną z dokładnością 10^-7

	Obliczenia przerywamy, gdy dwa kolejne wyrazy diagonalne staną się sobie równe (z zadaną dokładnością)


	A powinno wyjść 17 według starego wykładu


	1. Obliczenie całki za pomocą złożonego wzoru trapezów (IN ~ h(...)) i metody Pomberga --> aż do równości wyrazów z określonym błędem
		czyli A0k > A0k+1 o mniej niż błąd --> szybsze niż stosowanie samego wzoru trapezów

*/

#include <iostream>
#include <math.h>
#include <vector>
#include <functional>

// M_PI = pi 
// M_E = e

double trapezoidalIntegral(double a, double b, int n, const std::function<double (double)> &f) {        // metoda trapezów

    double h, x1, x2, trapezoidal_integral = 0;
    h = (b-a)/n;

    for(int step = 0; step < n; step++) {
        x1 = a + step*h;
        x2 = a + (step+1)*h;
        trapezoidal_integral += 0.5*(x2-x1)*(f(x1) + f(x2));
    }

    return trapezoidal_integral;
}

double rombergIntegral(double a, double b, int max, const std::function<double (double)> &f) {

    std::vector<std::vector<double>> romberg_integral(max, std::vector<double>(max));
    romberg_integral[0][0] = trapezoidalIntegral(a, b, 1, f);

    double h, trapezodial, number1, number2;
    int end;
    h = b-a;

    for(int step = 1; step < max; step++) {                     // metoda Romberga

        h *= 0.5;

        trapezodial = 0;
        int stepEnd = pow(2, step - 1);
        for(int step1 = 1; step1 <= stepEnd; step1++) {         // metoda trapezów
            number1 = (2*step1 - 1)*h;
            trapezodial += f(a + number1);
        }
        romberg_integral[step][0] = 0.5*romberg_integral[step-1][0] + trapezodial*h;

        for(int step2 = 1; step2 <= step; step2++) {            // ekstrapolacja
            number2 = pow(4, step2);
            romberg_integral[step][step2] = (number2*romberg_integral[step][step2-1] - romberg_integral[step-1][step2-1])/(number2-1);
        }

        // przerwanie działań przy określonym zakresie błędu
        if(abs(romberg_integral[step][step] - romberg_integral[step-1][step-1]) <= 0.0000001){
            end = step;
            printf("\nWynik metody Romberga, przy k = %i wynosi: %11.8f \n\n", step, romberg_integral[step][step]);
            break;
        }
    }

    return romberg_integral[end][end];
}

int main() {

    int max = 30;                                           // ustalenie górnej granicy liczby obliczeń funkcji
    double granica_gorna = 1;

    auto func = [](double x) {                              // funkcja z której wyliczamy całkę
        return sin(M_PI*(1 + std::sqrt(x))/(1 + (x*x)))*std::pow(M_E, -x);
    };

    while(pow(M_E, -granica_gorna) >= 0.0000001) {          //onliczenie A, czyli górnej granicy całki I(1)
        granica_gorna = granica_gorna + 1;
    }
    std::cout<<"\nA = "<<granica_gorna<<std::endl;

    double romberg = rombergIntegral(0, granica_gorna, max, func);  // metoda Romberga (z użyciem metody trapezów)

    std::cout<<"Przy e^(-A) = I(ogon) ~ 0.414 * 10^(-7), czyli jest mniejsze od zakresu błędu i przez to pomijalne --> I = I(1) + I(ogon) = I(1) = ";
    printf("%11.8f \n\n", romberg);

	return 0;
}