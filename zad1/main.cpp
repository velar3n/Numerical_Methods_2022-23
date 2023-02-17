/*
	pierwsza macierz - faktoyzacja LU - metoda Thomasa, ponieważ wtedy faktoryzacja LU ma złożoność O(N) --> evtl QR z obrotami Givensa, wtedy faktoryzacja również jest liniowa
		Normlanie faktoryzacja LU ma złożoność O(N3), ale uwzględniając jej strukturę, czyli to, że jest trójdiagonalna, da się to zrobić w czasie liniowym.
		Można też Dolittlea --> tworzymy dwie macierze trójkątne, i rozwiązujemy 2 równania - pierwsze przez forward substitution, drugie przez backsubstitution

	druga macierz - tu nie możemy zastosować algorytmu Thomasa, ponieważ n.ie jest to macierz trójdiagonalna. Jest to jednak w dalszym ciągu macierz symetryczna i rzadka,
		więc zastosowany zostanie faktoryzacja LU (algprytm Dolittle'a) --> używamy biblioteki eigen --> sparse LU, czyli algorytm LU dla macierzy rz
*/


// Natalia Kiełbasa - metody numeryczne 
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/SparseLU> 
#include <Eigen/SparseCore>
#include "thomas_algorithm.h"

using Eigen::VectorXi;

int main(){

	int  size = 7;									// rozmiar macierzy
	double a[7] = { 0, 1, 1, 1, 1, 1, 1 };				// diagonala górna
	double b[7] = { 3, 4, 4, 4, 4, 4, 3 };				// diagonala
	double c[7] = { 1, 1, 1, 1, 1, 1, 0 };				// diagonala dolna
	double d[7] = { 1, 2, 3, 4, 5, 6, 7 };				// kolumna wyników równań
	std::string XNumbers[7] = { "x1: ", "x2: ", "x3: ", "x4: ", "x5: ", "x6: ", "x7: " };


	// rozwiązanie macierzy algorytmem Thomasa z pliku thomas_algorithm.h
	thomas_algorithm(a, b, c, d, size);
	std::cout<<"\nWyniki dla macierzy a) rozwiąznej metodą Thomasa: "<<std::endl;
	for (int i = 0; i < size; i++) {
		std::cout<<XNumbers[i]<<d[i]<<std::endl;
	}




	// rozwiązanie macierzy z zastosowaniem faktoryzacji LU dla macierzy rzadkich z Eigen
	Eigen::VectorXd x(7);							//stworzenie wektora z szukanymi x-ami
	Eigen::VectorXd b2(7);
	Eigen::SparseMatrix<double> A(7,7);				//stworzenie macierzy rzadkiej o wymiarach 7x7

	A.reserve(VectorXi::Constant(7,4));				//rezerwacja pamięci dla elementów niezeroych w macierzy
	double count = 1.0;
	for(int i = 0; i < 7; i++){
		b2(i) = count;								//zapełnienie wektora b2 wynikami równań
		A.coeffRef(i, i) += 4.0;					//zapełnienie macierzy rzadkiej
		if(i != 6){
			A.coeffRef(i+1, i) += 1.0;
			A.coeffRef(i, i+1) += 1.0;
		} else {
			A.coeffRef(0, i) += 1.0;
			A.coeffRef(i, 0) += 1.0;
		}
		count++;
	}
	A.makeCompressed(); 

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;					//rozwiązanie równania

	solver.analyzePattern(A);
	solver.factorize(A);
	x = solver.solve(b2);

	std::cout<<"\nWyniki dla macierzy b) rozwiąznej faktoryzacją LU dla macierzy rzadkich: "<<std::endl;
	for (int i = 0; i < size; i++) {
		std::cout<<XNumbers[i]<<x(i)<<std::endl;
	}

	return 0;
}