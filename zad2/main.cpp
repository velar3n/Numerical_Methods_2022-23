#include <iostream>
#include <fstream>
#include <cstdlib>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>

using Eigen::VectorXi;
using namespace Eigen;

int main(){


	//roziwązywanie metodą Gaussa-Seidlera
	VectorXd diagonalA(124);					// stworzenie macierzy jako wektorów
	VectorXd diagonalB(127);
	VectorXd diagonalC(128);
	VectorXd diagonalD(127);
	VectorXd diagonalE(124);
	VectorXd b(128);
	VectorXd x(128);

	diagonalA.setOnes();							// zapełnienie macierzy
	diagonalB.setOnes();
	diagonalD.setOnes();
	diagonalE.setOnes();
	b.setOnes();

	for(int i=0; i<128; i++){
		diagonalC(i) = 4;
	}

	double norm1;
	double norm2 = 0;
	double convergence;

	std::ofstream file ("convergenceA.txt");
	if(!file.is_open()){
       std::cout<<"Error in creating file!!!";
       return 0;
   }

	for (int i=1; i<98; i++){						// rozwiązanie metodą Gaussa-Seidela - tyle samo iteracji co w metodzie z podbunktu b)
		
		x[0] = b[0] - x[1] - x[4];
		x[0] /= diagonalC[0];
		x[1] = b[1] - x[0] - x[2] - x[5];
		x[1] /= diagonalC[1];
		x[2] = b[2] - x[1] - x[3] - x[6];
		x[2] /= diagonalC[2];
		x[3] = b[3] - x[2] - x[4] - x[7];
		x[3] /= diagonalC[3];
		
		for (int j=4; j<124; j++) {
			x[j] = b[j] - x[j-4] - x[j-1] - x[j+1] - x[j+4];
			x[j] /= diagonalC[j];
		}
		
		x[124] = b[124] - x[120] - x[123] - x[125];
		x[124] /= diagonalC[124];
		x[125] = b[125] - x[121] - x[124] - x[126];
		x[125] /= diagonalC[125];
		x[126] = b[126] - x[122] - x[125] - x[127];
		x[126] /= diagonalC[126];
		x[127] = b[127] - x[123] - x[126];
		x[127] /= diagonalC[127];

		norm1 = norm2;
		norm2 = x.squaredNorm();

		if(i != 1){
			convergence = abs(norm2 - norm1);
			file<<i<<" "<<convergence<<"\n";
		}
	}
	file.close();


	std::cout<<"Wynik dla metody Gaussa-Seidela: "<<std::endl;
	for (int i=0; i<128; i++) {
		std::cout<<x(i)<<" ";
	}
	std::cout<<std::endl<<std::endl;






   Eigen::VectorXd x2(128);                     // storzenie wektorów x oraz b
   Eigen::VectorXd b2(128);

   Eigen::SparseMatrix<double> A(128,128);            //stworzenie podanej macierzy rzadkiej
   A.reserve(VectorXi::Constant(128,6));

   for(int i = 0; i < 128; i++){
      A.coeffRef(i, i) = 4.0;
      if(i != 127){
         A.coeffRef(i+1, i) = 1.0;
         A.coeffRef(i, i+1) = 1.0;
      }
      if(i < 124){
         A.coeffRef(i+4, i) = 1.0;
      }
      if(i > 3){
      	A.coeffRef(i-4, i) = 1.0;
      }
      b2(i) = 1;                             // wektor b - ustawienie wartości na "1"
   }
   A.makeCompressed();

   norm2 = 0;

   std::ofstream file2 ("convergenceB.txt");
	if(!file2.is_open()){
       std::cout<<"Error in creating file!!!";
       return 0;
   }

   LeastSquaresConjugateGradient<SparseMatrix<double> > cg;
   for(int i = 0; i < 98; i++){
	   cg.setMaxIterations(i);
		cg.compute(A);
		x2 = cg.solve(b2);

		norm1 = norm2;
		norm2 = x2.squaredNorm();

		if(i != 0){
			convergence = abs(norm2 - norm1);
			file2<<i<<" "<<convergence<<"\n";
		}
		/*b2.setOnes();
		x2 = cg.solve(b2);*/
   }

	file.close();

	std::cout<<"Wynik dla metody gradientów sprężonych: "<<std::endl;
	for (int i=0; i<128; i++) {
		std::cout<<x2(i)<<" ";
	}
	std::cout<<std::endl<<std::endl;

	return 0;

}