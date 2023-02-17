#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

int main(){

	Eigen::MatrixXd A(6, 6);
	Eigen::MatrixXd B(6, 6);

	A << 19/12,  13/12,  5/6, 5/6,  13/12, -17/12,
		 13/12,  13/12,  5/6, 5/6,  -11/12, 13/12,
		 5/6,    5/6,    5/6, -1/6, 5/6,    5/6,
		 5/6,    5/6,   -1/6, 5/6,  5/6,    5/6,
		 13/12,  -11/12, 5/6, 5/6,  13/12,  13/12,
		 -17/12, 13/12,  5/6, 5/6,  13/12,  19/12;

	Eigen::Tridiagonalization<Eigen::MatrixXd> tri;
	tri.compute(A);
	std::cout<<"The matrix T in the tridiagonal decomposition of A is: "<<std::endl;
	std::cout<<tri.matrixT()<<std::endl;

	B = tri.matrixT();
	Eigen::EigenSolver<Eigen::MatrixXd> es(B);
	std::cout<<"The eigenvalues of A are:"<<std::endl<<es.eigenvalues()<<std::endl;

	return 0;
}


