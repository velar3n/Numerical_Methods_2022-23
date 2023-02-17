#include <iostream>			//użycie interpolacji lagrangea do obliczenia wykresu funckji
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

using namespace std;

vector<double>* lagrange_c(const vector<double> x, const vector<double> fx) {
    auto length = x.size();
    auto res = new vector<double> (length, 0);
    
    for (int i=0; i < 8; i++) {
        vector<double> tmpcoeffs (length, 0);
        tmpcoeffs[0] = fx[i];
        double prod = 1;
        for(int j=0; j < 8; j++) {
            if (x[i] == x[j]) {
            	continue;
            }
            prod *= x[i] - x[j];
            double precedent = 0;
            for (auto resptr = tmpcoeffs.begin(); resptr < tmpcoeffs.end(); resptr++) {
                double newres = (*resptr) * (-x[j]) + precedent;
                precedent = *resptr;
                *resptr = newres;
            }
        }
        transform(res->begin(), res->end(),
                  tmpcoeffs.begin(),
                  res->begin(),
                  [=] (double oldcoeff, double add) {return oldcoeff+add/prod;}
                  );
    }
    return res;
}

int main(){

	std::ofstream file ("points.txt");
	if(!file.is_open()){
       std::cout<<"Error in creating file!!!";
       return 0;
   }

	vector<double> x = {-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0}; 
	vector<double> fx = {1.13092041015625, 2.3203125, 1.92840576171875, 1.0, 0.05548095703125, -0.6015625, -0.75250244140625, 0.0};

	for(int i=0; i < 8; i++){
		file<<x[i]<<" "<<fx[i]<<endl;
	}
    file.close();

	cout<<"Współczynniki przy x: "<<endl;
	int i = 0;
	for (double coefficient : *lagrange_c(x,fx)) {
		if(abs(coefficient) < 0.00009999){
			coefficient = 0;
		}
        cout<<"x^"<<i<<" = "<<setprecision(4)<<coefficient<<endl;
        i++;
    }

	return 0;
}