#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
using namespace std;
typedef complex<double> dcomplex;

// derivative function prototype
void deriv(complex<double> (*f)(const complex<double>& ) ,const double, double&);

	   
// help function for output
void print(double x, double exdf,double df){
 
  static bool first=true;
  const unsigned w=25;
 
  if(first){
    cout << scientific << setprecision(16);
    cout << setw(w) << "x" << setw(w) << "exact f'(x)" << setw(w);
    cout << "num. f'(x)" << setw(w) << "error" << endl;
  }
 
  cout << setw(w) << x << setw(w) << exdf;
  cout << setw(w) << df << setw(w) << exdf-df << endl;
  first =false;
}


//-------MAIN-----------------
int main(){
 
  double x,der;
  x = 0.0;
 
  for(unsigned i = 0; i < 20; ++i){
    deriv(tan,x,der); // compute numerical derivative
    print(x,1.0/(cos(x)*cos(x)),der);
    deriv(cos,x,der); // compute numerical derivative
    print(x,-sin(x),der);
    x = x + 0.1;
  }
 
  return 0;
}



//-----Derivaattan laskeminen-----------
void deriv(complex<double> (*fun)(const complex<double>& ) , const double x, double& der){
	
	double h = 1.0E-15;
	dcomplex ch(0.0, h);
	static bool first=true;
	  
	if( first ) {cout << "ch = 0 + ih = " << ch << endl; first = false; }
	der = imag( fun(x + ch) ) / h;
	//der = real((fun(x + h) - fun(x) )) / h;
	
}
