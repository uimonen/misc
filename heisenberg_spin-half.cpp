/***************************************************************************
* Heisenberg  spin-1/2 model
**************************************************************************/
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/iteration.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>
#include <iostream>
#include <bitset>

typedef std::complex<double> value_type;
typedef boost::numeric::ublas::hermitian_matrix<value_type, boost::numeric::ublas::lower> Matrix;
typedef boost::numeric::ublas::vector<value_type> Vector;
typedef ietl::vectorspace<Vector> Vecspace;
typedef boost::lagged_fibonacci607 Gen;

void FillHermitianMatrix(Matrix &, int);
void FillHermitianMatrix2(Matrix &);
void PrintEigenvalues(ietl::lanczos<Matrix,Vecspace> & );
void PrintEigenvectors(std::vector<Vector> &, ietl::Info<double> &);

int main() {

    // Creation of a sample matrix:
    int length = 10;              // length of the chain
    int N=(int)pow(2,length);
    //Matrix mat2(length, length);
    Matrix mat(N, N);

    //Filling the matrix of size N
    //FillHermitianMatrix2(mat2);
    FillHermitianMatrix(mat, length);

    Vecspace vec(N);
    Gen mygen;
    ietl::lanczos<Matrix,Vecspace> lanczos(mat,vec);

    int max_iter = 50;
    std::cout << "Computation of eigenvalues with fixed size of T-matrix\n\n";
    std::cout << "-----------------------------------\n\n";

    std::vector<double> eigen;

    ietl::fixed_lanczos_iteration<double> iter(max_iter);
    try {
        lanczos.calculate_eigenvalues(iter,mygen);
        eigen = lanczos.eigenvalues();
    }
    catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }

    // Printing eigenvalues with error & multiplicities:
    PrintEigenvalues(lanczos);

    // call of eigenvectors function follows:
    std::cout << "\nEigen vectors computations for 3 lowest eigenvalues:\n\n";
    std::vector<double>::iterator start = eigen.begin();
    std::vector<double>::iterator end = eigen.begin()+3;
    std::vector<Vector> eigenvectors; // for storing the eigen vector.
    ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residual, status).

    try {
        lanczos.eigenvectors(start, end, std::back_inserter(eigenvectors),info,mygen);
    }
    catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }
    PrintEigenvectors(eigenvectors, info);

    return 0;
}

/***********************************************************
* Constructs spin-1/2 heisenberg hamiltonian.
***********************************************************/
void FillHermitianMatrix(Matrix &mat, int N){

    int k,n;
    int max = mat.size1();
    std::cout << max <<std::endl;

    std::bitset<std::numeric_limits<int>::digits> a;

    // initialize the matrix mat to zero
    for(int i=0; i<max; i++){
        for(int j=0; j<=i; j++){
            mat(i,j) = 0.0;
        }
    }
    // fill the matrix accodfing to spin-1/2 heisenberg chain
    for (int i=0; i<max; i++){
        for(int j=0; j<N; j++){

            k = (j+1) % N;
            a = std::bitset<std::numeric_limits<int>::digits>(i);
            if(a[j] == a[k]) mat(i,i) +=0.25;
            else{
                mat(i,i) -=0.25;
                a.flip(j); a.flip(k);
                n=static_cast<int>(a.to_ulong());
                mat(i,n) = 0.5;
            }
        }
    }
}
/***********************************************************
* Constructs test hamiltonian
***********************************************************/
void FillHermitianMatrix2(Matrix &mat)
{
    int i, j, n = 1;

    for(int i=0; i<mat.size1(); i++){
        for(int j=0; j<=i; j++){
            mat(i,j) = n++;
        }
    }
    std::cout << std::endl << "Printing matrix\n";
    std::cout << "--------------------------------\n\n";
    std::cout << mat << std::endl;
}
/***********************************************************
* Printing the eigenvalues
***********************************************************/
void PrintEigenvalues(ietl::lanczos<Matrix,Vecspace> &lanczos){

    std::vector<double> eigen;
    std::vector<double> err;
    std::vector<int> multiplicity;

    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();

    // Printing eigenvalues with error & multiplicities:
    std::cout << "#         eigenvalue         error         multiplicity\n";
    std::cout.precision(10);
    for (int i=0;i<eigen.size();++i)
        std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t"
            << multiplicity[i] << "\n";
}
/***********************************************************
* Printing the eigenvectors
***********************************************************/
void PrintEigenvectors(std::vector<Vector> &eigenvectors, ietl::Info<double> &info ){

    std::cout << "Printing eigenvectors:" << std::endl << std::endl;
    for(std::vector<Vector>::iterator it = eigenvectors.begin(); it != eigenvectors.end(); it++){
        std::copy((it)->begin(),(it)->end(),std::ostream_iterator<value_type>(std::cout,"\n"));
        std::cout << "\n\n";
    }
    std::cout << " Information about the eigenvectors computations:\n\n";
    for(int i = 0; i < info.size(); i++) {
        std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
            << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
                << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
                    << info.residual(i) << " error_info(" << i+1 << "): "
                        << info.error_info(i) << std::endl << std::endl;
    }
}
/************************************************************/
