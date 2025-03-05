#include <iostream>
#include <math.h>
#include <complex.h>

const double PI = 3.14159265358979323846;


int main() {
    int N = 4;
    double x[N] = { 1,2,3,4 };
    std::complex<double> X[N] = {0};
    
    for (int k = 0; k < N; k++)
    {
        for(int n = 0; n < N; n++)
            {
                double angle = (-2*PI*k*n)/N ; 

                //need a function ( math fun or something to differentiate real and imaginary values )
                X[k] += x[n] * ( std::complex<double>( cos(angle), sin(angle) ) );
            }
    }

    for (int k = 0; k < N; k++) {
        std::cout << "X(" << k << ") = " << X[k] << std::endl;
    }
    
    return 0;
}


// // example to illustrate the use of sin(), cos() and tan()
// #include <iostream>     
 
// // CPP program to illustrate
// // std::complex, std::cos, std::sin, std::tan
// #include <complex> 
// using namespace std;
 
// // driver program
// int main ()
// {    
//   // initializing the complex: (-1.0+0.0i)
//   complex<double> mycomplex (0.0, 1.0);
 
//   // use of cos()
//   cout << "The cos of " << mycomplex << " is "
//        << cos(mycomplex) <<endl;
       
//   // use of sin()
//   cout << "The sin of " << mycomplex << " is "
//        << sin(mycomplex) <<endl; 
       
//   // use of tan()
//   cout << "The tan of " << mycomplex << " is "
//        << tan(mycomplex) <<endl; 
 
//   return 0;
// }
