#include <iostream>
#include <cmath>
#include "functions_copy.h"


int solve (const int n, double *array, double *spector, const double eps, double *x, double *z)
{
    using std::abs;
    double epsilon = 1e-10;
    double matrix_norm = matrixNorm (n, array);

    for (int j=0; j<n-1; j++)
    {
        for (int i=1; i<n; i++)
        {
            if (abs(array[i*n +j] - array[j*n +i]) > epsilon*matrix_norm)
            {
                return -4;
            }
        }
    }

    // matrix is symmetric 
    // приводим матрицу к трехдиагональному виду унитарным подобием методом отражений 
    reflect (n, array, x, z);

    std::cout << "PASSED\n";

    print (array, n, n, 12);



    double b = matrixNorm (n, array);
    double a = -b;
    //all eighenvalues are within [a, b] 
    //std::cout << "PASSED   " << a << "   " << b <<"\n";

    while (n_(n, array, b) == -5 || n_(n, array, a) == -5)
    {
        a -= epsilon;
        b += epsilon;           
    }

    //all eighenvalues are within (a, b) 
    
    int index = 0;
    recur (n, array, eps, spector, &index, a, b);

    //n_(n, array, 0);
    //std::cout << n_(n, array, a) << "\n";
    //std::cout << n_(n, array, b) << "\n";
    //n_(n, array, 10);


    return 0;
}



