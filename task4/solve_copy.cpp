#include <iostream>
#include <cmath>
#include "functions_copy.h"


int solve (const int n, double *array, double *spector, const double eps)
{
    using std::abs;
    double epsilon = 1e-16;
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
    reflect (n, array);

    std::cout << "PASSED\n";

    print (array, n, n, 6);



    double b = matrixNorm (n, array);
    double a = -b;
    //all eighenvalues are within [a, b] 

    int index = 0;
    //recur (n, array, eps, spector, &index, a, b);

    //n_(n, array, 0);
    //std::cout << n_(n, array, a) << "\n";
    //std::cout << n_(n, array, b) << "\n";
    //n_(n, array, 10);

    /*если запускать на единичной матрице, то b=1 и матрица A-bE = 0, и теорема становится неприменимой
    по теореме alpha нельзя брать равным одному из собственных значений, а тут норма матрицы равна собственному 
    значению     что тогда делать?*/


    return 0;
}



