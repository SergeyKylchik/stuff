#include <iostream>
#include <string>
#include <ctime>
#include "solve_copy.h"
#include "functions_copy.h"


/* Аргументы командной строки:  n - размерность матрицы
                                m - количество выводимых значений в матрице 
                                eps - точночть нахождения собственных значений матрицы 
                                k - номер формулы для инициализиции матрицы
                                filename - имя файла (k==0), и отсутствует, если k!=0
                                */


int main(int argc, char *argv[])
{
    if (argc < 5 || argc > 6)
    {
        std::cout << "incorrect number of arguments\n";
        return -1;
    }

    int n, m, k;
    double eps;
    std::string filename = "";

    if (!sscanf(argv[1], "%d", &n))
    {
        std::cout << "invalid arguments\n";
        return -1;
    }
    if (!sscanf(argv[2], "%d", &m))
    {
        std::cout << "invalid arguments\n";
        return -1;
    }
    if (!sscanf(argv[3], "%lf", &eps))
    {
        std::cout << "invalid arguments\n";
        return -1;
    }
    if (!sscanf(argv[4], "%d", &k))
    {
        std::cout << "invalid arguments\n";
        return -1;
    }
    if(k<0 || k>4)
    {
        std::cout << "invalid arguments\n";
        return -1;
    }
    if (argc == 6) 
    {
        if (k == 0)
            filename = argv[5];
        else
        {
            std::cout << "invalid arguments\n";
            return -1;
        }
    }
    ////// arguments reading done

    double *array = new (std::nothrow) double[n*n];
    if (!array)
    {
        std::cout << "Couldn't allocate memory for matrix\n";
        return -1;
    }

    double *spector = new (std::nothrow) double[n];
    if (!spector)
    {
        std::cout << "Couldn't allocate memory for spector\n";
        delete[] array;
        return -1;
    }

    double *array_copy = new (std::nothrow) double[n*n];
    if (!array_copy)
    {
        std::cout << "Couldn't allocate memory for array_copy\n";
        delete[] array, spector;
        return -1;
    }

    /////  memory allocation done

    switch (fillMatrix(array, array_copy, n, k, filename))
    {
        case -2:
            std::cout << "File cannot be opened\n";
            delete[] array, spector, array_copy;
            return -2;
        case -3:
            std::cout << "File cannot be read / data is invalid\n";
            delete[] array, spector, array_copy;
            return -3;
    }

    std::cout << "The initial matrix:\n";
    print (array, n, n, m);


    double time;
    time = std::clock ();

    if (solve (n, array, spector, eps) == -4)
    {
        std::cout << "Matrix is not symmetric\n";
        delete[] array, spector, array_copy;
        return -4;
    }

    time = (std::clock () - time) / (double) CLOCKS_PER_SEC;
    std::cout << "Total time for solving: " << time << "\n";

    //std::cout << "Matrix spector:\n";
    //print (spector, 1, n, m);

    //std::cout << "Error #1: " << error1 (n, array_copy, spector) << "\n";
    //std::cout << "Error #2: " << error2 (n, array_copy, spector) << "\n";

    delete[] array;
    delete[] array_copy;
    delete[] spector;

    return 0;
}