#include <iostream>
#include <string>
#include "solve_copy.h"
#include "functions_copy.h"
#include "mpi.h"
//#include <sys/time.h>


/* Аргументы командной строки:  n - размерность матрицы
                                m - количество выводимых значений в матрице 
                                k - номер формулы для инициализиции матрицы
                                filename - имя файла (k==0), и отсутствует, если k!=0
                                */

int main(int argc, char *argv[])
{
    int my_rank;  // rank of current machine
    int total_p;  // total number of machines 

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &total_p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    if (argc < 4 || argc > 5)
    {
        if (my_rank == 0)
            std::cout << "incorrect number of arguments\n";
        MPI_Finalize();
        return -1;
    }

    int n, m, k;
    std::string filename = "";

    if (!sscanf(argv[1], "%d", &n))
    {
        if (my_rank == 0)
            std::cout << "invalid arguments\n";
        MPI_Finalize();
        return -1;
    }
    if (!sscanf(argv[2], "%d", &m))
    {
        if (my_rank == 0)
            std::cout << "invalid arguments\n";
        MPI_Finalize();
        return -1;
    }
    if (!sscanf(argv[3], "%d", &k))
    {
        if (my_rank == 0)
            std::cout << "invalid arguments\n";
        MPI_Finalize();
        return -1;
    }
    if(k<0 || k>5)
    {
        if (my_rank == 0)
            std::cout << "invalid arguments\n";
        MPI_Finalize();
        return -1;
    }
    if (argc == 5) 
    {
        if (k == 0)
            filename = argv[4];
        else
        {
            if (my_rank == 0)
                std::cout << "invalid arguments\n";
            MPI_Finalize();
            return -1;
        }
    }
    else
    {
        if (k == 0)
        {
            if (my_rank == 0)
                std::cout << "invalid arguments\n";
            MPI_Finalize();
            return -1;
        }
    }
    
    ////// arguments reading done

    //we need to store some part of n*(n+1) matrix, where its (n+1)-colomn is vector b
    // we divide matrix by colomns and store equal (+-1) number of them on each machine 
    int colomns;  //number of colomns for current machine 

    if (my_rank < (n+1) % total_p)              //  0   1   2   0   1   (example)   n == 4
    {                                           //  .   .   .   .   .               total_p == 3
        colomns = (n+1) / total_p +1;           //  .   .   .   .   .
    }                                           //  .   .   .   .   .
    else                                        //  .   .   .   .   .
    {                                           //              _ + _  - (n+1)%total_p
        colomns = (n+1) / total_p;              //
    }

    double *array = new (std::nothrow) double[n*colomns];
    if (!array)
    {
        std::cout << "Couldn't allocate memory for matrix\n";
        MPI_Finalize();
        return -1;
    }
    /*double *A = new (std::nothrow) double[n*colomns];
    if (!A)
    {
        std::cout << "Couldn't allocate memory for matrix copy\n";
        delete[] array;
        MPI_Finalize();
        return -1;
    }*/
    
    if (k == 0)
    {
        switch (readMatrix (n, array, filename, total_p, my_rank))
        {
            case -1:
                std::cout << "Memory allocation error\n";
                delete[] array;
                //delete[] A;
                MPI_Finalize();
                return -1;
            case -2:
                std::cout << "File cannot be opened, reading failed\n";
                delete[] array;
                //delete[] A;
                MPI_Finalize();
                return -2;
            case -3:
                std::cout << "File cannot be read / data is invalid\n";
                delete[] array;
                //delete[] A;
                MPI_Finalize();
                return -3;
        }

    }
    else
    {
        fillMatrix (n, array, k, total_p, my_rank);
    }

    //print (array, colomns, n, m);

    double *full_matrix = nullptr;  // full n*n - matrix 
    double *vector = nullptr;       // vector b 

    double matrix_norm;

    if (my_rank == 0)
    {
        full_matrix = new (std::nothrow) double[n*n];
        if (!full_matrix)
        {
            std::cout << "Couldn't allocate memory for full_matrix\n";
            delete[] array;
            MPI_Finalize();
            return -1;
        }
        vector = new (std::nothrow) double[n];
        if (!vector)
        {
            std::cout << "Couldn't allocate memory for vector\n";
            delete[] array, full_matrix;
            MPI_Finalize();
            return -1;
        }

        switch (fillFullMatrix (full_matrix, n, k, filename))
        {
            case -2:
                std::cout << "File cannot be opened, reading failed\n";
                delete[] array, full_matrix, vector;
                MPI_Finalize();
                return -2; 
            case -3:
                std::cout << "File cannot be read / data is invalid\n";
                delete[] array, full_matrix, vector;
                MPI_Finalize();
                return -3;
        }

        fillVector (vector, full_matrix, n);

        std::cout << " The initial matrix A:\n";
        print (full_matrix, n, n, m);

        std::cout << " The initial vector b\n";
        print (vector, 1, n, m);


        matrix_norm = norm (n, full_matrix);

        delete[] full_matrix;
    }


    // this vector we will share among all processes
    double *x = new (std::nothrow) double[n];  // x -- solution 
    if (!x)
    {
        std::cout << "Couldn't allocate memory for solution x\n";
        MPI_Finalize();
        return -1;
    }

    int err = 0;
    
    /*if (my_rank == 0)
    {
        matrix_norm = norm (n, full_matrix);
    }*/

    double *cosphi_copy = new (std::nothrow) double[n];
    double *sinphi_copy = new (std::nothrow) double[n];
    int *flag_copy = new (std::nothrow) int[n];

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&matrix_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    double time = MPI_Wtime();

    solve (n, array, x, &err, total_p, my_rank, matrix_norm, colomns, sinphi_copy, cosphi_copy, flag_copy);

    MPI_Barrier(MPI_COMM_WORLD);
    time = MPI_Wtime() - time;

    if (err == -4)
    {
        if (my_rank == 0)
            std::cout << "Matrix is degenerate\n";
        delete[] sinphi_copy, cosphi_copy;
        delete[] flag_copy;
        delete[] x;
        delete[] array;
        //delete[] full_matrix;
        delete[] vector;
        MPI_Finalize();
        return -4;
    }

    if (my_rank == 0)
    {

        full_matrix = new (std::nothrow) double[n*n];
        if (!full_matrix)
        {
            std::cout << "Couldn't allocate memory for full_matrix\n";
            delete[] x;
            delete[] flag_copy;
            delete[] sinphi_copy, cosphi_copy;
            delete[] array;
            delete[] vector;
            MPI_Finalize();
            return -1;
        }
        switch (fillFullMatrix (full_matrix, n, k, filename))
        {
            case -2:
                std::cout << "File cannot be opened, reading failed\n";
                delete[] array, full_matrix, vector, x;
                delete[] flag_copy;
                delete[] sinphi_copy, cosphi_copy;
                MPI_Finalize();
                return -2; 
            case -3:
                std::cout << "File cannot be read / data is invalid\n";
                delete[] array, full_matrix, vector, x;
                delete[] flag_copy;
                delete[] sinphi_copy, cosphi_copy;
                MPI_Finalize();
                return -3;
        }

        std::cout << "Total time: " << time << "\n";
        std::cout << " Solution (vector x): \n";
        print (x, 1, n, m);

        residual (n, full_matrix, vector, x);
        error (n, x);

        delete[] full_matrix;
    }


    delete[] sinphi_copy, cosphi_copy;
    delete[] flag_copy;
    delete[] array;
    delete[] vector;
    delete[] x;
    MPI_Finalize();
    return 0;
}
