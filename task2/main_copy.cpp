#include <iostream>
#include <string>
//#include <ctime>
#include "solve_copy.h"
#include "functions_copy.h"


/* Аргументы командной строки:  n - размерность матрицы
                                m - количество выводимых значений в матрице 
                                k - номер формулы для инициализиции матрицы
                                filename - имя файла (k==0), и отсутствует, если k!=0
                                p - число задач 
                                */

int main(int argc, char *argv[])
{
    if (argc < 5 || argc > 6)
    {
        std::cout << "incorrect number of arguments\n";
        return -1;
    }

    int n, m, k, p;
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
    if (!sscanf(argv[3], "%d", &k))
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
            filename = argv[4];
        else
        {
            std::cout << "invalid arguments\n";
            return -1;
        }
        if (!sscanf(argv[5], "%d", &p))
        {
            std::cout << "invalid arguments\n";
            return -1;
        }
    }
    else
    {
        if (!sscanf(argv[4], "%d", &p))
        {
            std::cout << "invalid arguments\n";
            return -1;
        }
    }

    /*if (p > 1000)
    {
        std::cout << "Too many threads";
        return -1;
    }*/

    ////// arguments reading done
    pthread_barrier_t barrier;
    pthread_barrier_init (&barrier, 0, p);


    //pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

    //pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    int err = 0;


    double *array = new (std::nothrow) double[n*n];
    if (!array)
    {
        std::cout << "Couldn't allocate memory for matrix\n";
        pthread_barrier_destroy (&barrier);
        return -1;
    }
    /*double *A = new (std::nothrow) double[n*n];
    if (!A)
    {
        std::cout << "Couldn't allocate memory for matrix copy\n";
        pthread_barrier_destroy (&barrier);
        delete[] array;
        return -1;
    }*/
    double *vector = new (std::nothrow) double[n];
    if (!vector)
    {
        std::cout << "Couldn't allocate memory for vector b\n";
        pthread_barrier_destroy (&barrier);
        delete[] array;
        return -1;
    }
    double *b = new (std::nothrow) double[n];
    if (!b)
    {
        std::cout << "Couldn't allocate memory for copy of vector b\n";
        pthread_barrier_destroy (&barrier);
        delete[] array;
        delete[] vector;
        return -1;
    }

    if (fillMatrix (array, n, k, filename) == -1)
    {
        pthread_barrier_destroy (&barrier);
        delete[] array;
        delete[] vector;
        delete[] b;
        return -2;
    }

    fillVector (vector, b, array, n);

    double *x = new (std::nothrow) double[n]; // x - solution 
    if (!x)
    {
        std::cout << "Couldn't allocate memory for solution x\n";
        pthread_barrier_destroy (&barrier);
        delete[] vector;
        delete[] array;
        delete[] b;
        return -1;
    }

    double *transp_array = new (std::nothrow) double[n*n];
    if (!transp_array)
    {
        std::cout << "Couldn't allocate memory for solution x\n";
        pthread_barrier_destroy (&barrier);
        delete[] vector, array, b, x;
        return -1;
    }
    int l=0;
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            transp_array[l] = array[i*n +j];
            l++;
        }
    }




    std::cout << " The initial matrix A:\n";
    print (array, n, n, m);

    std::cout << " Vector b:\n";
    print (vector, 1, n, m);

    pthread_t *threads = new (std::nothrow) pthread_t[p];
    if(!threads)
    {
        printf("Couldn't allocate memory for threads\n");
        pthread_barrier_destroy (&barrier);
        delete[] transp_array;
        delete[] vector;
        delete[] array;
        delete[] b;
        delete[] x;
        return -1;
    }

	Arg *args = new (std::nothrow) Arg[p];
	if(!args)
    {
        printf("Couldn't allocate memory for args\n");
        pthread_barrier_destroy (&barrier);
        delete[] transp_array;
        delete[] vector;
        delete[] array;
        delete[] b;
        delete[] x;
        delete[] threads;
        return -1;
    }


    ////////////////////
    ////////////////
    double *sinphi = new (std::nothrow) double[n];
    if (!sinphi)
    {
        std::cout << "Couldn't allocate memory for sinphi\n";
        pthread_barrier_destroy (&barrier);
        delete[] transp_array;
        delete[] vector;
        delete[] array;
        delete[] b;
        delete[] x;
        delete[] threads;
        delete[] args;
        return -1;
    }


    double *cosphi = new (std::nothrow) double[n];
    if (!cosphi)
    {
        std::cout << "Couldn't allocate memory for cosphi\n";
        pthread_barrier_destroy (&barrier);
        delete[] transp_array;
        delete[] vector;
        delete[] array;
        delete[] b;
        delete[] x;
        delete[] threads;
        delete[] args;
        delete[] sinphi;
        return -1;
    }


    int *flag = new (std::nothrow) int[n];
    if (!flag)
    {
        std::cout << "Couldn't allocate memory for flag\n";
        pthread_barrier_destroy (&barrier);
        delete[] transp_array;
        delete[] vector;
        delete[] array;
        delete[] b;
        delete[] x;
        delete[] threads;
        delete[] args;
        delete[] sinphi;
        delete[] cosphi;
        return -1;
    }
    for (int i = 0; i < n; i++)
    {
        flag[i] = 0;
    }
    ////////////////
    ////////////////////


    
    for (int i = 0; i < p; i++)
    {
        args[i].n = n;
        args[i].p = p;
        args[i].my_rank = i;
        args[i].err = &err;
        args[i].A = transp_array;
        args[i].b = vector;
        args[i].x = x;
        args[i].cosphi = cosphi;
        args[i].sinphi = sinphi;
        args[i].barrier = &barrier;
        //args[i].mutex = &mutex;
        args[i].flag = flag;
    }
    
    for(int i = 0; i < p; i++)
    {
        if (pthread_create(threads + i, 0, &solve, args+i))
        {
            printf ("Cannot create thread %d\n", i);
            /*for(int i=0; i<p; i++)
            {
                pthread_cancel(threads[i]);
                if (i%1000 == 0)
                    std::cout << i << " "; 
            }
            //pthread_barrier_destroy (&barrier);
            delete[] transp_array;
            delete[] vector;
            delete[] array;
            delete[] b;
            delete[] x;
            delete[] threads;
            delete[] args;
            delete[] sinphi;
            delete[] cosphi;
            delete[] flag;
            //pthread_barrier_destroy (&barrier);
            return -1;*/
        }  
        
    }
    
    for(int i=0; i<p; i++)
    {
        pthread_join(threads[i], 0);
    }
    
    pthread_barrier_destroy (&barrier); 

    if (err == -3)  // matrix is degenerate
    {
        std::cout << "matrix is degenerate\n";
        delete[] transp_array;
        delete[] vector;
        delete[] array;
        delete[] b;
        delete[] x;
        delete[] threads;
        delete[] args;
        delete[] sinphi;
        delete[] cosphi;
        delete[] flag;
        return -3;
    }

    
    std::cout << " Total time for solving: " << args[0].time_thr << "\n";
    
    std::cout << " Solution (vector x): \n";
    print (x, 1, n, m);

    residual (n, array, b, x);
    error (n, x);

    delete[] transp_array;
    delete[] array;
    delete[] vector;
    delete[] threads;
    delete[] args;
    delete[] b;
    delete[] x;
    delete[] sinphi;
    delete[] cosphi;
    delete[] flag;
    return 0;
}
