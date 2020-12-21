#include <iostream>
#include <cmath>
#include "functions_copy.h"

void* solve (void *pa) 
{
    using std::abs;


    
    Arg *a = (Arg*) pa;
    int n = a -> n;
    int p = a -> p;
    int my_rank = a -> my_rank;
	double *transp_array = a -> A;
	double *vector = a -> b;
    double *x = a -> x;
    //double *sinphi = a -> sinphi;
    //double *cosphi = a -> cosphi;
    double *sinphi = new (std::nothrow) double[n];
    double *cosphi = new (std::nothrow) double[n];


    int *flag = a -> flag;
    

    ///////
    pthread_barrier_wait (a -> barrier);
    double time;
    time = get_full_time();
    ///////

    const double epsilon = 1e-16;
    double eps = epsilon * norm (n, transp_array);
    double temp, norm;

    //std::cout << "thread # " << my_rank << " PASSED\n" << eps << "\n";


    int i, j, k;
    for (int index = 0; index < n-1; index++) // (n - index) - current size of submatrix 
    {
        
        if (*(a -> err) == 0)
        {
            if (my_rank == index%p)
            {
                for (i = index+1; i < n; i++)
                {
                    if (abs(transp_array[index*(n+1)]) < eps && abs(transp_array[i*n + index]) < eps)
                    {   // skip the rotation if the rotating vector is (0,0)
                        if (i == n-1)  // case in which the entire colomn is 0 
                        {
                            *(a -> err) = -3;
                        }
                        a->cosphi[i] = 0;
                        a->sinphi[i] = 0;
                        flag[i] = 1;
                    }
                    else
                    {
                        norm = sqrt (transp_array[index*(n+1)]*transp_array[index*(n+1)] + transp_array[index*n + i]*transp_array[index*n + i]);
                        cosphi[i] = a->cosphi[i] = transp_array[index*(n+1)]/norm;
                        sinphi[i] = a->sinphi[i] = - transp_array[index*n + i]/norm;
                        

                        // rotation of the "index" colomn (index - is its number)
                        temp = transp_array[index*n + index]*a->cosphi[i] - transp_array[index*n + i]*a->sinphi[i];
                        transp_array[index*n + i] = transp_array[index*n + index]*a->sinphi[i] + transp_array[index*n + i]*a->cosphi[i];
                        transp_array[index*n + index] = temp;
                    }

                }
            }
            
            pthread_barrier_wait (a -> barrier);  // barrier


            
            for (k = index+1; k<n; k++)
            {
                sinphi[k] = a->sinphi[k];
                cosphi[k] = a->cosphi[k];
            }

            if (*(a -> err) == 0)
            {
                
                for (j = index + 1 + my_rank; j < n; j += p)  // rotation of the colomns
                //for (int j = n-1-my_rank; j > index; j -= p)  
                {
                    for (i = index+1; i < n; i++)
                    {
                        if (flag[i] == 0) // cosphi and sinphi are well-defined 
                        {
                            temp = transp_array[j*n + index]*cosphi[i] - transp_array[j*n + i]*sinphi[i];
                            transp_array[j*n + i] = transp_array[j*n + index]*sinphi[i] + transp_array[j*n + i]*cosphi[i];
                            transp_array[j*n + index] = temp;
                        }
                    }
                }
                if (my_rank == n%p)
                {
                    for (i = index+1; i < n; i++)
                    {
                        if (flag[i] == 0) // cosphi and sinphi are well-defined 
                        {
                            temp = vector[index]*cosphi[i] - vector[i]*sinphi[i];
                            vector[i] = vector[index]*sinphi[i] + vector[i]*cosphi[i];
                            vector[index] = temp;
                        }
                    }
                }
            }
            
            pthread_barrier_wait (a -> barrier);  // barrier
        }
    }


    if (my_rank == 0)
    {
        if (abs(transp_array[n*n-1]) < eps) // Matrix is degenerate
        {
            *(a -> err) = -3;
        }
    }
    pthread_barrier_wait (a -> barrier);  // barrier

        // inverse gaussian elimination
    //print (array, n, n, 4);

    if (*(a -> err) == 0)
    {
        for (j = n-1; j>=0; j--)
        {
            
            if (my_rank == 0)
            {
                vector[j] /= transp_array[j*n +j];
                //array[j*n +j] = 1;
                x[j] = vector[j];
            }

            pthread_barrier_wait (a -> barrier);  // barrier

            for (i = j-1-my_rank; i>=0; i -= p)
            {
                vector[i] -= transp_array[j*n +i]*vector[j]; 
                //array[i*n +j] = 0;
            }
            
            pthread_barrier_wait (a -> barrier);  // barrier
        }

    }
    

    delete[] sinphi;
    delete[] cosphi;


    pthread_barrier_wait (a -> barrier);
    a -> time_thr = get_full_time() - time;
    return 0;
}
