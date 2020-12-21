#include <iostream>
#include <cmath>
#include "functions_copy.h"
#include "mpi.h"

void solve (int n, double *array, double *x, int *err, int total_p, int my_rank, int matrix_norm, int colomns)
{
    using std::abs;

    const double epsilon = 1e-16;
    double eps = epsilon * matrix_norm;
    double temp, norm;


    double *sinphi = new (std::nothrow) double[n];
    double *cosphi = new (std::nothrow) double[n];
    int *flag = new (std::nothrow) int[n];
    for (int i; i<n; i++)
    {
        flag[i] = 0;
    }
    sinphi[0] = 0;
    cosphi[0] = 0;


    for (int index = 0; index < n-1; index++)
    {
        if (*err == 0)
        {
            if (my_rank == index%total_p)  // "index" colomn belongs to my_rank process
            {
                // cosphi and sinphi evaluation 

                // array[(index/total_p)*n]  --  1st element of the "index" colomn
                int begin = (index/total_p)*n;
                // array[begin + 0]  --  1st element of the "index" colomn

                for (int i = index+1; i < n; i++)
                {
                    if (abs(array[begin + index]) < eps && abs(array[begin + i]) < eps)
                    {  // skip the rotation if the rotating vector is (0,0)

                        if (i == n-1)  // case in which the entire colomn is 0
                        {
                            //std::cout << "PASSED 1\n";
                            //fflush(stdout);
                            *err = -4;
                        }
                        cosphi[i] = 0;
                        sinphi[i] = 0;
                        flag[i] = 1;
                    }
                    else
                    {
                        norm = sqrt (array[begin + index]*array[begin + index] + array[begin + i]*array[begin + i]);
                        cosphi[i] = array[begin + index]/norm;
                        sinphi[i] = - array[begin + i]/norm;

                        // rotation of the "index" colomn (index - is its number)
                        /*
                        temp = array[begin + index]*cosphi[i] - array[begin + i]*sinphi[i];
                        array[begin + i] = array[begin + index]*sinphi[i] + array[begin + i]*cosphi[i];
                        array[begin + index] = temp;
                        */
                        array[begin + index] = array[begin + index]*cosphi[i] - array[begin + i]*sinphi[i];


                    }

                }

            }

            /*if (my_rank == (n-1)%total_p)
            {
                std::cout << array[n*colomns-1] << '\n';
                fflush(stdout);
            }*/




            //MPI_Barrier(MPI_COMM_WORLD);

            // Broadcasting 
            MPI_Bcast(sinphi, n, MPI_DOUBLE, index%total_p, MPI_COMM_WORLD);
            MPI_Bcast(cosphi, n, MPI_DOUBLE, index%total_p, MPI_COMM_WORLD);
            MPI_Bcast(flag, n, MPI_INT, index%total_p, MPI_COMM_WORLD);
            MPI_Bcast(err, 1, MPI_INT, index%total_p, MPI_COMM_WORLD);
            // 

            //MPI_Barrier(MPI_COMM_WORLD);

            if (*err == 0)
            {
                if (my_rank == index%total_p)
                {
                    for (int j = 0; j < colomns; j++)
                    {
                        if (j*total_p + my_rank != index)
                        {
                            for (int i = index+1; i < n; i++)               
                            {
                                if (flag[i] == 0) // cosphi and sinphi are well-defined 
                                {
                                    temp = array[j*n + index]*cosphi[i] - array[j*n + i]*sinphi[i];
                                    array[j*n + i] = array[j*n + index]*sinphi[i] + array[j*n + i]*cosphi[i];
                                    array[j*n + index] = temp;
                                }
                            }

                        }
                    }

                }
                else
                {
                    for (int j = 0; j < colomns; j++)
                    {
                        for (int i = index+1; i < n; i++)               
                        {
                            if (flag[i] == 0) // cosphi and sinphi are well-defined 
                            {
                                temp = array[j*n + index]*cosphi[i] - array[j*n + i]*sinphi[i];
                                array[j*n + i] = array[j*n + index]*sinphi[i] + array[j*n + i]*cosphi[i];
                                array[j*n + index] = temp;
                            }
                        }

                    }

                }


            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    if (my_rank == (n-1)%total_p)
    {
        if (total_p == 1)
        {
            if (abs(array[n*n-1]) < eps)
                *err = -4;
        }
        else
        {
            if (abs(array[n*colomns-1]) < eps) // Matrix is degenerate
            {
                //std::cout << "PASSED 2\n";
                //fflush(stdout);
                *err = -4;
            }
        }

    }

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(err, 1, MPI_INT, (n-1)%total_p, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);

    //double *var = new double[n];  // temporary 


    if (*err == 0)
    {
        // inverse gaussian elimination
        if (my_rank == n%total_p)
        {
            for (int i=0; i<n; i++)
            {
                x[i] = array[(n/total_p)*n +i];
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(x, n, MPI_DOUBLE, n%total_p, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);

        for (int j = n-1; j>=0; j--)
        {
            if (my_rank == j%total_p)
            {
                x[j] /= array[(j/total_p)*n +j];
                for (int i = j-1; i>=0; i--)
                {
                    x[i] -= array[(j/total_p)*n +i]*x[j];
                }

            }

            //MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(x, n, MPI_DOUBLE, j%total_p, MPI_COMM_WORLD);
            //MPI_Barrier(MPI_COMM_WORLD);
        }


    }

    delete[] sinphi;
    delete[] cosphi;
    delete[] flag;
}
