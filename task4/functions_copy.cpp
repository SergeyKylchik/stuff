#include <iostream>
#include <string>
#include <cmath>
#include "functions_copy.h"


double euclidNorm (const int n, double *x) // euclid norm in R^n
{
    double value = 0;
    for (int i=0; i<n; i++)
    {
        value += x[i]*x[i];
    }
    return sqrt(value);
}

double matrixNorm (const int n, double *array)
{
    using std::abs;
    double value = 0;
    double max = 0;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            value += abs(array[i*n +j]);
        }
        if (value > max)   // ???
            max = value;
        value = 0;
    }
    return max;
}

double f(int k, int n, int i, int j)  // i,j = 1, ..., n (NOT from 0!)
{
    switch (k)
    {
        case 1:
            return n - ((i<j)?j:i) + 1;
        case 2:
            if (i == j)
                return 2;
            if (i - j == 1 || j - i == 1)
                return -1;
            return 0;
        case 3:
            if (i == j && i<n)
                return 1;
            if (i == n)
                return j;
            if (j == n)
                return i;
            return 0;
        case 4:
            return 1./(i+j-1);
    }
    return 0;
}

int fillMatrix (double *array, double *array_copy, int n, int k, std::string filename)
{
    if (k == 0)
    {
        FILE *fp = fopen (filename.c_str(), "r");
        if (!fp)
            {
                return -2;
            }
        int i;
        for (i=0; i < n*n; i++)
        {
            if (!fscanf (fp, "%lf", array + i))
            {
                return -3;
            }
            array_copy[i] = array[i];
        }
        fclose(fp);
    }
    else
    {
        for (int i=0; i < n*n; i++)
        {
            array[i] = f(k,n,(i/n)+1,(i%n)+1);
            array_copy[i] = array[i];
        }
    }
    return 0;
}



int reflect (int n, double *array, double *x, double *z)
{

    //double *x = new (std::nothrow) double[n];
    //double *z = new (std::nothrow) double[n];

    for (int index=0; index < n-2; index++)
    {
        double epsilon = 1e-16;
        // сначала вычисляем вектор х, по которому строится матрица отражений 
        // 
        // array [(index+1)*n +index]
        // ..
        // ..
        // ..
        // ..
        // array [(n-1)*n +index]
        double s = 0;
        double norm_a, norm_x;
        for (int i = index+2; i<n; i++)
        {
            s += array[i*n + index]*array[i*n + index];
        }
        if (std::abs(s) <= epsilon)
        {
            if (array[(index+1)*n +index] < 0)
            {
                array[(index+1)*n +index] = - array[(index+1)*n +index];
                array[index*n + index+1] = - array[index*n + index+1];
            }
            continue;
        }
        norm_a = sqrt ((array[(index+1)*n +index] * array[(index+1)*n +index]) + s);

        x[index+1] = array[(index+1)*n +index] - norm_a;
        for (int i = index+2; i<n; i++)
        {
            x[i] = array[i*n + index];
        }
        norm_x = sqrt (x[index+1]*x[index+1] + s);
        /*if (norm_x == 0)
            continue;*/
        for (int i = index+1; i<n; i++)
        {
            x[i] /= norm_x;
        }


        /*if (index < 4)
        {
            std::cout << "IN REFLECT  " << index << "  " << s << "  " << norm_a << "  " << norm_x << "\n";
        }*/

        // остается перемножить U*A*U, где U=U(x), а A:
        // array[(index+1)*n + index+1]    ...   array[(index+1)*n + n-1]
        //   .                                          .
        //   .                                          .
        //   .                                          .
        // array[(n-1)*n + index+1]        ...   array[(n-1)*n + n-1]

        double dot_pr = 0;
        double temp;
        for (int i = index+1; i<n; i++)
        {
            temp = 0;
            for (int j = index+1; j<n; j++)
            {
                temp += array[i*n + j]*x[j];
            }
            z[i] = temp;
        }
        for (int i = index+1; i<n; i++)
        {
            dot_pr += x[i]*z[i];
        }
        for (int i = index+1; i<n; i++)
        {
            z[i] = 2*(z[i] - dot_pr*x[i]);
        }
        for (int i = index+1; i<n; i++)
        {
            for (int j = index+1; j <= i; j++)
            {
                array[i*n + j] = array[i*n + j] - z[i]*x[j] - z[j]*x[i];
                array[j*n + i] = array[i*n + j];
            }
        }
        
        array[(index+1)*n + index] = norm_a;
        array[index*n + index+1] = norm_a;
        for (int i = index+2; i<n; i++)
        {
            array[i*n + index] = 0;
            array[index*n + i] = 0;
        }

    }

    //delete[] x;
    //delete[] z;
    return 0;
}


int n_(int n, double *array, double alpha)
{
    double epsilon = 1e-16;
    //std::cout << "PASSED alpha:" << alpha << "\n";
    int count = 0;
    double l1, l2; 
    l1 = array[0] - alpha;
    //std::cout << "PASSED l1:" << l1 << "\n";
    if (l1 < 0)
        count++;
    for (int k=1; k<n; k++)
    {
        if (std::abs(l1) < epsilon)
            return -5;
        l2 = (array[k*n + k] - alpha) - array[k*n + k-1] * (array[(k-1)*n + k]/l1);
        if (l2 < 0)
            count++;
        l1 = l2;
        //std::cout << "PASSED l1:" << l1 << "\n";
    }

    return count;
}


/*
void recur (int n, double *array, double eps, double *spector, int *index, double a, double b)
{
    double epsilon = 1e-16;
    
    if ( n_(n, array, b)-n_(n, array, a) > 0 )
    {
        if (b-a < eps)
        {
            std::cout << "HERE\n";
            for (int i=0; i < (n_(n, array, b)-n_(n, array, a)); i++)
            {
                //std::cout << "PASSED";
                //fflush(stdout);
                //std::cout << "PASSED: " << *index << "\n" << "a + b = " << (a+b)/2 << " a = " << a << " b = "<< b << "\n";
                //std::cout << (n_(n, array, b)-n_(n, array, a)) << "\n";
                
                spector[*index] = (a+b)/2;
                (*index)++;
                std::cout << *index << "\n";
            }
        }
        else
        {
            double p = (a+b)/2;
            while (n_(n, array, p) == -5)
            {
                p += epsilon;//*(b-a);
                //std::cout <<  "HERE\n";
            }
            //std::cout << "FUCK1\n";
            recur (n, array, eps, spector, index, a, p);
            std::cout << "FUCK2\n";
            recur (n, array, eps, spector, index, p, b);
            //recur (n, array, eps, spector, index, a, p);
        }
    }
}*/


void recur (int n, double *array, double eps, double *spector, int *index, double a, double b)
{
    double epsilon = 1e-14;
    //std::cout << "PASSED: " << *index << "\n" << " a = " << a << " b = "<< b << "\n";
    //std::cout << (n_(n, array, b)-n_(n, array, a)) << "\n";
    if ( n_(n, array, b)-n_(n, array, a) > 0 )
    {
        std::cout.precision(30);
        std::cout << "PASSED: " << *index << "\n" << " a = " << a << " b = "<< b << "\n";
        //std::cout << "PASSED\n";
        if (b-a > eps)
        {
            //std::cout << "PASSED\n";
            double p = (a+b)/2;
            std::cout << "PASSED:   " << p << "\n";
            while (n_(n, array, p) == -5)
            {
                p += epsilon;
            }
            //std::cout << "PASSED:   " << p << "\n";
            recur (n, array, eps, spector, index, a, p);
            recur (n, array, eps, spector, index, p, b);
        }
        else
        {
            for (int i=0; i < (n_(n, array, b)-n_(n, array, a)); i++)
            {
                spector[*index] = (a+b)/2;
                (*index)++;
                //std::cout << *index << "\n";
                //std::cout.precision(30);
                //std::cout << "PASSED: " << *index << "\n(a + b)/2 = " << (a+b)/2. << " a = " << a << " b = "<< b << "\n";
                //double p = (a+b);
                //std::cout << p << "\n";
                //std::cout << p/2. << "\n";
            }
        }
    }


}








void print (double *array, int l, int n, int m)
{
    for (int i=0; i < ((l<m)?l:m); i++)
    {
        for (int j=0; j < ((n<m)?n:m); j++)
        {
            printf(" %10.3e", array[i*n + j]); 
        }
        printf("\n");
    }
    printf("\n");
}


double error1 (const int n, double *array_copy, double *spector)
{
    double trace = 0, sum = 0;
    for (int i=0; i<n; i++)
    {
        trace += array_copy[i*n +i];
        sum += spector[i];
    }
    return (trace<sum)?(sum-trace):(trace-sum);
}

double error2 (const int n, double *array_copy, double *spector)
{
    double sum = 0;
    for (int i=0; i<n; i++)
    {
        sum += spector[i]*spector[i];
    }
    return std::abs(euclidNorm(n*n, array_copy) - sqrt(sum));
}

