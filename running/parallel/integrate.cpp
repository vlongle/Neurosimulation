// https://www.youtube.com/watch?v=OuzYICZUthM
#include <iostream>
#include <omp.h>

using namespace std;

#define NUM_THREADS 4
static long num_steps = 20;


/*
 There are 2 problems with this code
 1. for (i = 0, ...) should have been for(i = id, ...)
 2. num_threads = omp_get_num_threads(); should have been put inside #pragma omp parallel
 --> At the beginning, there's only 1 threads. omp_get_num_threads() only returns the correct
 result once we are actually in the parallel region. However, the actual constant num_threads
 should be put outside #pragma so that we can access it later. 
 
 */
int parallel_integrate(){

    omp_set_num_threads(NUM_THREADS);


    double sum[NUM_THREADS]; 
    int num_threads = omp_get_num_threads();
    cout << "num threads outside " << num_threads << endl;
    double step = 1.0/(double)num_steps;
    #pragma omp parallel
    {
//        num_threads = omp_get_num_threads();
//        #pragma omp critical
//        cout << "num threads " << num_threads << endl;
        
        int id = omp_get_thread_num();
        double x;
        int i; 
        for (i=0, sum[id] = 0.0; i <= num_steps; i = i + num_threads){
            x = i*step;
            double res = 4.0/(1.0+x*x);
                #pragma omp critical
             cout << "id " << id << " | res| " << res << "| at x |" << x << endl;
            // cout << "id" << id << " " << sum[id] << endl;
            sum[id] += res; 

        }
        }
    double totalSum = 0.0;
    // sum the sum array
    for (int i=0; i < NUM_THREADS; i++){
        cout << "sum[" << i << "] is " << sum[i] << endl;
        totalSum += sum[i];
    }

    double pi = totalSum*step;
    cout << "pi: " << pi << endl;
    return pi;
    }


int main(){
    double step = 1.0/(double)num_steps;
    double sum = 0.0;
    double x = 0.0;

    for (int i=0; i<num_steps; i++){
        x += step;
        sum += 4.0/(1.0+x*x);
    }

    double pi = sum*step;
    cout << "pi: " << pi << endl;


    parallel_integrate();
}

