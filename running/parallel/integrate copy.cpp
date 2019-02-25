#include <iostream>
#include <omp.h>

using namespace std;

#define NUM_THREADS 4
static long num_steps = 20;
int parallel_integrate(){

    double sum[NUM_THREADS]; 
    int num_threads;
    // cout << "num_threads master " << num_threads << endl;

    double step = 1.0/(double)num_steps;

    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
//        cout << "Parallel num_threads " << num_threads << endl;
        int id = omp_get_thread_num();
        double x;
        int i; 
        sum[id] = 0.0;
        for (i=id; i <= num_steps; i = i + num_threads){
            x = i*step;
            double res = 4.0/(1.0+x*x);
            #pragma omp critical
            cout << "Parallel x " << x << " res " << res << endl;
//            cout << "id " << id << " sum " << sum[id] << " at x " << x << " with res " << res << endl;

            // cout << "id " << id << " | res| " << res << "| at x |" << x << endl;
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
    cout << "step: " << step << endl;
    cout << "parallel pi: " << pi << endl;
    return pi;
    }


int main(){
    double step = 1.0/(double)num_steps;
    double sum = 0.0;
    double x = 0.0;

    for (int i=0; i<=num_steps; i++){
        double res =4.0/(1.0+x*x);
        cout << "Serial x " << x << " is " << res << endl;
        x += step;
        sum += res;
    }

    double pi = sum*step;
    cout << " serial pi: " << pi << endl;


    parallel_integrate();
}

