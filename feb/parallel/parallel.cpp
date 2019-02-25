#include <omp.h>
#include <sys/time.h>
#include <iostream>
#include <random>

using namespace std;

double* parallel_random(int numRandoms, double from, double to) {
  double* randomArr = (double*)malloc(numRandoms * sizeof(double));
  int i;
  uniform_real_distribution<double> u(0, 1);

  int numThreads;

#pragma omp parallel
  {
    numThreads = omp_get_num_threads();

    // seconds since January 1, 1970
    time_t seconds = time(NULL);

    int threadNum = omp_get_thread_num();

    // seed = t + threadNum
    mt19937 mt(seconds + (double)threadNum);

#pragma omp for
    for (i = 0; i < numRandoms; i++) {
      randomArr[i] = from + to * u(mt);
    }
  }

  cout << "number of threads: " << numThreads << endl;
  return randomArr;
}

int main() {
  double start = omp_get_wtime();

  int numRandoms = 10000000;
  double from = 2.0;
  double to = 6.0;

  cout << "Num randoms " << numRandoms << " from " << from << " to " << to
       << endl;
  double* randomArr = parallel_random(numRandoms, from, to);

  //    for(int i = 0; i < numRandoms; ++i){
  //        cout << "arr " << i << " " << randomArr[i] << endl;
  //    }

  double end = omp_get_wtime();

  cout << "Process takes " << end - start << endl;
}
