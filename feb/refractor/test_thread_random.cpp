#include <omp.h>
#include <iostream>
#include <random>

using namespace std;


void parallel_random(){
        uniform_real_distribution<double> u(0, 1);
    
    
        #pragma omp parallel
    {
        
        time_t seconds = time(NULL);

        int myseed = omp_get_thread_num();
        mt19937 mt(myseed + seconds);
        
        
        #pragma omp for
        for(int i=0; i<12; i++){
            #pragma omp critical
            {
                cout << "thread num " << myseed << "  " << u(mt) << endl;
        }
        }


    }

    
}
// the same thread will use the same u(mt) again. 
int main(){
    
    
    parallel_random();
    
    
//    uniform_real_distribution<double> u(0, 1);
//
//    cout << "parallel thread unsafe" << endl;
//    #pragma omp parallel for
//    for(int i=0; i < 12; i++){
//        int threadNum = omp_get_thread_num();
//        time_t seconds = time(NULL);
//        mt19937 mt(seconds + (double)threadNum);
//
//#pragma omp critical
//        {
//            cout << "ThreadNum " << threadNum << " " << u(mt) << endl;
//        }
//
//    }
//
//
//    cout << endl << "sequential random generate" << endl;
//
//        time_t seconds = time(NULL);
//    for(int i=0; i < 12; i++){
//            int threadNum = omp_get_thread_num();
//            mt19937 mt(seconds);
//            cout << "ThreadNum " << threadNum << " " << u(mt) << endl;
//    }
//
//
    
    
// This is very subtle. BUT YOU WILL NEED TO INITIALIZE mt OUTSIDE OF FOR LOOP!!!
// the code below generates distinct random numbers while the code above does NOT.
// the reason is that u(mt) generates based on "step" so if you re-define mt again
// inside for loop then u(mt) will be step 1 and you'll get the same number over
// and over again!!
    
//    time_t seconds = time(NULL);
//    mt19937 mt(seconds);
//
//    cout << "seconds " << seconds << endl;
//    for (int i=0; i < 12; i++){
//        cout << u(mt) << endl;
//    }
    
    
}
