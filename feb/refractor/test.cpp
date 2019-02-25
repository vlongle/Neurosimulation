#include "Vector.h"
#include "simulator.h"
#include <iostream>
#include <omp.h>


using namespace std;


void test_vector(){
    cout << "From test.cpp: test vector" << endl;
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> u(0, 1);
    
    Vector<double> test;
    test.reserve(10);
    
    struct RandomEngine randEng = {mt, u};
    
    
    
    for (int i=0; i < 5; i++){
        test.put(i, i+10);
        cout << "test[" << i << "]: " << test[i] << endl;
    }
}

void test_simulator(int rank){
    #pragma omp critical
    cout << "test.cpp: test simulator from thread " <<  rank << endl;
    
    struct Parameters p = {10.0, -5.0, 7.0, 1.0, 69.5, 69.0};
    
//    random_device rd;
//    mt19937 mt(rd());
    mt19937 mt(rank);
    uniform_real_distribution<double> u(0, 1);
    
    
    struct RandomEngine randEng = {mt, u};
    

    
    
    Population pop(&p, &randEng);
    
    struct Result res = pop.result;
//    cout << "population result " << res.E_spike << " , " << res.I_spike << endl;
//    cout << "population result " << pop.parameters->SEE << " , " << pop.parameters->SEI << endl;
//
//    cout << "population result " << pop.parameters->SIE << " , " << pop.parameters->SII << endl;
//
//
//    cout << "population result " << pop.parameters->kickE << " , " << pop.parameters->kickI << endl;
//
    
    p.kickE = -100.69;
    
//    cout << "population result " << pop.parameters->kickE << " , " << pop.parameters->kickI << endl;
    
    
    pop.simulate(10);
    
}

int main(){
    cout << "main of test.cpp" << endl;
    #pragma omp parallel for
    for(int i=0; i < 2; i++){
        int rank = omp_get_thread_num();
        test_simulator(rank);
    }
    cout << "Number of thread: " << omp_get_num_threads() << endl;
}
