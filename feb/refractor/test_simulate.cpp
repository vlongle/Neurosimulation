#include <omp.h>
#include "Vector.h"
#include "simulator.h"

#include <chrono>
#include <ctime>
#include <fstream>


#define NUM_SAMPLES 1

int NE = 300;
int NI = 100;


using namespace std;




void generate_random_connections(int numRandoms, Vector<Parameters> *pi, double SEE_low, double SEE_high,
                                 double SIE_low, double SIE_high, double SEI_low, double SEI_high,
                                 double SII_low, double SII_high){
    uniform_real_distribution<double> u(0, 1);
    
    int numThreads;
    
    // threads are re-using the same random number ==> very problematic.
//    #pragma omp parallel for
//    for(int i=0; i < 3*NUM_SAMPLES; i++){
#pragma omp parallel
    {
    
        time_t seconds = time(NULL);
        int thread_num = omp_get_thread_num();
        numThreads = omp_get_num_threads();
        
        
        mt19937 mt(seconds + (double)thread_num);


    #pragma omp for
    for(int i=0; i < 3*NUM_SAMPLES; i++){
        

        
//        #pragma omp critical
//        cout << "thread " << threadNum << endl;
        
        

        double kickE, kickI;

        if (i % 3 == 0){
            kickI = 1000.0;
            kickE = 1000.0;
        }

        else if (i % 3 == 1) {
            kickI = 3000.0;
            kickE = 3000.0;
        }

        else{
            kickI = 5000.0;
            kickE = 5000.0;
        }

        
        
        
        double SEE = SEE_low + SEE_high * u(mt);
        double SIE = SIE_low + SIE_high * u(mt);
        double SEI = SEI_low + SEI_high * u(mt);
        double SII = SII_low + SII_high * u(mt);

        
//#pragma omp critical
//        {
//        cout << "thread " << threadNum << " " << kickE << " " << SEE << endl;
////            cout << u(mt) << endl;
//            }
        
        (*pi)[i].kickE = kickE;
        (*pi)[i].kickI = kickI;

        (*pi)[i].SEE = SEE;
        (*pi)[i].SIE = SIE;
        (*pi)[i].SEI = SEI;
        (*pi)[i].SII = SII;

        
        (*pi)[i].NE = NE;
        (*pi)[i].NI = NI;
        
        
    }

}
    
        cout << "Number of threads: " << numThreads << endl;

}



void parallel_simulate(Vector<Population> *populations, const double terminate_time){
    
    #pragma omp parallel for
    for(int i=0; i < 3*NUM_SAMPLES; i++){
        int thread_num = omp_get_thread_num();
        (*populations)[i].set_random_engine(thread_num);
        
//        (*populations)[i].simulate(terminate_time);
        (*populations)[i].test_simulate(terminate_time);

        
    }
}

void save_to_disk(const char* file_name, Vector<Population> *pops){
    ofstream myfile;
    // open in append mode
    myfile.open(file_name, ofstream::app);
    
    auto end = chrono::system_clock::now();

    time_t end_time = chrono::system_clock::to_time_t(end);
    myfile << endl << endl << "Data generated at " << ctime(&end_time) << endl << endl;
    
    for (int i=0; i < 3*NUM_SAMPLES; i++){
        Parameters* pi = (*pops)[i].parameters;
        Result res    = (*pops)[i].result;
        
        myfile << "SEE, SEI, SIE, SII " << pi->SEE << " " << pi->SEI << " " << pi->SIE << " "
        << pi->SII << " kickE, kickI " << pi->kickE << " " << pi->kickI << " | E_spike, I_spike " <<
        res.E_spike << " " << res.I_spike << endl;
    }
    
    myfile.close();
}
int main(){

    Vector<Parameters> pi;
    pi.reserve(3*NUM_SAMPLES);


    double SEE_low = 2.0;
    double SEE_high = 6.0;

    double SIE_low = 2.0;
    double SIE_high = 6.0;

    double SEI_low = -5.0;
    double SEI_high = -1.0;

    double SII_low = -5.0;
    double SII_high = -1.0;

    
    
    
    
    
    
    // randomly generate parameters pi in parallel
    generate_random_connections(3*NUM_SAMPLES, &pi, SEE_low, SEE_high, SIE_low,
                                SIE_high, SEI_low, SEI_high, SII_low, SII_high);



    cout << "Test_simulate: test done generate_random_connections" << endl;

    
//        Vector<fakeClass> pops;

    Vector<Population> pops;
    pops.reserve(3*NUM_SAMPLES);
    
    
    

    
    // initialize each population with its parameters
    for(int i=0; i<3*NUM_SAMPLES; i++){
        pops[i].set_parameters(&pi[i]);
    }
    
    
    // simulation starts in here !
    const double terminate_time = 5;
    parallel_simulate(&pops, terminate_time);
    
    
    
    save_to_disk("training_data.txt", &pops);



    
    
    
//
//
//
//
//
//
//        for (int i=0; i < 3*NUM_SAMPLES; i++){
//            cout << "result " << pops[i].result.E_spike << " " << pops[i].result.I_spike << endl;
//            cout << "population[" << i << "]" << endl;
//
//            Population p = pops[i];
//            Parameters* par = p.parameters;
//
//
//            cout << par->SEE << " " << par->SIE << " " << par->SEI << " " << par->SII << endl;
//            cout << par->kickI << " " << par->kickE << endl << endl;
//
//        }
    
    
    
//    const double terminate_time = 10.0;

    //    parallel_simulate(&pops, terminate_time);

    
//    for (int i=0; i < 3*NUM_SAMPLES; i++){
//        cout << "pi[" << i << "]" << endl;
//        cout << pi[i].SEE << " " << pi[i].SIE << " " << pi[i].SEI << " " << pi[i].SII << endl;
//        cout << pi[i].kickI << " " << pi[i].kickE << endl << endl;
//
//    }
//
//
//    uniform_real_distribution<double> u(0, 1);
//    time_t seconds = time(NULL);
//    mt19937 mt(seconds);
//
//    cout << "seconds " << seconds << endl;
//    for (int i=0; i < 3*NUM_SAMPLES; i++){
//        cout << u(mt) << endl;
//    }
//
}


