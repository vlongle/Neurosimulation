#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <cstdlib>
#include <float.h>
#include <random>
#include <stdint.h>
#include <climits>
#include <vector>

#include <omp.h>
#include <chrono>
#include <ctime>


#define NUM_SAMPLES 5

#define _USE_MATH_DEFINES
using namespace std;

int NE = 300;
int NI = 100;
int Level = 100;
double PEE = 0.15;
double PIE = 0.5;
double PEI = 0.5;
double PII = 0.4;

double Ref = 250.0;
double HitE = 1000.0;
double HitI = 1000.0;
int Reverse = -66;


struct Result{
    int E_spike;
    int I_spike;
};

struct Parameters{
    // strength of connetion: synaptic weight
    double SEE, SEI, SIE, SII;
    // Poisson process mean of external current
    double kickE, kickI;
};



int real2int(const double x, mt19937& mt, uniform_real_distribution<double>& u)
{
    int xf = floor(x);
    double q = x - (double)xf;
    double y = 0;
    if( u(mt) < q)
        y = xf + 1;
    else
        y = xf;
    return y;
}

template <class T>
class Vector : public vector<T> {
private:
    T sum;
public:

    double get_sum()
    {
        return sum;
    }

    void maintain()
    {
        T tmp_sum = 0;
        for(auto && it : (*this))
            tmp_sum += it;
 //       cout<<tmp_sum<<endl;
        sum = tmp_sum;
    }

    void switch_element(const int index, T element)
    //switch [index] by element and maintain sum
    {
        sum = sum - (*this)[index] + element;
        (*this)[index] = element;
    }
    //The first three function is only for Clock[]

    bool remove_elt(unsigned index)
    //remove element [index]
    {
        if(this -> empty() || index < 0 || index >= this -> size()) {
            cout<<"Empty vector"<<endl;
            return false;
        }
        if(index == this->size() - 1) {
            this->pop_back();
        }
        else {
            T elt = (*this)[(*this).size() - 1];
            (*this)[index] = elt;
            (*this).pop_back();
        }
        return true;
    }

    int select(mt19937& mt, uniform_real_distribution<double>& u)
    {
        if( this->size() == 0)
        {
            cout<<"Size cannot be zero "<<endl;
            return -1;
        }
        else
        {
            int index = floor(u(mt)*this->size());
            T tmp =  (*this)[index];
            remove_elt(index);
            return tmp;
        }
    }

    int find_min() {
        if(!this -> empty()) {
            auto min = (*this)[0];
            auto index = 0;
            auto iter = 0;
            for(auto&& item : (*this)) {
                if(item < min) {
                    min = item;
                    index = iter;
                }
                iter++;
            }
            return index;
        }
        return -1;
    }

    void print_vector() {
        for(auto && it : (*this)) {
            cout<<it<<" ";
        }
        cout<<endl;
    }

};


int find_index(Vector<double>& array, mt19937& mt, uniform_real_distribution<double>& u)
{
    double sum = array.get_sum();
    double tmp = u(mt);
    int index = 0;
    int size = array.size();
    double tmp_double = array[0]/sum;
//    cout<<tmp_double<<" "<<tmp<<endl;
    while(tmp_double < tmp && index < size - 1)
    {
        index++;
        tmp_double+=array[index]/sum;
    }
    return index;
}

void spikeE(const int whichHit, int* E_spike, Vector<double>& Clock, vector<int> &VE, Vector<int> &HEE, Vector<int> &HIE, Vector<int> &Eref, vector<int> &awakeE, vector<int> &awakeI, mt19937& mt, uniform_real_distribution<double>& u)
{

//    cout << "spikeE!!" << endl;

    (*E_spike) ++;
    VE[whichHit] = 0;
    awakeE[whichHit] = 0;
    Eref.push_back(whichHit);
    Clock.switch_element(6, Ref*Eref.size());
    for(int i = 0; i < NE; i++)
    {
        // independently flipping coin for every neuron.
        if(u(mt) < PEE && awakeE[i])
        {
            HEE.push_back(i);
        }
    }
    for(int i = 0; i < NI; i++)
    {
        if(u(mt) < PIE && awakeI[i])
        {
            HIE.push_back(i);
        }
    }
    Clock.switch_element(2, HitE*HEE.size());
    Clock.switch_element(3, HitE*HIE.size());

}



void spikeI(const int whichHit, int* I_spike, Vector<double>& Clock, vector<int> &VI, Vector<int> &HEI, Vector<int> &HII, Vector<int> &Iref, vector<int> &awakeE, vector<int> &awakeI, mt19937& mt, uniform_real_distribution<double>& u)
{

//    cout << "spikeI!!" << endl;


    (*I_spike) ++;
    VI[whichHit] = 0;
    awakeI[whichHit] = 0;
    Iref.push_back(whichHit);
    Clock.switch_element(7, Ref*Iref.size());
    for(int i = 0; i < NE; i++)
    {
        if(u(mt) < PEI && awakeE[i])
        {
            HEI.push_back(i);
        }
    }
    for(int i = 0; i < NI; i++)
    {
        if(u(mt) < PII && awakeI[i])
        {
            HII.push_back(i);
        }
    }
    Clock.switch_element(4, HitI*HEI.size());
    Clock.switch_element(5, HitI*HII.size());
}

void update(Result *res, vector<double>& time_spike, vector<int>& num_spike, Vector<double>& Clock, vector<int> &VE, vector<int> &VI, Vector<int> &HEE, Vector<int> &HEI, Vector<int> &HIE, Vector<int> &HII, Vector<int> &Eref, Vector<int> &Iref,
    vector<int> &awakeE, vector<int> &awakeI, const double terminate_time, mt19937& mt, uniform_real_distribution<double>& u, double SEE, double SIE, double SEI, double SII)
{
    double current_time = 0.0;
    int count = 0;
    while(current_time < terminate_time)
    {
        current_time += -log(1 - u(mt))/Clock.get_sum();
        int index =         find_index(Clock, mt, u);
        int whichHit;
        count ++;
        int local_index;
//            cout<<"time "<<current_time <<" index "<<index<<endl;
        switch (index)
        {
            case 0:
                whichHit = floor(u(mt)*NE);
                if(awakeE[whichHit])
                {
                    VE[whichHit]++;
                    if(VE[whichHit] >= Level)
                    {
                        spikeE(whichHit, &(*res).E_spike, Clock, VE, HEE, HIE, Eref, awakeE, awakeI, mt, u);
                        time_spike.push_back(current_time);
                        num_spike.push_back(whichHit);
                    }
                }
                break;
            case 1:
                whichHit = floor(u(mt)*NI);
                if(awakeI[whichHit])
                {
                    VI[whichHit]++;
                    if(VI[whichHit]>= Level)
                    {
                        spikeI(whichHit, &(*res).I_spike, Clock, VI, HEI, HII, Iref, awakeE, awakeI, mt, u);
                        time_spike.push_back(current_time);
                        num_spike.push_back(whichHit + NE);
                    }
                }
                break;
            case 2:
                whichHit = HEE.select(mt, u);
//                cout<<"ID = "<<whichHit<<" V = "<<VE[whichHit]<<endl<<" status = "<<awakeE[whichHit]<<endl;
                if(awakeE[whichHit])
                {

                    VE[whichHit] += real2int(SEE, mt, u);
                    if( VE[whichHit] >= Level)
                    {
                        spikeE(whichHit, &(*res).E_spike, Clock, VE, HEE, HIE, Eref, awakeE, awakeI, mt, u);
                        time_spike.push_back(current_time);
                        num_spike.push_back(whichHit);
                    }
                }
//                cout<<" after "<<VE[whichHit]<<endl;
                Clock.switch_element(2, HitE*HEE.size());
                break;
            case 3:
                whichHit = HIE.select(mt, u);
//                cout<<HEE.size()<<endl;
                if(awakeI[whichHit])
                {
                    VI[whichHit] += real2int(SIE, mt, u);
                    if( VI[whichHit] >= Level)
                    {
                        spikeI(whichHit, &(*res).I_spike, Clock, VI, HEI, HII, Iref, awakeE, awakeI, mt, u);
                        time_spike.push_back(current_time);
                        num_spike.push_back(whichHit + NE);
                    }
                }
                Clock.switch_element(3, HitE*HIE.size());
                break;
            case 4:
                whichHit = HEI.select(mt, u);
                if(awakeE[whichHit])
                {
                    VE[whichHit] += real2int(SEI, mt, u);
                    if(VE[whichHit] < Reverse)
                        VE[whichHit] = Reverse;
                }
                Clock.switch_element(4, HitI*HEI.size());
                break;
            case 5:
                whichHit = HII.select(mt, u);
                if(awakeI[whichHit])
                {
                    VI[whichHit] += real2int(SII, mt, u);
                    if(VE[whichHit] < Reverse)
                        VE[whichHit] = Reverse;
                }
                Clock.switch_element(5, HitI*HII.size());
                break;
            case 6:
                whichHit = Eref.select(mt, u);
                awakeE[whichHit] = 1;
                Clock.switch_element(6, Ref*Eref.size());
                break;
            case 7:
                whichHit = Iref.select(mt, u);
                awakeI[whichHit] = 1;
                Clock.switch_element(7, Ref*Iref.size());
                break;
        }
    }
}


Result simulate(int thread_num, int terminate_time, double kickE, double kickI, double SEE, double SIE, double SEI, double SII){

    time_t seconds = time(NULL);
    mt19937 mt(seconds + (double)thread_num);

    uniform_real_distribution<double> u(0, 1);



    struct Result res = {0, 0};
    // preparing methods
    Vector<int> VE, VI;
    VE.reserve(NE);
    VI.reserve(NI);


    for(auto i : VE)
        i = 0;
    for(auto i : VI)
        i = 0;

        Vector<int> HEE, HEI, HII, HIE, Eref, Iref;
        HEE.reserve(100000);
        HEI.reserve(100000);
        HII.reserve(100000);
        HIE.reserve(100000);
        Eref.reserve(NE);
        Iref.reserve(NI);
        vector<int> awakeE(NE);
        for(auto & i : awakeE)
            i = 1;
        vector<int> awakeI(NI);
        for(auto & i : awakeI)
            i = 1;
        Vector<double> Clock;
        Clock.reserve(8);
        //0 Edrive, 1 Idrive, 2 HEE, 3, HEI, 4, HIE, 5, HII, 6, Eref, 7, Iref
        Clock.push_back(NE*kickE);
        Clock.push_back(NI*kickI);
        for(int i = 2; i < 8; i++)
            Clock.push_back(0);
        Clock.maintain();

        vector<double> time_spike;
        time_spike.reserve(100000);
        vector<int> num_spike;
        num_spike.reserve(100000);

        update(&res, time_spike, num_spike, Clock, VE, VI, HEE, HEI, HIE, HII, Eref, Iref, awakeE,awakeI, terminate_time, mt, u,
         SEE,  SIE,  SEI,  SII);

        return res;
}





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
    for(int i=0; i < NUM_SAMPLES; i++){


        double kickE, kickI;

            kickI = 1000.0;
            kickE = 1000.0;
            
        double SEE = SEE_low + (SEE_high-SEE_low) * u(mt);
        double SIE = SIE_low + (SIE_high-SIE_low) * u(mt);
        double SEI = SEI_low + (SEI_high -SEI_low)* u(mt);
        double SII = SII_low + (SII_high-SII_low) * u(mt);
        
        
            (*pi)[i].kickE = kickE;
            (*pi)[i].kickI = kickI;
            
            (*pi)[i].SEE = SEE;
            (*pi)[i].SIE = SIE;
            (*pi)[i].SEI = SEI;
            (*pi)[i].SII = SII;
            
            
//            cout << "i%3 == 0 --> new stuff " << "SEE " << (*pi)[i].SEE << endl;
        
    }

}

    
    // can be parallelized too
    
    
    #pragma omp parallel for
    for(int i=NUM_SAMPLES; i < 2*NUM_SAMPLES; i++){
        // generate for kickE = 3000
        
            double kickI = 3000.0;
            double kickE = 3000.0;
            
            (*pi)[i].kickE = kickE;
            (*pi)[i].kickI = kickI;
            
            
            (*pi)[i].SEE = (*pi)[i-NUM_SAMPLES].SEE;
            (*pi)[i].SIE = (*pi)[i-NUM_SAMPLES].SIE;
            (*pi)[i].SEI = (*pi)[i-NUM_SAMPLES].SEI;
            (*pi)[i].SII = (*pi)[i-NUM_SAMPLES].SII;
            
//             cout << "kickI 3000.0 --> reusing stuff " << "SEE " << (*pi)[i].SEE << endl;

        }

    
    
    #pragma omp parallel for
    for(int i= 2*NUM_SAMPLES; i < 3*NUM_SAMPLES; i++){
        // generate for kickE = 5000
        
        double kickI = 5000.0;
        double kickE = 5000.0;
        
        (*pi)[i].kickE = kickE;
        (*pi)[i].kickI = kickI;
        
        
        (*pi)[i].SEE = (*pi)[i-2*NUM_SAMPLES].SEE;
        (*pi)[i].SIE = (*pi)[i-2*NUM_SAMPLES].SIE;
        (*pi)[i].SEI = (*pi)[i-2*NUM_SAMPLES].SEI;
        (*pi)[i].SII = (*pi)[i-2*NUM_SAMPLES].SII;
        
//        cout << "kickI 5000.0 --> reusing stuff " << "SEE " << (*pi)[i].SEE << endl;
        
    }
    
    
        
        
  
        
    
    
    
        cout << "Number of threads: " << numThreads << endl;

}



void test_serial(){
    cout << "test serial." << endl;

    double SEE = 7.33249;
    double SIE = 6.84878;
    double SEI = -5.33498;
    double SII = -5.07308;
    double kickE = 1000.0;
    double kickI = 1000.0;
    int terminate_time = 1.0;

    int thread_num = 0;

    struct Result res = simulate(thread_num, terminate_time, kickE, kickI, SEE, SIE, SEI, SII);

    double E_spike = res.E_spike;
    double I_spike = res.I_spike;


    cout << endl << endl << "E_spike, I_spike " << E_spike << " " << I_spike << endl;


    cout<<"E spike rate= "<<(double)E_spike/(terminate_time*NE)<<endl;



    cout<<"I spike rate = "<<(double)I_spike/(terminate_time*NI)<<endl << endl;

}

void generate_data(double terminate_time , Vector<Parameters> &pi, Vector<Result> &results){




        pi.reserve(3*NUM_SAMPLES);


        double SEE_low = 5.0;
        double SEE_high = 5.5;

        double SIE_low = 5.0;
        double SIE_high = 5.5;

        double SEI_low = -3.5;
        double SEI_high = -3.0;

        double SII_low = -3.5;
        double SII_high = -3.0;


//    double SEE_low = 5.0;
//    double SEE_high = 5.5;
//
//    double SIE_low = 5.0;
//    double SIE_high = 5.5;
//
//    double SEI_low = -3.5;
//    double SEI_high = -3.0;
//
//    double SII_low = -3.5;
//    double SII_high = -3.0;
    
    
        // randomly generate parameters pi in parallel
        generate_random_connections(3*NUM_SAMPLES, &pi, SEE_low, SEE_high, SIE_low,
                                    SIE_high, SEI_low, SEI_high, SII_low, SII_high);



        cout << "Main: test done generate_random_connections" << endl;






    // myfile.open("spike_info.txt");



    results.reserve(3*NUM_SAMPLES);


    #pragma omp parallel for
    for(int i=0; i < 3*NUM_SAMPLES; i++){
        int thread_num = omp_get_thread_num();

        // unpacking
        double SEE = pi[i].SEE;
        double SIE = pi[i].SIE;
        double SEI = pi[i].SEI;
        double SII= pi[i].SII;
        double kickE = pi[i].kickE;
        double kickI = pi[i].kickI;

        results[i] = simulate(thread_num, terminate_time, kickE, kickI, SEE, SIE, SEI, SII);

    }

// saving training data
    // cout << "results.size() " << results.size() << endl; // return 0 due to the indexing operator



}


// human-readable file
void save_to_disk_human(const char* file_name, double terminate_time , Vector<Parameters> &pi, Vector<Result> &results){
    cout << "writing to disk ..." << endl;
    ofstream myfile;
    
    myfile.open(file_name, ofstream::app);
    
    auto end = chrono::system_clock::now();
    
    time_t end_time = chrono::system_clock::to_time_t(end);
    
    myfile << "Data generated at " << ctime(&end_time) << endl << endl;
    
    
    // for (int i=0; i < results.size(); i++){
    for (int i=0; i < 3*NUM_SAMPLES; i++){
        
        int E_spike = results[i].E_spike;
        int I_spike = results[i].I_spike;
        
        
        
        myfile << "SEE, SIE, SEI, SII | " << pi[i].SEE << " " << pi[i].SIE << " " << pi[i].SEI << " "
        << pi[i].SII << " | kickE, kickI " << pi[i].kickE << " " << pi[i].kickI << endl;
        
        myfile << "E_spike, I_spike | " << E_spike << " " << I_spike << endl;
        
        
        myfile<<"E spike rate= "<<(double)E_spike/(terminate_time*NE)<<endl;
        
        
        
        myfile<<"I spike rate = "<<(double)I_spike/(terminate_time*NI)<<endl << endl;
        
    }
    
    cout << endl << endl;

    
}


// machine-readable data for neural net
void save_to_disk_machine(const char* file_name, double terminate_time , Vector<Parameters> &pi, Vector<Result> &results){
    cout << "writing to disk ..." << endl;
    ofstream myfile;
    
    myfile.open(file_name, ofstream::app);
    
    auto end = chrono::system_clock::now();
    
    time_t end_time = chrono::system_clock::to_time_t(end);
    
    myfile << "Data generated at " << ctime(&end_time) << endl << endl;
    
    
    // for (int i=0; i < results.size(); i++){
    for (int i=0; i < 3*NUM_SAMPLES; i++){
        
        int E_spike = results[i].E_spike;
        int I_spike = results[i].I_spike;
        
        
        
        myfile << pi[i].SEE << " " << pi[i].SIE << " " << pi[i].SEI << " "
        << pi[i].SII << " " << pi[i].kickE << " " << pi[i].kickI;
        
        myfile << " " << E_spike << " " << I_spike;
        
        
        myfile<<" "<<(double)E_spike/(terminate_time*NE);
        
        
        
        myfile<<" "<<(double)I_spike/(terminate_time*NI)<< endl;
        
    }
    
    cout << endl << endl;
}





int main()
{
    struct timeval t1, t2;
    gettimeofday(&t1,NULL);



    char* file_name_machine = "machine_range.txt";
    char* file_name_human = "human_range.txt";


    // generate data in parallel

    // test_serial();

    
    
    
    
    // parallel simulation
    double terminate_time = 1.0;

    
    
    
    Vector<Parameters> pi;
    Vector<Result> results;

    
    generate_data(terminate_time, pi, results);

    
    
    save_to_disk_machine(file_name_machine, terminate_time, pi, results);
    save_to_disk_human(file_name_human, terminate_time, pi, results);

    
    gettimeofday(&t2, NULL);
    double delta = ((t2.tv_sec  - t1.tv_sec) * 1000000u +
                    t2.tv_usec - t1.tv_usec) / 1.e6;

    cout << "total CPU time = " << delta <<endl;

}
