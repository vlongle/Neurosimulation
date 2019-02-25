#ifndef Simulator_H_
#define Simulator_H_

#include "Vector.h"
#include "randomEngine.h"


#include <random>

struct Parameters{
    // strength of connetion: synaptic weight
    double SEE, SEI, SIE, SII;
    // Poisson process mean of external current
    double kickE, kickI;
    
    int NE, NI;
};

struct Result{
    int E_spike;
    int I_spike;
};

class Population{
//    private:
//        struct Result result;
//        struct Parameters* parameters;
    public:
        struct Result result = {0, 0};
        struct Parameters* parameters;
    
        struct RandomEngine randEng;
    
        // internal data structures
        // HEE: the index of the receiving neurons in pending neuron-to-neuron kicks.
        // Eref: store index of excitatory neurons in refraction
        // VE: voltage of each excitatory neurons
        Vector<int> HEE, HEI, HII, HIE, Eref, Iref, VE, VI;
        // mean of Poisson processes:
        //      0 Edrive, 1 Idrive, 2 HEE, 3, HEI, 4, HIE, 5, HII, 6, Eref, 7, Iref
        Vector<double> Clock;
    
        // boolean array of neurons population: awake 1, not awake 0
        Vector<int> awakeE, awakeI;

        vector<double> time_spike;
        vector<int> num_spike;
    
        // maximum voltage range
        int level = 100;
        // minimum voltage
        int reverse = -66;

        int Level = 100;
        int Reverse = -66;
    
    
        // tau_r the mean of exponential r.v. of refractory length
        double Ref = 250.0;
    
        // connectivity: probability of a random connection, determined by coin flips
        // in each spike --> why don't pre-determine? Memory issue?
        double PEE = 0.15;
        double PIE = 0.5; // probability that an pre-synaptic E connected to post-synaptic I is 0.5
        double PEI = 0.5;
        double PII = 0.4;
    
        // delay in neuron-neuron kicks' effects. Also exponentially distributed
        double HitE = 1000.0; // tauE. This is exhibitory kick delay
        double HitI = 1000.0; // tauI. This is inhibitory kick delay
    
    void set_parameters(struct Parameters* p){
        parameters = p;
    }
    
    void set_random_engine(int thread_num);
    
        void simulate(const double terminate_time);
        void test_simulate(const double terminate_time);

        void prepare_clock();
        int get_voltage_load(const double x);
    
    void spikeE(const int neuron_index);
    void spikeI(const int neuron_index);
    
    void current_at_excitatory(const int current_time);
    void current_at_inhibitory(const int current_time);
    
    // HEE: current from excitatory to excitatory neurons
    void E_to_E(const int current_time);
    void E_to_I(const double current_time); // HIE
    void I_to_E(); // HEI
    void I_to_I(); // HII
    
    // Eref, Iref cases. Excitatatory/inhibitory goes into refractory period
    void excitatory_awake();
    void inhibitory_awake();

    
    void update(const double terminate_time);

         
};

#endif
