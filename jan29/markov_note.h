/* jan 29, 2019. Annotation of Prof. Yao Li's code for generating training data */

// parameters
// number of neurons
int NE = 300;
int NI = 100;


// voltage scale
int level = 100;
int Reverse = -66; // What is M_r: the inhibitory reversal potential? Why is it lower-bounded?

// connectivity
// prob of connection
double PEE = 0.15;
double PIE = 0.5;
double PEI = 0.5;
double PII = 0.4;
// strength of connection (synaptic weight)
double SEE = 5;
double SIE = 3;
double SEI = -2;
double SII = -2;

// random variables
// mean of external current (Poisson process)
double kickE = 1000.0;
double kickI = 1000.0;
// mean of the delay in neuron-neuron kicks' effects. Exponentially dist
double HitE = 1000.0;
double HitI = 1000.0;
// mean of the refractory length. Exponentially dist
double Ref = 250.0;

// summary statistics to calculate the mean firing rates
int E_spike = 0;
int I_spike = 0;


int real2int(double neuron-neuron_voltage, bernoulli_dist)
return resulting_voltage;



// helper vector
class vector{
    private sum;
    
    // switch_element
    void put(index, element); // put an element at the specified index
    bool remove_elt(index);
    
    // uniformly select a neuron stored in vector and remove it
    int select(index, uniform_dist);
}

// ?
int find_index(clock, uniform_dist_on_sum_of_mean);

void spikeE();
