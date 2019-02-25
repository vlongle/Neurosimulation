// Annotate Prof Yao Li's code


/* parameters */
int NE = 300; // number of excitatory neurons
int NI = 100; // number of inhibitory neurons

// connectivity: probability of a random connection, determined by coin flips
// in each spike --> why don't pre-determine? Memory issue?
double PEE = 0.15;
double PIE = 0.5; // probability that an pre-synaptic E connected to post-synaptic I is 0.5
double PEI = 0.5;
double PII = 0.4;

// the voltage scale
int Reverse = -66; // M_R: the reversal voltage?
// OLD: M (maximal) voltage that spikes the neuron in the paper is divided into M_E and M_I

// what is E_spike, I_spike??
// more like summary statistics to calculate the E and I spike rate.
int E_spike = 0;
int I_spike = 0;
// level can be the maximal voltage here.



// strength of connection. Synaptic weight
// note that pre E tends to increase post neurons' voltage.
// Weird: while pre E after spike increases V by 1, the V of affected post neuron
// can increase as high as 5.
double SEE = 5;
double SIE = 3;

// the effect of pre I depends on V of post.
// Formula: the V of post after effect = (V_post + M_r)/(M+M_R)*S_{,I}
double SEI = -2;
double SII = -2;


// Poisson process of arriving external current
double kickE = 1000.0; // lambda_E
double kickI = 1000.0; // lambda_I


// delay in neuron-neuron kicks' effects. Also exponentially distributed
double HitE = 1000.0; // tauE
double HitI = 1000.0; // tauI



/* unknown parameters */
int level = 100; // ???
// where is tau_R: the refractory length of a neuron as exponential r.v.???
double Ref = 250.0;

/* functions */
// calculate the neuron-neuron voltage deliveryd
// Rule: round up if the fractional part of x exceeds a random Bernoulli flip
int real2int(const double x, mt19937& mt, uniform_real_distribution<double>& u){
    int xf = floor(x); // integral part
    double q = x - (double)xf; // fractional part
    if (u(mt) < q){ // randomly generate a fraction
        return xf + 1;
        
    }
    else{
        return xf;
    }
}


// a variable-length array (vector) to store incoming external current.
template <class T>
class Vector<T>{
private:
    T sum;
public:
    double get_sum();
    
    void maintain(); // update the sum by looping through array.
    
    void switch_element(index, T element); // more like replace the val at index with
    // element, and update the sum (maintain).
    void remove_elt(index);
    
    // returning index. Also REMOVE that element as well.
    int select(uniform_real_distribution); // uniformly select a neuron
    int find_min();
    
}

// main functions

// used to randomly choose an action to take during each timestep.
int find_index(Vector<double>& array, uniform_real_distribution);

// seems like clock stores a bunch of mean of 8 different exponential distributions.
int main(){
    
    uniform_real_distribution u(0,1); // uniform from 0 to 1
    
    // various tables to keep info
    Vector<int> VE(NE), VI(NI); // voltage of E and I neurons at any time
    
    // awakened neurons. False means the neuron is in refractory stage.
    Vector<bool> awakeE(NE), awakeI(NI);
    
    // the number of pending neuron-to-neuron kicks. What's int??
    // discrete-time?
    Vector<int> HEE, HEI, HII, HIE;
    
    // refractory neurons?? What's different from awakeI????
    Vector<int> Eref, Iref;
    
    // what are HEE, HEI, HIE, HII ??
    Vector<double> Clock(8); // ????
    //      0 Edrive, 1 Idrive, 2 HEE, 3, HEI, 4, HIE, 5, HII, 6, Eref, 7, Iref
    // init: NE*kickE, NI*kichI, 0,       0,      0,      0,       0,      0
    
    // Summary stats?
    Vector<double> time_spike;
    Vector><int> num_spike;
    
    // time simulation
    double terminate_time = 10.0;
    
    
    update(); // most of the ACTION here.
    
}


void spikeI(whichHit); // transfer current to other neurons
void spikeE(whichHit); // transfer current to other neurons


void update(...){
    
    double current_time = 0.0;
    
    // we are using one random variable for multiple thing??
    while (current_time < terminate_time){
        // some archaic way to increment time?
        current_time += -log(1-u(mt))/Clock.get_sum();
        // inverse CDF method for generating exponential dist with mean Clock.get_sum
        // https://www.youtube.com/watch?v=rnBbYsysPaU
        index = random(0, 7);
        
        switch(index){
        // current arrived at Excitatory neurons
            case Edrive:
                whichHit = random(0, NE)
                if (whichHit.isAwake){
                    VE[whichHit] ++; // increment voltage by 1
                    
                    // spike
                    if (VE[whichHit] >= level){
                        spikeE(whichHit);
                        time_spike.append(current_time);
                        num_spike.append(whichHit);
                    }
                    
                    
                }
                
            // current arrived at Inhibitory neurons
            case IDrive:
                whichHit = random(0, NI);
                if (whicHit.isAwake){
                    VI[whichHit] ++;
                    if (VI[whichHit] >= level){
                        spikeI(whichHit);
                        time_spike.append(current_time);
                        num_spike.append(whichHit);
                    }
                }
                clock.put(2, HitE*HEE.size());
                
            case HEE:
                whichHit = HEE.select(mt, u);
                if (whicHit.isAwake){
                    VE[whichHit] += real2int(SEE); // this adds some randomness to
                    // the resulting neuron-to-neuron current. Why?
                    if (VI[whichHit] >= level){
                        spikeI(whichHit);
                        time_spike.append(current_time);
                        num_spike.append(whichHit);
                    }
                                }
            case HEI:
                ...
            case Eref:
                whichHit = Eref.random()
                awakeE[whichHit] = True;
                ...
                
                
        }
        
    }
    
}
