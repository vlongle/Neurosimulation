#include "simulator.h"
#include <omp.h>

// spikeE and spikeI might be BUGGY



// calculate the neuron-neuron voltage delivery load
// Rule: round up if the fractional part of voltage_load exceeds a random Bernoulli flip
int Population::get_voltage_load(const double voltage_load){
    int load_int = floor(voltage_load); // integral part
    double q = voltage_load - (double)load_int; // fractional part
    if (randEng.get_random_number() < q){ // randomly generate a fraction
        return load_int + 1;
        
    }
    else{
        return load_int;
    }
}


void Population::set_random_engine(int thread_num){
    time_t seconds = time(NULL);
    mt19937 mt(seconds + (double)thread_num);
    uniform_real_distribution<double> u(0, 1);
    randEng.set_random_generator(mt, u);
    
}

void Population::spikeE(const int neuron_index){
//    #pragma omp critical
//    cout << "spikeE !!!" << endl;
    result.E_spike ++;
    VE[neuron_index] = 0;
    awakeE[neuron_index] = 0;
    Eref.push_back(neuron_index);
    Clock.put(6, Ref*Eref.size());
    
    for(int i=0; i < parameters->NE; i++){
        // independently flipping coin for every neuron.
        if(randEng.get_random_number() < PEE && awakeE[i]){
            HEE.push_back(i);
            
        }
    }
    
    for(int i = 0; i < parameters->NI; i++)
    {
        if(randEng.get_random_number()  < PIE && awakeI[i])
        {
            HIE.push_back(i);
        }
    }
    Clock.put(2, HitE*HEE.size());
    Clock.put(3, HitE*HIE.size());
    
    
}

void Population::spikeI(const int neuron_index){
//#pragma omp critical
//    cout << "spikeI !!!" << endl;
    result.I_spike ++;
    VI[neuron_index] = 0;
    awakeI[neuron_index] = 0;
    Iref.push_back(neuron_index);
    Clock.put(7, Ref*Iref.size());
    
    for(int i = 0; i < parameters->NE; i++)
    {
        if(randEng.get_random_number() < PEI && awakeE[i])
        {
            HEI.push_back(i);
        }
    }
    for(int i = 0; i < parameters->NI; i++)
    {
        if(randEng.get_random_number() < PII && awakeI[i])
        {
            HII.push_back(i);
        }
    }
    Clock.put(4, HitI*HEI.size());
    Clock.put(5, HitI*HII.size());
    
    
}

void Population:: prepare_clock(){
    int NE = parameters->NE;
    int NI = parameters->NI;
    double kickE = parameters->kickE;
    double kickI = parameters->kickI;

    HEE.reserve(100000);
    HEI.reserve(100000);
    HII.reserve(100000);
    HIE.reserve(100000);
    time_spike.reserve(100000);
    num_spike.reserve(100000);

    awakeE.reserve(NE);
    awakeI.reserve(NI);

    VE.reserve(NE);
    VI.reserve(NI);
    
    
    Clock.reserve(8);
    
    Eref.reserve(NE);
    Iref.reserve(NI);
    // initialization
    for (auto& i : awakeE) i = 1;
    for (auto& i : awakeI) i = 1;
    
    for(auto i : VE) i = 0;
    for(auto i : VI) i = 0;
    
    Clock.push_back(NE * kickE);
    Clock.push_back(NI * kickI);
    
#pragma omp critical
    cout << "push back " << NE << " " << NI << " " << kickE << " " << kickI << endl;
    // init: NE*kickE, NI*kichI, 0,       0,      0,      0,       0,      0

    for (int i = 2; i < 8; i++) Clock.push_back(0);
    Clock.update_sum();
    
    
    
}


void Population::current_at_excitatory(const int current_time){
    int receiving_neuron_index = floor(randEng.get_random_number()*parameters->NE);
    
    if (awakeE[receiving_neuron_index]){
        VE[receiving_neuron_index] ++;
        
        
        if(VE[receiving_neuron_index] >= level)
        {
            spikeE(receiving_neuron_index);
            time_spike.push_back(current_time);
            num_spike.push_back(receiving_neuron_index);

        }
        
        
    }
}

void Population::current_at_inhibitory(const int current_time){
    int receiving_neuron_index = floor(randEng.get_random_number()*parameters->NI);
    if (awakeE[receiving_neuron_index]){
        VI[receiving_neuron_index] ++;
        
        
        if(VI[receiving_neuron_index] >= level)
        {
            spikeI(receiving_neuron_index);
            time_spike.push_back(current_time);
            num_spike.push_back(receiving_neuron_index + parameters->NE);
            
        }
        
}
}

void Population::E_to_E(const int current_time){
    
    int index = floor(randEng.get_random_number() * HEE.size());
    
    int receiving_neuron_index = HEE.pop(index);
    if (awakeE[receiving_neuron_index]){
        VI[receiving_neuron_index] ++;
        
        if( VI[receiving_neuron_index] >= level)
        {
            spikeI(receiving_neuron_index);
            time_spike.push_back(current_time);
            num_spike.push_back(receiving_neuron_index + parameters->NE);
        }
        
        }
    Clock.put(2, HitE * HEE.size());

    
}


void Population::E_to_I(const double current_time){
    int index = floor(randEng.get_random_number() * HIE.size());
    int receiving_neuron_index = HIE.pop(index);
    
    if (awakeI[receiving_neuron_index]) {
        VI[receiving_neuron_index] += get_voltage_load(parameters->SIE);
        if (VI[receiving_neuron_index] >= level) {
            spikeI(receiving_neuron_index);
            time_spike.push_back(current_time);
            num_spike.push_back(receiving_neuron_index + parameters->NE);
        }
    }
    Clock.put(3, HitE * HIE.size());

}


void Population::I_to_E(){
    
    int index = floor(randEng.get_random_number() * HEI.size());
    int receiving_neuron_index = HEI.pop(index);
    
    
    if (awakeE[receiving_neuron_index]) {
        VE[receiving_neuron_index] += get_voltage_load(parameters->SEI);
        // note that there is no spiking here. This is because Inhibitory current
        // would depress our receiving_neuron's voltage.
        if (VE[receiving_neuron_index] < reverse){
            VE[receiving_neuron_index] = reverse;
        }
    }
    Clock.put(4, HitI * HEI.size());
    
}



void Population::I_to_I(){
    
    int index = floor(randEng.get_random_number() * HII.size());
    int receiving_neuron_index = HII.pop(index);
    
    if (awakeI[receiving_neuron_index]) {
        
        VI[receiving_neuron_index] += get_voltage_load(parameters->SII);
        if (VE[receiving_neuron_index] < reverse){
            VE[receiving_neuron_index] = reverse;
        }
    }
    Clock.put(5, HitI * HII.size());
}

void Population::excitatory_awake(){
    
    int index = floor(randEng.get_random_number() * Eref.size());
    int neuron_index = Eref.pop(index);
    
    awakeE[neuron_index] = 1;
    Clock.put(6, Ref * Eref.size());

    
}

void Population::inhibitory_awake(){
    
    int index = floor(randEng.get_random_number() * Iref.size());
    int neuron_index = Iref.pop(index);
    
    awakeE[neuron_index] = 1;
    Clock.put(7, Ref * Iref.size());

    
    
}


void Population:: simulate(const double terminate_time){
    double current_time = 0.0;
    
//    mt19937 mt = randEng->mt;
//    uniform_real_distribution<double> u = randEng->u;
    
    prepare_clock();
    
    while (current_time < terminate_time){

        double mean_sum = Clock.get_sum();

        // inverse CDF method for generating exponential dist with mean mean_sum
        
        double random_num = randEng.get_random_number();
//        #pragma omp critical
//        cout << "random " << random_num << " mean_sum " << mean_sum << endl;

        current_time += -log(1 - random_num)/mean_sum;
        
        int index = randEng.weighted_random_select(&Clock);
        
        switch (index) {
            case 0:
                current_at_excitatory(current_time);
                break;

            case 1:
                current_at_inhibitory(current_time);
                break;

            case 2:
                E_to_E(current_time);
                break;

            case 3:
                E_to_I(current_time);
                break;


            case 4:
                I_to_E();
                break;

            case 5:
                I_to_I();
                break;

            case 6:
                excitatory_awake();
                break;

            case 7:
                inhibitory_awake();
                break;

        }
        
        
//        #pragma omp critical
//        cout << "current_time:" << current_time << endl;
    }
    
}



void Population::test_simulate(const double terminate_time){
    double current_time = 0.0;
    prepare_clock();
    
    while (current_time < terminate_time){
        
        double mean_sum = Clock.get_sum();
        double random_num = randEng.get_random_number();
        int thread_num = omp_get_thread_num();
        
        
        
        //        double time_incr =  -log(1 - random_num)/mean_sum;
        double time_incr =  -log(1 - random_num);
        
        current_time += time_incr;
        
#pragma omp critical
        cout << "thread_num " << thread_num << " time increment " << time_incr <<
        " mean_sum " << mean_sum << endl;
        
        
    }
    
}
    
    void Population::update(const double terminate_time)
    {
        
        
        // This MIGHT be BUGGY !!!!!
        prepare_clock();

        
        
        int NE = parameters->NE;
        int NI = parameters->NI;
        double SEE = parameters->SEE;
        double SEI = parameters->SEI;
        double SIE = parameters->SIE;
        double SII = parameters->SII;

//        int NE = 300;
//        int NI = 100;
//        double SEE = 5;
//        double SIE = 3;
//        double SEI = -2;
//        double SII = -2;
        
        
        
        double current_time = 0.0;
        int count = 0;
        
        cout << "in update " << endl;
        
        mt19937 mt = randEng.randGen.mt;
        uniform_real_distribution<double> u = randEng.randGen.u;
        
        
        while(current_time < terminate_time)
        {
            
//            cout << "random number: " << randEng.get_random_number() << endl;
//            current_time += -log(1 - randEng.get_random_number())/Clock.get_sum();
            
            
            current_time += -log(1 - u(mt))/Clock.get_sum();
            
            
            
//            cout << "current_time " << current_time << endl;
//            cout << "clock.get_sum " << Clock.get_sum() << endl;
//
//            cout << "randEng " << &randEng << endl;
            int index =         randEng.weighted_random_select(&Clock);
            
//            cout << "index " << index << endl;

            
            
            int whichHit;
            int local_index;
            
            int index_neuron;
            
//            cout << "ready to index current_time " << current_time << endl;


            switch (index)
            {
                case 0:
//                    whichHit = floor(u(mt)*NE);
                    
                    whichHit = floor(randEng.get_random_number()*NE);
                    if(awakeE[whichHit])
                    {
                        VE[whichHit]++;
                        if(VE[whichHit] >= Level)
                        {
                            spikeE(whichHit);
                            time_spike.push_back(current_time);
                            num_spike.push_back(whichHit);
                        }
                    }
                    break;
                case 1:
//                    whichHit = floor(u(mt)*NI);
                    whichHit = floor(randEng.get_random_number()*NI);

                    
                    
                    if(awakeI[whichHit])
                    {
                        VI[whichHit]++;
                        if(VI[whichHit]>= Level)
                        {
                            spikeI(whichHit);
                            time_spike.push_back(current_time);
                            num_spike.push_back(whichHit + NE);
                        }
                    }
                    break;
                case 2:
//                    whichHit = HEE.select(mt, u);
                    //                cout<<"ID = "<<whichHit<<" V = "<<VE[whichHit]<<endl<<" status = "<<awakeE[whichHit]<<endl;
                    
                    index_neuron = floor(randEng.get_random_number() * HEE.size());
                    
                    whichHit = HEE.pop(index);
                    
                    if(awakeE[whichHit])
                    {
                        
                        VE[whichHit] += get_voltage_load(SEE);
                        if( VE[whichHit] >= Level)
                        {
                            spikeE(whichHit);
                            time_spike.push_back(current_time);
                            num_spike.push_back(whichHit);
                        }
                    }
                    //                cout<<" after "<<VE[whichHit]<<endl;
                    Clock.put(2, HitE*HEE.size());
                    break;
                case 3:
                    whichHit = HIE.select(mt, u);
                    //                cout<<HEE.size()<<endl;
                    
//                    index_neuron = floor(randEng.get_random_number() * HIE.size());
//
//                    whichHit = HIE.pop(index);
                    
                    if(awakeI[whichHit])
                    {
                        VI[whichHit] += get_voltage_load(SIE);
                        if( VI[whichHit] >= Level)
                        {
                            spikeI(whichHit);
                            time_spike.push_back(current_time);
                            num_spike.push_back(whichHit + NE);
                        }
                    }
                    Clock.put(3, HitE*HIE.size());
                    break;
                case 4:
                    whichHit = HEI.select(mt, u);
                    
                    
//                    index_neuron = floor(randEng.get_random_number() * HEI.size());
//
//                    whichHit = HEI.pop(index);
                    
                    
                    if(awakeE[whichHit])
                    {
                        VE[whichHit] += get_voltage_load(SEI);
                        if(VE[whichHit] < Reverse)
                            VE[whichHit] = Reverse;
                    }
                    Clock.put(4, HitI*HEI.size());
                    break;
                case 5:
                    
                    whichHit = HII.select(mt, u);
                    
//                    index_neuron = floor(randEng.get_random_number() * HII.size());
//
//                    whichHit = HII.pop(index);
                    
                    
                    if(awakeI[whichHit])
                    {
                        VI[whichHit] += get_voltage_load(SII);
                        if(VE[whichHit] < Reverse)
                            VE[whichHit] = Reverse;
                    }
                    Clock.put(5, HitI*HII.size());
                    break;
                    
                case 6:
                    whichHit = Eref.select(mt, u);
                    
//                    index_neuron = floor(randEng.get_random_number() * Eref.size());
//
//                    whichHit = Eref.pop(index);
                    
                    awakeE[whichHit] = 1;
                    Clock.put(6, Ref*Eref.size());
                    break;
                case 7:
                    whichHit = Iref.select(mt, u);
                    
                    
//                    index_neuron = floor(randEng.get_random_number() * Iref.size());
//
//                    whichHit = Iref.pop(index);
                    
                    
                    awakeI[whichHit] = 1;
                    Clock.put(7, Ref*Iref.size());
                    break;
            }
        }
    
    }
    
    
    





// refractor every class that has to do with random stuff into a randomEng class

