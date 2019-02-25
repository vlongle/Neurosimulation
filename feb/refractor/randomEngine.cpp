#include "RandomEngine.h"

void RandomEngine::set_random_generator(mt19937 set_mt, uniform_real_distribution<double> set_u){
    randGen.mt = set_mt;
    randGen.u = set_u;
    
}


int RandomEngine::weighted_random_select(Vector<double>* array){
    double sum = (*array).get_sum();
    double tmp = get_random_number();
    int index = 0;
    int size = (*array).size();
    double tmp_double = (*array)[0]/sum;
    while(tmp_double < tmp && index < size - 1)
    {
        index++;
        tmp_double += (*array)[index]/sum;
    }
    return index;
}



// always return 0 for some reasons. WHY????
double RandomEngine::get_random_number(){
//#pragma omp critical
//    cout << "randomEngie random_num " <<randGen.u(randGen.mt) << endl;
    return randGen.u(randGen.mt);
}
