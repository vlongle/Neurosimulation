#include <iostream>
using namespace std;

struct Result{
    int E_spike;
    int I_spike;
};

void spikeE(int* E_spike){
    (*E_spike)++;
}

void update(Result *res){
    spikeE(&(*res).E_spike);    
}
// return pointer will have a segfault
Result simulate(){
    struct Result res = {0,0};
    update(&res);    
    return res;
}

int main(){
    struct Result res = simulate();
    cout << "res.E_spike" << res.E_spike << endl;
    cout << "res.I_spike" << res.I_spike << endl;
}
