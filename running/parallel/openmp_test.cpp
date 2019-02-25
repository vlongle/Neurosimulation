#include <omp.h>
#include <iostream>

using namespace std;

int main(){
    #pragma omp parallel
    {    
    int ID = omp_get_thread_num();
    cout << "hello " << ID << endl;
    cout << "world " << ID << endl;
    }

}
