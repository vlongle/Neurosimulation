#include "Vector.h"
#include <iostream>
using namespace std;

int main(){

    Vector<double> test;
    test.reserve(10);
    cout << "main" << endl;
    cout << test.size() << endl;
    for (int i=0; i < 10; i++){
//        cout << "loop" << endl;
        
        test.put(i, i*2); 
        cout << "test[" << i << "]: " << test[i] << endl;
    }
    
    cout << "end" << endl;
}
