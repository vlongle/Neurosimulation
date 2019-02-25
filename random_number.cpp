#include <random>
#include <iostream>
#include <ctime>
using namespace std;

int main(){
    /* roll dice 0-5.
     */
    
    default_random_engine random_generator(time(0));
    uniform_int_distribution<int> dice_roll(1,6); // int range from 1 to 6
    
    cout << "Dice roll: " << dice_roll(random_generator) << endl;
    cout << "Dice roll: " << dice_roll(random_generator) << endl;

}
