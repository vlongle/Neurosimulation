#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

/* Accept 4 parameters SEE, SIE, SEI, SEE
 ./a.out 5 4 -3 -3

 Output: 3 pairs of (E_mean_firing, I_mean_firing)
 for \lambda_E = \lambda_I = 1000, 3000, 5000 i.e. kichE, kickI
 (mean of the external current time in Poisson time process)
 */

#define _USE_MATH_DEFINES
using namespace std;

int NE = 300;
int NI = 100;
int Level = 100;
double PEE = 0.15;
double PIE = 0.5;
double PEI = 0.5;
double PII = 0.4;
double SEE = 5;
double SIE = 3;
double SEI = -2;
double SII = -2;
double kickE = 1000.0;
double kickI = 1000.0;
double Ref = 250.0;
double HitE = 1000.0;
double HitI = 1000.0;
int Reverse = -66;
int E_spike = 0;
int I_spike = 0;

int real2int(const double x, mt19937& mt,
             uniform_real_distribution<double>& u) {
  int xf = floor(x);
  double q = x - (double)xf;
  double y = 0;
  if (u(mt) < q)
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
  double get_sum() { return sum; }

  void maintain() {
    T tmp_sum = 0;
    for (auto&& it : (*this)) tmp_sum += it;
    //       cout<<tmp_sum<<endl;
    sum = tmp_sum;
  }

  void switch_element(const int index, T element)
  // switch [index] by element and maintain sum
  {
    sum = sum - (*this)[index] + element;
    (*this)[index] = element;
  }
  // The first three function is only for Clock[]

  bool remove_elt(unsigned index)
  // remove element [index]
  {
    if (this->empty() || index < 0 || index >= this->size()) {
      cout << "Empty vector" << endl;
      return false;
    }
    if (index == this->size() - 1) {
      this->pop_back();
    } else {
      T elt = (*this)[(*this).size() - 1];
      (*this)[index] = elt;
      (*this).pop_back();
    }
    return true;
  }

  int select(mt19937& mt, uniform_real_distribution<double>& u) {
    if (this->size() == 0) {
      cout << "Size cannot be zero " << endl;
      return -1;
    } else {
      int index = floor(u(mt) * this->size());
      T tmp = (*this)[index];
      remove_elt(index);
      return tmp;
    }
  }

  int find_min() {
    if (!this->empty()) {
      auto min = (*this)[0];
      auto index = 0;
      auto iter = 0;
      for (auto&& item : (*this)) {
        if (item < min) {
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
    for (auto&& it : (*this)) {
      cout << it << " ";
    }
    cout << endl;
  }
};

int find_index(Vector<double>& array, mt19937& mt,
               uniform_real_distribution<double>& u) {
  double sum = array.get_sum();
  double tmp = u(mt);
  int index = 0;
  int size = array.size();
  double tmp_double = array[0] / sum;
  //    cout<<tmp_double<<" "<<tmp<<endl;
  while (tmp_double < tmp && index < size - 1) {
    index++;
    tmp_double += array[index] / sum;
  }
  return index;
}

void spikeE(const int whichHit, Vector<double>& Clock, vector<int>& VE,
            Vector<int>& HEE, Vector<int>& HIE, Vector<int>& Eref,
            vector<int>& awakeE, vector<int>& awakeI, mt19937& mt,
            uniform_real_distribution<double>& u) {
  E_spike++;
  VE[whichHit] = 0;
  awakeE[whichHit] = 0;
  Eref.push_back(whichHit);
  Clock.switch_element(6, Ref * Eref.size());
  for (int i = 0; i < NE; i++) {
    if (u(mt) < PEE && awakeE[i]) {
      HEE.push_back(i);
    }
  }
  for (int i = 0; i < NI; i++) {
    if (u(mt) < PIE && awakeI[i]) {
      HIE.push_back(i);
    }
  }
  Clock.switch_element(2, HitE * HEE.size());
  Clock.switch_element(3, HitE * HIE.size());
}
void spikeI(const int whichHit, Vector<double>& Clock, vector<int>& VI,
            Vector<int>& HEI, Vector<int>& HII, Vector<int>& Iref,
            vector<int>& awakeE, vector<int>& awakeI, mt19937& mt,
            uniform_real_distribution<double>& u) {
  I_spike++;
  VI[whichHit] = 0;
  awakeI[whichHit] = 0;
  Iref.push_back(whichHit);
  Clock.switch_element(7, Ref * Iref.size());
  for (int i = 0; i < NE; i++) {
    if (u(mt) < PEI && awakeE[i]) {
      HEI.push_back(i);
    }
  }
  for (int i = 0; i < NI; i++) {
    if (u(mt) < PII && awakeI[i]) {
      HII.push_back(i);
    }
  }
  Clock.switch_element(4, HitI * HEI.size());
  Clock.switch_element(5, HitI * HII.size());
}

void update(vector<double>& time_spike, vector<int>& num_spike,
            Vector<double>& Clock, vector<int>& VE, vector<int>& VI,
            Vector<int>& HEE, Vector<int>& HEI, Vector<int>& HIE,
            Vector<int>& HII, Vector<int>& Eref, Vector<int>& Iref,
            vector<int>& awakeE, vector<int>& awakeI,
            const double terminate_time, mt19937& mt,
            uniform_real_distribution<double>& u) {
  double current_time = 0.0;
  int count = 0;
  while (current_time < terminate_time) {
    current_time += -log(1 - u(mt)) / Clock.get_sum();
    int index = find_index(Clock, mt, u);
    int whichHit;
    count++;
    int local_index;
    //            cout<<"time "<<current_time <<" index "<<index<<endl;
    switch (index) {
      case 0:
        whichHit = floor(u(mt) * NE);
        if (awakeE[whichHit]) {
          VE[whichHit]++;
          if (VE[whichHit] >= Level) {
            spikeE(whichHit, Clock, VE, HEE, HIE, Eref, awakeE, awakeI, mt, u);
            time_spike.push_back(current_time);
            num_spike.push_back(whichHit);
          }
        }
        break;
      case 1:
        whichHit = floor(u(mt) * NI);
        if (awakeI[whichHit]) {
          VI[whichHit]++;
          if (VI[whichHit] >= Level) {
            spikeI(whichHit, Clock, VI, HEI, HII, Iref, awakeE, awakeI, mt, u);
            time_spike.push_back(current_time);
            num_spike.push_back(whichHit + NE);
          }
        }
        break;
      case 2:
        whichHit = HEE.select(mt, u);
        //                cout<<"ID = "<<whichHit<<" V =
        //                "<<VE[whichHit]<<endl<<" status =
        //                "<<awakeE[whichHit]<<endl;
        if (awakeE[whichHit]) {
          VE[whichHit] += real2int(SEE, mt, u);
          if (VE[whichHit] >= Level) {
            spikeE(whichHit, Clock, VE, HEE, HIE, Eref, awakeE, awakeI, mt, u);
            time_spike.push_back(current_time);
            num_spike.push_back(whichHit);
          }
        }
        //                cout<<" after "<<VE[whichHit]<<endl;
        Clock.switch_element(2, HitE * HEE.size());
        break;
      case 3:
        whichHit = HIE.select(mt, u);
        //                cout<<HEE.size()<<endl;
        if (awakeI[whichHit]) {
          VI[whichHit] += real2int(SIE, mt, u);
          if (VI[whichHit] >= Level) {
            spikeI(whichHit, Clock, VI, HEI, HII, Iref, awakeE, awakeI, mt, u);
            time_spike.push_back(current_time);
            num_spike.push_back(whichHit + NE);
          }
        }
        Clock.switch_element(3, HitE * HIE.size());
        break;
      case 4:
        whichHit = HEI.select(mt, u);
        if (awakeE[whichHit]) {
          VE[whichHit] += real2int(SEI, mt, u);
          if (VE[whichHit] < Reverse) VE[whichHit] = Reverse;
        }
        Clock.switch_element(4, HitI * HEI.size());
        break;
      case 5:
        whichHit = HII.select(mt, u);
        if (awakeI[whichHit]) {
          VI[whichHit] += real2int(SII, mt, u);
          if (VE[whichHit] < Reverse) VE[whichHit] = Reverse;
        }
        Clock.switch_element(5, HitI * HII.size());
        break;
      case 6:
        whichHit = Eref.select(mt, u);
        awakeE[whichHit] = 1;
        Clock.switch_element(6, Ref * Eref.size());
        break;
      case 7:
        whichHit = Iref.select(mt, u);
        awakeI[whichHit] = 1;
        Clock.switch_element(7, Ref * Iref.size());
        break;
    }
  }
}

void simulate(double kickE, double kickI, double SEE, double SIE, double SEI,
              double SII, double terminate_time) {
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> u(0, 1);
  vector<int> VE(NE);
  vector<int> VI(NI);
  for (auto i : VE) i = 0;
  for (auto i : VI) i = 0;
  Vector<int> HEE, HEI, HII, HIE, Eref, Iref;
  HEE.reserve(100000);
  HEI.reserve(100000);
  HII.reserve(100000);
  HIE.reserve(100000);
  Eref.reserve(NE);
  Iref.reserve(NI);
  vector<int> awakeE(NE);
  for (auto& i : awakeE) i = 1;
  vector<int> awakeI(NI);
  for (auto& i : awakeI) i = 1;
  Vector<double> Clock;
  Clock.reserve(8);
  // 0 Edrive, 1 Idrive, 2 HEE, 3, HEI, 4, HIE, 5, HII, 6, Eref, 7, Iref
  Clock.push_back(NE * kickE);
  Clock.push_back(NI * kickI);
  for (int i = 2; i < 8; i++) Clock.push_back(0);
  Clock.maintain();
  vector<double> time_spike;
  time_spike.reserve(100000);
  vector<int> num_spike;
  num_spike.reserve(100000);
  update(time_spike, num_spike, Clock, VE, VI, HEE, HEI, HIE, HII, Eref, Iref,
         awakeE, awakeI, terminate_time, mt, u);
  //    for(auto& i : VE)
  //        cout<<i<<" ";

  //    cout<<"E spike rate= "<<(double)E_spike/(terminate_time*NE)<<endl;
  //    cout<<"I spike rate = "<<(double)I_spike/(terminate_time*NI)<<endl;
  //
  //
  //
  //    // should be synchronize here.
  //    myfile<< "E spike rate= " << (double)E_spike/(terminate_time*NE)<<" ";
  //    myfile<< "I spike rate= " << (double)I_spike/(terminate_time*NI)<<endl;
  //

  /*
  int spike_count = time_spike.size();
  for(int i = 0; i < spike_count; i++)
  {
      myfile<<time_spike[i]<<"  "<<num_spike[i]<<endl;
  }
  myfile.close();
  */
}

double* parallel_random(int numRandoms, double from, double to) {
  double* randomArr = (double*)malloc(numRandoms * sizeof(double));
  int i;
  uniform_real_distribution<double> u(0, 1);

  int numThreads;

#pragma omp parallel
  {
    numThreads = omp_get_num_threads();

    // seconds since January 1, 1970
    time_t seconds = time(NULL);

    int threadNum = omp_get_thread_num();

    // seed = t + threadNum
    mt19937 mt(seconds + (double)threadNum);

#pragma omp for
    for (i = 0; i < numRandoms; i++) {
      randomArr[i] = from + to * u(mt);
    }
  }

  cout << "number of threads: " << numThreads << endl;
  return randomArr;
}

int main(int argc, char* argv[]) {
  struct timeval t1, t2;
  gettimeofday(&t1, NULL);

  int numRandoms = 100;
  vector<vector<int> > Rates_Sample(numRandoms, vector<int>(6));
  // each row is E rate 1000, I rate 1000, ... E rate 5000, I rate 5000
  // This is in shared memory

  // hard-coded
  double SEE_low = 2.0;
  double SEE_high = 6.0;

  double SIE_low = 2.0;
  double SIE_high = 6.0;

  double SEI_low = -5.0;
  double SEI_high = -1.0;

  double SII_low = -5.0;
  double SII_high = -1.0;

  double* SEEArr = parallel_random(numRandoms, SEE_low, SEE_high);
  double* SEIArr = parallel_random(numRandoms, SEI_low, SEI_high);
  double* SIEArr = parallel_random(numRandoms, SIE_low, SIE_high);
  double* SIIArr = parallel_random(numRandoms, SII_low, SII_high);

  double terminate_time = 1.0;
  ofstream myfile;
  // append mode
  myfile.open("training_set.txt", ofstream::app);
  myfile << endl
         << endl
         << "New training set. Num samples " << numRandoms << endl
         << endl;
  // generate for kickI = 3000 and 5000 as well
#pragma omp parallel for
  for (int i = 0; i < numRandoms; i++) {
    int rank = omp_get_thread_num();
    simulate(kickE, kickI, SEE, SIE, SEI, SII, terminate_time);

#pragma omp critical
    {
      cout << "thread " << omp_get_thread_num() << " / "
           << omp_get_num_threads() - 1 << " done with " << i << endl;
      myfile << "SEE, SIE, SEI, SII: " << SEE << " " << SIE << " " << SEI << " "
             << SII << " " << endl;
      myfile << "kichE, kichI: " << kickE << " " << kickI << " " << endl;
      myfile << "E spike rate= " << (double)E_spike / (terminate_time * NE)
             << " ";
      myfile << "I spike rate= " << (double)I_spike / (terminate_time * NI)
             << endl
             << endl;
    }
  }

  cout << "Terminate_time: " << terminate_time << endl << endl;

  gettimeofday(&t2, NULL);
  double delta =
      ((t2.tv_sec - t1.tv_sec) * 1000000u + t2.tv_usec - t1.tv_usec) / 1.e6;

  cout << "total CPU time = " << delta << endl;
  myfile << "total CPU time = " << delta << endl;

  myfile.close();
}
