#include "FMsketch/FMsketch.h"
#include "MurmurHash3/MurmurHash3.h"
#include "build/FMSConfig.h"
#include "noise/noise.h"
#include "FMSsketch.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <math.h>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <time.h>
#include <vector>
#include "LLSketch/LLSketch.h"

using namespace std;
using namespace chrono;

extern double output[7];
extern double batch_fm_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m,double union_inter[]);
extern double batch_noisy_fm_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, int noisy_type ,double epsilon,double union_inter[]);

extern double batch_ll_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m,double union_inter[]);
extern double batch_noisy_ll_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, int noisy_type, double epsilon, double union_inter[]);

extern double batch_fms_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, double union_inter[]); //  n1: car of set S1;  n2: car of set S2;  n3: car of set S1 union S2;  n4: car of set S1 inter S2
extern double batch_noisy_fms_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, double res[], int noisy_type, double epsilon); //  n1: car of set S1;  n2: car of set S2;  n3: car of set S1 union S2;  n4: car of set S1 inter S2
