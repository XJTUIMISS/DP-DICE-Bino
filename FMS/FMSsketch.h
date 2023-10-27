#ifndef _SKETCHES_H_
#define _SKETCHES_H_

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <mutex>
#include <time.h>
#include <random>
#include <chrono>
#include <atomic>
#include "MurmurHash3/MurmurHash3.h"
using namespace std;
using namespace chrono;
class FMS
{
   public:  //  class member variables
      unsigned * sk;  //  the sketch
      unsigned * temp_sk;  //  temporary sketch maintained for each of the 24 threads
      std::thread * threads;
      unsigned seed;  //  the hash function seed for hashing elements inserted
      unsigned number_of_threads;  //  the actual number of threads in this class instance
      unsigned M;
      unsigned r;  //  log_2(M)
      unsigned w;
      unsigned real_w;
      unsigned long long n;
      unsigned long long the_n_base;

    public:  //  class member functions
      FMS(unsigned t_M, unsigned t_w, unsigned t_seed, unsigned t_number_of_threads);
      FMS(FMS & sk1, FMS & sk2);
      ~FMS();
      void clear_sketch();
      unsigned long long rank_bits(unsigned long long bits);
      double generate_sketch(unsigned long long t_n, unsigned long long t_n_base);
      void thread_gene_sketch(unsigned index);
      void thread_sketch_combine();
      double estimate();
      void print_sketch(unsigned t_lines);
      double func_ev(double t_n);
      double d_func_ev(double t_n);  //  function E(V)
      double noisy_func_ev(double t_n);
      double noisy_d_func_ev(double t_n);  //  function E(V)
      unsigned noisy_count_zero();  //  count the number of zero bits
      double px(double x);
      unsigned count_zero();  //  count the number of zero bits
      double noisy_estimate(double noise, double real_n);
      void change_seed(unsigned t_seed);
};

#endif