#include <iostream>
#include <fstream>
#include <cmath>
#include <iterator>
#include <ctime>
#include <thread>
#include <atomic>
#include <vector>
#include <typeinfo>
#include <bitset>
#include <assert.h>
#include <math.h>
#include "../MurmurHash3/MurmurHash3.h"
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <time.h>
#include <random>
#include <chrono>
#include <thread>

using namespace std;
using namespace chrono;


class FmSketch
{
   public:  //  class member variables
      unsigned * sk;  //  the sketch
      unsigned * temp_sk;  //  temporary sketch maintained for each of the 24 threads
      std::thread * threads;
      unsigned * seeds;  //  the hash function seed for hashing elements inserted
      unsigned number_of_threads;  //  the actual number of threads in this class instance
      unsigned M;
      unsigned r;  //  log_2(M)
      unsigned w;
      unsigned * R;  //  parameters used for cardinality estimation
      unsigned long long n;
      unsigned long long the_n_base;

    public:  //  class member functions
      FmSketch(unsigned t_M, unsigned t_w, unsigned t_seeds[], unsigned t_number_of_threads);
      FmSketch(FmSketch & sk1, FmSketch & sk2);
      ~FmSketch();
      void clear_sketch();
      unsigned long long rank_bits(unsigned long long bits);
      double generate_sketch(unsigned long long t_n, unsigned long long t_n_base);
      void thread_gene_sketch(unsigned index);
      void thread_sketch_combine();
      double estimate();
      void print_sketch(unsigned t_lines);
      void generate_R();
      double noisy48_estimate(double noise, double real_n);
      double noisy_estimate(double noise, double real_n);
      void change_seed(unsigned t_seeds[]);
};
