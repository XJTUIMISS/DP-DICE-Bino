#include <cmath>
#include <iostream>
#include <random>
#include <time.h>
#include <vector>
using namespace std;

std::vector<double> generate_gaussian_noise(int size, double epsilon, double sensibility);
// Function to generate discrete Gaussian noise
std::vector<int> generate_discrete_gaussian_noise(int size, double sigma);
double cdp_epsilon_gen(double epsilon);
double cdp_sigma_gen(double cdp_epsilon,double sensibility);
double EV_epsilon(double sigma, double sensibility);
std::vector<double> generate_cdp_noise(int size, double epsilon, double sensibility);
long long mpcdp_t_gen(double epsilon, double sensibility);
std::vector<double> generate_binomial_noise(int size, double epsilon, double sensibility);