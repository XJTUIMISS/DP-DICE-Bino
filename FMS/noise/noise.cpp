#include "noise.h"

#define PI 3.14159265358
#define delta 1e-12
#define d_parties 20
#define s_mpcdp 10

double cdp_epsilon_gen(double epsilon)
{
    double cdp_epsilon, temp, b;
    b = pow(-2 * log(delta), 0.5);
    temp = b * b - 4 * 0.5 * (-epsilon);
    if (temp < 0)
    {
        cout << "temp error" << endl;
    }
    else if (temp == 0)
    {
        cdp_epsilon = -b;
    }
    else
    {
        cdp_epsilon = pow(temp, 0.5) - b;
    }
    return cdp_epsilon;
}

double EV_epsilon(double sigma, double sensibility)
{
    double temp1 = 0, temp2 = 0;
    for (int i = 1; i < d_parties; i++)
    {
        temp1 += 5 * exp(-2 * PI * PI * sigma * sigma * i / (i + 1));
        temp2 += 10 * exp(-2 * PI * PI * sigma * sigma * i / (i + 1));
    }
    temp1 += (sensibility * sensibility) / (d_parties * sigma * sigma);
    temp1 = pow(temp1, 0.5);
    temp2 += (abs(sensibility)) / (pow(d_parties, 0.5) * sigma);
    if (temp1 < temp2)
    {
        return temp1;
    }
    return temp2;
}

double cdp_sigma_gen(double epsilon, double sensibility)
{
    double temp1 = 0, temp2 = 0, temp, b, cdp_sigma = 0.5, ev_epsilon, last_ev_epsilon, delta_epsilon, cdp_epsilon;
    cdp_epsilon = cdp_epsilon_gen(epsilon);
    if (sensibility > 100)
    {
        while ((EV_epsilon(cdp_sigma, sensibility) - cdp_epsilon) * (EV_epsilon(cdp_sigma + 1000, sensibility) - cdp_epsilon) > 0)
        {
            cdp_sigma += 1000;
        }
    }
    else
    {
        while ((EV_epsilon(cdp_sigma, sensibility) - cdp_epsilon) * (EV_epsilon(cdp_sigma + 0.001, sensibility) - cdp_epsilon) > 0)
        {
            cdp_sigma += 0.0001;
        }
    }
    if (cdp_sigma <= 0.5)
    {
        cout << "cdp_sigma error" << endl;
    }
    return cdp_sigma;
}

long long mpcdp_t_gen(double epsilon, double sensibility)
{
    long long T;
    double t, a_, b_, a, b, c, temp;
    a_ = 2 * s_mpcdp * sensibility * pow(log(1.25 / delta), 0.5);
    b_ = (15 * s_mpcdp * sensibility * pow(2 * log(10 / delta), 0.5) + 2 * s_mpcdp) / (3 * (1 - delta / 10));
    b_ += (4 * s_mpcdp * sensibility * (log(1.25 / delta) + log(20 / delta) * log(10 / delta))) / 3;
    a = epsilon * epsilon;
    b = -(2 * b_ * epsilon + a_ * a_);
    c = b_ * b_;

    temp = b * b - 4 * a * c;
    if (temp < 0)
    {
        cout << "temp error" << endl;
    }
    else if (temp == 0)
    {
        t = -b / (2 * a);
    }
    else
    {
        t = (pow(temp, 0.5) - b) / (2 * a);
    }
    T = ceil(t);
    return T;
}
// Function to generate discrete Gaussian noise

std::vector<double> generate_gaussian_noise(int size, double epsilon, double sensibility)
{
    // Initialize a vector to store the noise
    std::vector<double> noise(size, 0);
    double sigma = (pow((2.0 * log(1.25 / delta)), 0.5) * sensibility) / (double)epsilon;
    // Set up a random number generator to generate standard normal random variables
    std::default_random_engine random(time(NULL) + rand());
    // normal_distribution<double> distribution(0.0, 1.0);

    std::normal_distribution<double> dis1(0.0, sigma);
    for (int i = 0; i < size; i++)
    {
        noise[i] = dis1(random);
    }
    return noise;
}
std::vector<int> generate_discrete_gaussian_noise(int size, double sigma)
{
    // Initialize a vector to store the noise
    std::vector<int> noise(size, 0);

    // Set up a random number generator to generate standard normal random variables
    std::default_random_engine random(time(NULL));

    std::uniform_real_distribution<double> dis1(0.0, 1.0);
    std::uniform_real_distribution<double> dis2(0.0, 1.0);
    // Generate the noise using Box-Muller transform
    for (int i = 0; i < size; i++)
    {
        double u1 = dis1(random);
        double u2 = dis2(random);
        double z = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
        noise[i] = round(z * sigma); // Round to the nearest integer
    }
    return noise;
}
// Function to generate discrete Gaussian noise satisfied cdp
std::vector<double> generate_cdp_noise(int size, double epsilon, double sensibility)
{
    double cdp_sigma, u1, u2, z;
    // Initialize a vector to store the noise
    std::vector<double> noises(size, 0);
    int noise;
    cdp_sigma = cdp_sigma_gen(epsilon, sensibility);
    std::default_random_engine random(time(NULL) + rand());
    std::uniform_real_distribution<double> dis1(0.0, 1.0);
    std::uniform_real_distribution<double> dis2(0.0, 1.0);
    double calcu;
    // Generate the noise using Box-Muller transform
    for (int i = 0; i < size; i++)
    {
        noise = 0;
        for (int j = 0; j < d_parties; j++)
        {
            u1 = dis1(random);
            u2 = dis2(random);
            z = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
            noise += round(z * cdp_sigma);
        }                 
        noises[i] = noise; // Round to the nearest integer
        calcu += noise * noise;
    }
    calcu = pow(calcu / size, 0.5);
    return noises;
}

// Function to generate noise with Binomial mechanism
std::vector<double> generate_binomial_noise(int size, double epsilon, double sensibility)
{

    double mpcdp_epsilon, mpcdp_sigma, u1, u2, z;
    // Initialize a vector to store the noise
    std::vector<double> noises(size, 0);
    long long noise, T;
    T = mpcdp_t_gen(epsilon, sensibility);
    //  Set up a random number generator
    std::default_random_engine random(time(NULL) + rand());

    // Create a binomial distribution
    std::binomial_distribution<long long> distribution(2 * T, 0.5);
    // Generate random numbers
    for (int i = 0; i < size; i++)
    {
        long long randomValue = distribution(random);
        noises[i] = (double)((randomValue - T) / s_mpcdp);
    }
    return noises;
}
