#include "batch_est.h"

using namespace std;
using namespace chrono;

double batch_fm_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m)
{
    unsigned seed = clock(), M = m, w = 32, number_of_threads = 24;
    unsigned seeds[M];
    srand(seed);
    for (int j = 0; j < M; j++)
    {
        seed = rand();
        seeds[j] = seed;
        srand(seed);
    }
    FmSketch a = FmSketch(M, w, seeds, number_of_threads);
    FmSketch b = FmSketch(M, w, seeds, number_of_threads);
    vector<double> errors;
    double error;
    for (int i = 0; i < repeat_times; ++i)
    {
        seed = rand() * i;
        srand(seed);
        for (int j = 0; j < M; j++)
        {
            seed = rand();
            seeds[j] = seed;
        }
        a.change_seed(seeds);
        b.change_seed(seeds);
        a.generate_sketch(n1, 0);
        b.generate_sketch(n2, n1 - n4);
        FmSketch c(a, b);
        double estimation = a.estimate() + b.estimate() - c.estimate();
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors.push_back(error);
    }
    error = 0.0;
    for (int i = 0; i < errors.size(); ++i)
    {
        error += errors[i];
    }
    error = pow(error / (double)errors.size(), 0.5);
    a.clear_sketch();
    b.clear_sketch();
    return error;
}

double batch_noisy_fm_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, int noisy_type, double epsilon, double union_inter[])
{
    unsigned seed = clock(), M = m, w = 32, number_of_threads = 24;
    unsigned long long n3 = n1 + n2 - n4;
    unsigned seeds[M];
    double sensibility[3], w_sensibility[3];
    w_sensibility[0] = ceil(log2(n1)) + 6;
    sensibility[0] = w_sensibility[0] * (double)M;
    w_sensibility[1] = ceil(log2(n2)) + 6;
    sensibility[1] = w_sensibility[1] * (double)M;
    w_sensibility[2] = ceil(log2(n3)) + 6;
    sensibility[2] = w_sensibility[2] * (double)M;
    std::vector<double> noises1(repeat_times, 0);
    std::vector<double> noises2(repeat_times, 0);
    std::vector<double> noises3(repeat_times, 0);
    if (noisy_type == 1)
    {
        noises1 = generate_gaussian_noise(repeat_times, epsilon, sensibility[0]);
        noises2 = generate_gaussian_noise(repeat_times, epsilon, sensibility[1]);
        noises3 = generate_gaussian_noise(repeat_times, epsilon, sensibility[2]);
    }
    else if (noisy_type == 2)
    {
        noises1 = generate_cdp_noise(repeat_times, epsilon, sensibility[0]);
        noises2 = generate_cdp_noise(repeat_times, epsilon, sensibility[1]);
        noises3 = generate_cdp_noise(repeat_times, epsilon, sensibility[2]);
    }
    else if (noisy_type == 3)
    {
        noises1 = generate_binomial_noise(repeat_times, epsilon, sensibility[0]);
        noises2 = generate_binomial_noise(repeat_times, epsilon, sensibility[1]);
        noises3 = generate_binomial_noise(repeat_times, epsilon, sensibility[2]);
    }
    double error = 0.0;
    FmSketch a = FmSketch(M, w, seeds, number_of_threads);
    FmSketch b = FmSketch(M, w, seeds, number_of_threads);
    vector<double> errors;
    vector<double> errors_union;
    double error_union;
    for (int i = 0; i < repeat_times; ++i)
    {
        FmSketch c(a, b);
        double a_estimation, b_estimation, c_estimation;
        a_estimation = a.noisy48_estimate(noises1[i], n1);
        b_estimation = b.noisy48_estimate(noises2[i], n2);
        c_estimation = c.noisy48_estimate(noises3[i], n3);
        double estimation = a_estimation + b_estimation - c_estimation;
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors.push_back(error);
        error_union = pow(((c_estimation - (double)n3) / (double)n3), 2.0);
        errors_union.push_back(error_union);
    }
    error_union = 0.0;
    for (int i = 0; i < errors_union.size(); ++i)
    {
        error_union += errors_union[i];
    }
    error_union = pow(error_union / (double)errors_union.size(), 0.5);
    union_inter[0] = error_union;
    error = 0.0;
    for (int i = 0; i < errors.size(); ++i)
    {
        error += errors[i];
    }
    error = pow(error / (double)errors.size(), 0.5);
    union_inter[1] = error;
    a.clear_sketch();
    b.clear_sketch();
    return error;
}
double batch_ll_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m)
{
    unsigned seed = clock(), M = m, w = 32, number_of_threads = 24;
    srand(seed);
    seed = rand();
    LLSketch a = LLSketch(M, w, seed, number_of_threads);
    LLSketch b = LLSketch(M, w, seed, number_of_threads);
    vector<double> errors;
    double error;
    for (int i = 0; i < repeat_times; ++i)
    {
        seed = rand();
        srand(seed);
        seed = rand();
        a.change_seed(seed);
        b.change_seed(seed);
        a.generate_sketch(n1, 0);
        b.generate_sketch(n2, n1 - n4);
        LLSketch c(a, b);
        double a_est, b_est, c_est;
        if (n1 < 5 * M / 2)
            a_est = a.s_estimate();
        else
        {
            a_est = a.estimate();
        }
        if (n2 < 5 * M / 2)
            b_est = b.s_estimate();
        else
        {
            b_est = b.estimate();
        }
        if (n1 + n2 - n4 < 5 * M / 2)
            c_est = c.s_estimate();
        else
        {
            c_est = c.estimate();
        }
        double estimation = a_est + b_est - c_est;
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors.push_back(error);
    }
    error = 0.0;
    for (int i = 0; i < errors.size(); ++i)
    {
        error += errors[i];
    }
    error = pow(error / (double)errors.size(), 0.5);
    a.clear_sketch();
    b.clear_sketch();
    return error;
}

double batch_noisy_ll_est(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, int noisy_type, double epsilon, double union_inter[])
{
    unsigned seed = clock(), M = m, w = 32, number_of_threads = 24;
    srand(seed);
    double real_w[3];
    unsigned long long n3 = n1 + n2 - n4;
    default_random_engine gen(seed);
    normal_distribution<double> dis(0.0, 39.0);
    LLSketch a = LLSketch(M, w, seed, number_of_threads);
    LLSketch b = LLSketch(M, w, seed, number_of_threads);
    vector<double> errors;
    vector<double> errors_union;
    real_w[0] = ceil(log2((double)(n1) / (double)(M))) + 6;
    real_w[1] = ceil(log2((double)(n2) / (double)(M))) + 6;
    real_w[2] = ceil(log2((double)(n3) / (double)(M))) + 6;
    std::vector<double> noises1(repeat_times, 0);
    std::vector<double> noises2(repeat_times, 0);
    std::vector<double> noises3(repeat_times, 0);
    if (noisy_type == 1)
    {
        noises1 = generate_gaussian_noise(repeat_times, epsilon, real_w[0]);
        noises2 = generate_gaussian_noise(repeat_times, epsilon, real_w[1]);
        noises3 = generate_gaussian_noise(repeat_times, epsilon, real_w[2]);
    }
    else if (noisy_type == 2)
    {
        noises1 = generate_cdp_noise(repeat_times, epsilon, real_w[0]);
        noises2 = generate_cdp_noise(repeat_times, epsilon, real_w[1]);
        noises3 = generate_cdp_noise(repeat_times, epsilon, real_w[2]);
    }
    else if (noisy_type == 3)
    {
        noises1 = generate_binomial_noise(repeat_times, epsilon, real_w[0]);
        noises2 = generate_binomial_noise(repeat_times, epsilon, real_w[1]);
        noises3 = generate_binomial_noise(repeat_times, epsilon, real_w[2]);
    }
    double error, error_union;
    for (int i = 0; i < repeat_times; ++i)
    {
        seed = rand();
        srand(seed);
        seed = rand();
        a.change_seed(seed);
        b.change_seed(seed);
        a.generate_sketch(n1, 0);
        b.generate_sketch(n2, n1 - n4);
        LLSketch c(a, b);
        double a_estimation, b_estimation, c_estimation;
        a_estimation = a.noisy_estimate(noises1[i], n1);
        b_estimation = b.noisy_estimate(noises2[i], n2);
        c_estimation = c.noisy_estimate(noises3[i], n3);
        error_union = pow(((c_estimation - (double)n3) / (double)n3), 2.0);
        errors_union.push_back(error_union);
        double estimation = a_estimation + b_estimation - c_estimation;
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors.push_back(error);
    }
    error_union = 0.0;
    for (int i = 0; i < errors_union.size(); ++i)
    {
        error_union += errors_union[i];
    }
    error_union = pow(error_union / (double)errors_union.size(), 0.5);
    union_inter[0] = error_union;
    error = 0.0;
    for (int i = 0; i < errors.size(); ++i)
    {
        error += errors[i];
    }
    error = pow(error / (double)errors.size(), 0.5);
    union_inter[1] = error;
    a.clear_sketch();
    b.clear_sketch();
    return error;
}
double batch_fms_intersection(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m) //  n1: car of set S1;  n2: car of set S2;  n3: car of set S1 union S2;  n4: car of set S1 inter S2
{
    unsigned i = 0;
    unsigned seed = clock(), M = 4096, w = 32, number_of_threads = 24;
    vector<double> errors, errors3;
    double error, estimation;
    FMS a = FMS(M, w, seed, number_of_threads);
    FMS b = FMS(M, w, seed, number_of_threads);
    srand(seed);
    for (i = 0; i < repeat_times; ++i)
    {
        seed = rand();
        srand(seed);
        seed = rand();
        a.change_seed(seed);
        b.change_seed(seed);
        a.generate_sketch(n1, 0);
        b.generate_sketch(n2, n1 - n4);
        FMS c(a, b);
        double c_est = c.estimate();
        estimation = a.estimate() + b.estimate() - c_est;
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors3.push_back(error);
    }
    error = 0.0;
    for (i = 0; i < errors3.size(); ++i)
    {
        error += errors3[i];
    }
    error = pow(error / (double)errors3.size(), 0.5);
    a.clear_sketch();
    b.clear_sketch();

    return error;
}
//  n1: car of set S1;  n2: car of set S2;  n3: car of set S1 union S2;  n4: car of set S1 inter S2
// noisy_type: gaussian_noise(1) CDP(2) Binomial mechanism(3)
double batch_noisy_fms_intersection(unsigned long long n1, unsigned long long n2, unsigned long long n4, unsigned repeat_times, unsigned m, double res[], int noisy_type, double epsilon)
{
    unsigned i = 0;
    unsigned seed = clock(), M = m, w = 32, number_of_threads = 24;
    unsigned long long n3 = n1 + n2 - n4;
    vector<double> errors;
    vector<double> errors3;
    vector<double> errors_union;
    double error, estimation, error_union;
    FMS a = FMS(M, w, seed, number_of_threads);
    FMS b = FMS(M, w, seed, number_of_threads);
    srand(seed);
    std::vector<double> noises(3 * repeat_times, 0);
    if (noisy_type == 1)
    {
        noises = generate_gaussian_noise(3 * repeat_times, epsilon, 1);
    }
    else if (noisy_type == 2)
    {
        noises = generate_cdp_noise(3 * repeat_times, epsilon, 1);
    }
    else if (noisy_type == 3)
    {
        noises = generate_binomial_noise(3 * repeat_times, epsilon, 1);
    }
    double noisy_error = 0;
    for (int i = 0; i < 3 * repeat_times;i++){
        noisy_error += noises[i] * noises[i];
    }
    noisy_error = pow(noisy_error / (3 * repeat_times), 0.5);
    for (i = 0; i < repeat_times; ++i)
    {
        seed = rand();
        srand(seed);
        seed = rand();
        a.change_seed(seed);
        b.change_seed(seed);
        a.generate_sketch(n1, 0);
        b.generate_sketch(n2, n1 - n4);
        FMS c(a, b);
        double a_estimation, b_estimation, c_estimation;
        a_estimation = a.noisy_estimate(noises[i], n1);
        b_estimation = b.noisy_estimate(noises[i + repeat_times], n2);
        c_estimation = c.noisy_estimate(noises[i + 2 * repeat_times], n3);
        error_union = pow(((c_estimation - (double)n3) / (double)n3), 2.0);
        errors_union.push_back(error_union);
        estimation = a_estimation + b_estimation - c_estimation;
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors3.push_back(error);
        estimation = n1 + n2 - c_estimation;
        error = pow(((estimation - (double)n4) / (double)n4), 2.0);
        errors.push_back(error);
        }
    error = 0.0;
    for (i = 0; i < errors3.size(); ++i)
    {
        error += errors3[i];
    }
    error = pow(error / (double)errors3.size(), 0.5);
    res[0] = error;
    error = 0.0;
    for (i = 0; i < errors.size(); ++i)
    {
        error += errors[i];
    }
    error = pow(error / (double)errors.size(), 0.5);
    res[1] = error;
    error_union = 0.0;
    for (i = 0; i < errors_union.size(); ++i)
    {
        error_union += errors_union[i];
    }
    error_union = pow(error_union / (double)errors_union.size(), 0.5);
    res[2] = error_union;
    // ct_errors output
    a.clear_sketch();
    b.clear_sketch();
    return error;
}
