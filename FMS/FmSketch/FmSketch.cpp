#include "FmSketch.h"

FmSketch::FmSketch(unsigned t_M, unsigned t_w, unsigned t_seeds[], unsigned t_number_of_threads)
{
    M = t_M;
    r = log(M) / log(2) + 0.1;
    w = t_w;
    number_of_threads = t_number_of_threads;
    sk = new unsigned[M * w];
    temp_sk = new unsigned[number_of_threads * M * w];
    threads = new thread[number_of_threads];
    R = new unsigned[M];
    seeds = new unsigned[M];
    copy(t_seeds, t_seeds + M, seeds);
}

FmSketch::FmSketch(FmSketch &sk1, FmSketch &sk2)
{
    if (sk1.M != sk2.M || sk1.w != sk2.w || sk1.seeds[0] != sk2.seeds[0])
    {
        cout << "cannot merge!" << endl;
        return;
    }
    unsigned i, j;
    M = sk1.M;
    w = sk1.w;
    seeds = new unsigned[M];
    copy(sk1.seeds, sk1.seeds + M, seeds);
    r = sk1.r;
    number_of_threads = (sk1.number_of_threads > sk2.number_of_threads) ? sk1.number_of_threads : sk2.number_of_threads;
    sk = new unsigned[M * w];
    temp_sk = new unsigned[number_of_threads * M * w];
    threads = new thread[number_of_threads];
    R = new unsigned[M];
    for (i = 0; i < M; ++i)
    {
        for (j = 0; j < w; ++j)
        {
            sk[i * w + j] = sk1.sk[i * w + j] | sk2.sk[i * w + j];
        }
    }
}

FmSketch::~FmSketch()
{
    delete[] sk;
    delete[] temp_sk;
    delete[] threads;
    delete[] R;
    delete[] seeds;
}

void FmSketch::clear_sketch()
{
    memset(sk, 0, M * w * sizeof(sk[0]));
    memset(temp_sk, 0, number_of_threads * M * w * sizeof(temp_sk[0]));
    memset(R, 0, M * sizeof(R[0]));
}

void FmSketch::thread_gene_sketch(unsigned index)
{
    unsigned long long n_base, t_n, i, x, y;
    unsigned long long hash_uni[2] = {0};
    n_base = the_n_base + (unsigned long long)index * (n / (unsigned long long)number_of_threads);
    if (index == number_of_threads - 1)
    {
        t_n = (n % (unsigned long long)number_of_threads) + (n / (unsigned long long)number_of_threads);
    }
    else
    {
        t_n = (n / (unsigned long long)number_of_threads);
    }
    for (i = n_base; i < n_base + t_n; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            MurmurHash3_x64_128(&i, 8, seeds[j], hash_uni);
            y = rank_bits(hash_uni[1]);
            temp_sk[index * M * w + j * w + y] = 1;
        }
    }
}

void FmSketch::thread_sketch_combine()
{
    unsigned long long x, y, i;
    for (i = 0; i < number_of_threads; ++i)
    {
        for (x = 0; x < M; ++x)
        {
            for (y = 0; y < w; ++y)
            {
                sk[x * w + y] |= temp_sk[i * M * w + x * w + y];
            }
        }
    }
}

unsigned long long FmSketch::rank_bits(unsigned long long bits)
{
    unsigned i = 0;
    unsigned long long musk = 1;
    for (i = 0; i < w - 1; ++i)
    {
        if ((musk & bits) != 0)
        {
            return i;
        }
        bits = bits >> 1;
    }
    return w - 1;
}

double FmSketch::generate_sketch(unsigned long long t_n, unsigned long long t_n_base)
{
    n = t_n;
    the_n_base = t_n_base;
    clear_sketch();
    double time_consumption = 0.0;
    unsigned long long i = 0, x = 0, y = 0;
    unsigned long long hash_uni[2] = {0};
    if (n < pow(10, 4) - 10)
    {
        for (i = t_n_base; i < n + t_n_base; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                MurmurHash3_x64_128(&i, 8, seeds[j], hash_uni);
                y = rank_bits(hash_uni[1]);
                if (y > w - 1)
                {
                    y = 31;
                }
                sk[j * w + y] = 1;
            }
        }
    }
    else //  multi thread from 10^6 - 10
    {
        for (i = 0; i < number_of_threads; ++i)
        {
            threads[i] = thread(&FmSketch::thread_gene_sketch, this, i);
        }
        for (i = 0; i < number_of_threads; ++i)
        {
            threads[i].join();
        }
        thread_sketch_combine();
    }
    return time_consumption;
}

void FmSketch::generate_R()
{
    unsigned i, j;
    for (i = 0; i < M; ++i)
    {
        for (j = 0; j < w; ++j)
        {
            if (sk[i * w + j] == 0)
            {
                R[i] = j;
                break;
            }
            if (j == w - 1)
            {
                cout << "sk full!" << endl;
                R[i] = w-1;
            }
        }
    }
}

double FmSketch::estimate()
{
    generate_R();
    double Phi = 0.77351, average_r = 0, t_n;
    for (int i = 0; i < M; i++)
    {
        average_r += R[i];
    }
    average_r = average_r / M;
    t_n = 1.0 / Phi * (pow(2, average_r));
    return t_n;
}
double FmSketch::noisy_estimate(double noise, double real_n)
{
    generate_R();
    double Phi = 0.77351, average_r = 0, t_n, z_real;
    z_real = log2(real_n * Phi);
    for (int i = 0; i < M; i++)
    {
        average_r += R[i];
    }
    average_r += noise;
    average_r = (average_r + noise) / M;
    if (average_r > ceil(log2(real_n)) + 6 - 1)
    {
        average_r = ceil(log2(real_n)) + 6 - 1;
    }
    else if (average_r < 0)
    {
        average_r = 0;
    }
    t_n = 1.0 / Phi * (pow(2, average_r));
    return t_n;
}

double FmSketch::noisy48_estimate(double noise, double real_n)
{
    double Phi = 0.77351, average_r = 0, t_n, z_real, z_noisy ;
    z_real = log2(real_n * Phi);
    z_noisy = z_real + (noise / M);
    average_r += z_noisy;
    if (average_r > ceil(log2(real_n)) + 6 - 1)
    {
        average_r = ceil(log2(real_n)) + 6 - 1;
    }
    else if (average_r < 0)
    {
        average_r = 0;
    }
    t_n = 1.0 / Phi * (pow(2, average_r));
    return t_n;
}
void FmSketch::print_sketch(unsigned t_lines)
{
    unsigned i = 0, j = 0;
    for (i = 0; i < t_lines; ++i)
    {
        for (j = 0; j < w; ++j)
        {
            cout << sk[i * w + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void FmSketch::change_seed(unsigned t_seeds[])
{
    copy(t_seeds, t_seeds + M, seeds);
}
