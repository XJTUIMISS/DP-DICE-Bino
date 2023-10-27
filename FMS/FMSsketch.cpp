#include "FMSsketch.h"


//  hyper parameters in this file
unsigned max_threads = 24;
unsigned max_M = 16384;
unsigned max_w = 40;

FMS::FMS(unsigned t_M, unsigned t_w, unsigned t_seed, unsigned t_number_of_threads)
{
	M = t_M;
	r = log(M) / log(2) + 0.1;
	w = t_w;
	seed = t_seed;
	number_of_threads = t_number_of_threads;
	sk = new unsigned[M * w];
	temp_sk = new unsigned[number_of_threads * M * w];
	threads = new thread[number_of_threads];
}


FMS::FMS(FMS & sk1, FMS & sk2)
{
	if (sk1.M != sk2.M || sk1.w != sk2.w || sk1.seed != sk2.seed)
	{
		cout << "cannot merge!" << endl;
		return;
	}
	unsigned i, j;
	M = sk1.M;
	w = sk1.w;
	seed = sk1.seed;
	r = sk1.r;
	number_of_threads = (sk1.number_of_threads > sk2.number_of_threads) ? sk1.number_of_threads : sk2.number_of_threads;
	sk = new unsigned[M * w];
	temp_sk = new unsigned[number_of_threads * M * w];
	threads = new thread[number_of_threads];
	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < w; ++j)
		{
			sk[i * w + j] = sk1.sk[i * w + j] | sk2.sk[i * w + j];
		}
	}
}


FMS::~FMS()
{
	delete [] sk;
	delete [] temp_sk;
	delete [] threads;
}


void FMS::clear_sketch()
{
	memset(sk, 0, M * w * sizeof(sk[0]));
	memset(temp_sk, 0, number_of_threads * M * w * sizeof(temp_sk[0]));
}

void FMS::thread_gene_sketch(unsigned index)
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
		MurmurHash3_x64_128(&i, 8, seed, hash_uni);
		x = hash_uni[0] >> (64 - r);
		y = rank_bits(hash_uni[1]);
		temp_sk[index * M * w + x * w + y] = 1;
	}
}


void FMS::thread_sketch_combine()
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


unsigned long long FMS::rank_bits(unsigned long long bits)
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


double FMS::generate_sketch(unsigned long long t_n, unsigned long long t_n_base)
{
	n = t_n;
	the_n_base = t_n_base;
	clear_sketch();
	double time_consumption = 0.0;   
	unsigned long long i = 0, x = 0, y = 0;
	unsigned long long hash_uni[2] = {0};
	if (n < pow(10, 6) - 10)
	{
		for (i = t_n_base; i < n + t_n_base; ++i)
		{
			MurmurHash3_x64_128(&i, 8, seed, hash_uni);
			x = hash_uni[0] >> (64 - r);
			y = rank_bits(hash_uni[1]);
			sk[x * w + y] = 1;
		}
	}
	else  //  multi thread from 10^6 - 10
	{
		for (i = 0; i < number_of_threads; ++i)
		{
			threads[i] = thread(&FMS::thread_gene_sketch, this, i);
		}
		for (i = 0; i < number_of_threads; ++i)
		{
			threads[i].join();
		}
		thread_sketch_combine();
	}
	return time_consumption;
}


double FMS::px(double x)
{
	if (x < w - 1.5)
	{
		return 1.0 / ((double)M * pow(2.0, x + 1.0));
	}
	else
	{
		return 1.0 / ((double)M * pow(2.0, x));
	}
	
}


unsigned FMS::count_zero()
{
	unsigned i, j, count = M * w;
	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < w; ++j)
		{
			count -= sk[i * w + j];
		}
	}
	return count;
}


double FMS::func_ev(double t_n)
{
	unsigned i = 0;
	double re = 0.0;
	for (i = 0; i < w; ++i)
	{
		re += pow(1.0 - px((double)i), t_n);
	}
	return re / (double)w;
}


double FMS::d_func_ev(double t_n)
{
	unsigned i = 0;
	double re = 0.0;
	for (i = 0; i < w; ++i)
	{
		re += pow(1.0 - px((double)i), t_n) * log(1.0 - px((double)i));
	}
	return re / (double)w;
}

unsigned FMS::noisy_count_zero()
{
	unsigned i, j, count = M * real_w;
	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < real_w; ++j)
		{
			count -= sk[i * w + j];
		}
	}
	return count;
}


double FMS::noisy_func_ev(double t_n)
{
	unsigned i = 0;
	double re = 0.0;
	for (i = 0; i < real_w; ++i)
	{
		re += pow(1.0 - px((double)i), t_n);
	}
	return re / (double)real_w;
}


double FMS::noisy_d_func_ev(double t_n)
{
	unsigned i = 0;
	double re = 0.0;
	for (i = 0; i < real_w; ++i)
	{
		re += pow(1.0 - px((double)i), t_n) * log(1.0 - px((double)i));
	}
	return re / (double)real_w;
}

double FMS::estimate()
{
	double v, t_n;
	unsigned i = 0;
	v = (double)count_zero() / (double)(M * w);
	t_n = 10.0;
	while((func_ev(t_n) - v) * (func_ev(t_n + 0.5) - v) > 0)
	{
		t_n = t_n - (func_ev(t_n) - v) / d_func_ev(t_n);
	}
	return t_n;
}

double FMS::noisy_estimate(double noise, double real_n)
{
	double v, t_n;
	real_w = ceil(log2((double)(real_n) / (double)(M))) + 6;
	unsigned i = 0;
	v = ((double)noisy_count_zero() + noise) / (double)(M * real_w);
	if(v>1){
		v = 1;
	}
	else if(v<0){
		v = 0;
	}
	t_n = 10.0;
	while((noisy_func_ev(t_n) - v) * (noisy_func_ev(t_n + 0.5) - v) > 0)
	{
		t_n = t_n - (noisy_func_ev(t_n) - v) / noisy_d_func_ev(t_n);
	}
	return t_n;
}


void FMS::print_sketch(unsigned t_lines)
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


void FMS::change_seed(unsigned t_seed)
{
	seed = t_seed;
}
