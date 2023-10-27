#include "LLSketch.h"

LLSketch::LLSketch(unsigned t_M, unsigned t_w, unsigned t_seed, unsigned t_number_of_threads)
{
	M = t_M;
	r = log(M) / log(2) + 0.1;
	w = t_w;
	seed = t_seed;
	number_of_threads = t_number_of_threads;
	sk = new unsigned[M];
	temp_sk = new unsigned[number_of_threads * M];
	threads = new thread[number_of_threads];
}

LLSketch::LLSketch(LLSketch &sk1, LLSketch &sk2)
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
	sk = new unsigned[M];
	temp_sk = new unsigned[number_of_threads * M];
	threads = new thread[number_of_threads];
	for (i = 0; i < M; ++i)
	{
		sk[i] = (sk1.sk[i] > sk2.sk[i]) ? sk1.sk[i] : sk2.sk[i];
	}
}

LLSketch::~LLSketch()
{
	delete[] sk;
	delete[] temp_sk;
	delete[] threads;
}

void LLSketch::clear_sketch()
{
	memset(sk, 0, M * sizeof(sk[0]));
	memset(temp_sk, 0, number_of_threads * M * sizeof(temp_sk[0]));
}

void LLSketch::thread_gene_sketch(unsigned index)
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
		y = (y > w) ? w : y;
		if (temp_sk[index * M + x] < y)
			temp_sk[index * M + x] = y;
	}
}

void LLSketch::thread_sketch_combine()
{
	unsigned long long x, y, i;
	for (i = 0; i < number_of_threads; ++i)
	{
		for (x = 0; x < M; ++x)
		{
			sk[x] = (sk[x] > temp_sk[i * M + x]) ? sk[x] : temp_sk[i * M + x];
		}
	}
}

unsigned long long LLSketch::rank_bits(unsigned long long bits)
{
	unsigned i = 0;
	unsigned long long musk = 1;
	for (i = 0; i < w - 1; ++i)
	{
		if ((musk & bits) != 0)
		{
			return i + 1;
		}
		bits = bits >> 1;
	}
	return w;
}

double LLSketch::generate_sketch(unsigned long long t_n, unsigned long long t_n_base)
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
			if (y > w)
			{
				y = w;
				cout << "w-1 limitation" << endl;
			}
			if (sk[x] < y)
			{
				sk[x] = y;
			}
		}
	}
	else //  multi thread from 10^6 - 10
	{
		for (i = 0; i < number_of_threads; ++i)
		{
			threads[i] = thread(&LLSketch::thread_gene_sketch, this, i);
		}
		for (i = 0; i < number_of_threads; ++i)
		{
			threads[i].join();
		}
		thread_sketch_combine();
	}
	return time_consumption;
}

double LLSketch::s_estimate()
{
	double count = 0, t_n;
	for (int i = 0; i < M; i++)
	{
		if (sk[i] == 0)
			count++;
	}
	if (count == 0)
	{
		return estimate();
	}
	else
	{
		t_n = M * log((double)M / count);
	}
	return t_n;
}

double LLSketch::estimate()
{
	double Arpha = 0.39701, average_r = 0, t_n;
	for (int i = 0; i < M; i++)
	{
		average_r += sk[i];
	}
	average_r = average_r / M;
	t_n = Arpha * (double)M * (pow(2, average_r));
	return t_n;
}

double LLSketch::s_noisy_estimate(double noise, double real_n)
{
	double count = 0, t_n;
	for (int i = 0; i < M; i++)
	{
		if (sk[i] == 0)
			count++;
	}
	if (count == 0)
	{
		return noisy_estimate(noise, real_n);
	}
	else
	{
		count += noise;
		if (count < 1)
		{
			cout << "low limit" << endl;
			count = 1;
		}
		else if (count > M)
		{
			cout << "high limit" << endl;
			count = M;
		}
		t_n = M * log((double)M / count);
	}
	return t_n;
}

double LLSketch::noisy_estimate(double noise, double real_n)
{
	double Arpha = 0.39701, average_r = 0, t_n, real_w = ceil(log2((double)(real_n)/(double)(M))) + 6;
	for (int i = 0; i < M; i++)
	{
		if (sk[i] > real_w)
		{
			average_r += real_w;
		}
		else
		{
			average_r += sk[i];
		}
	}
	average_r = (average_r + noise) / M;
	if (average_r < 0)
	{
		//cout << "low limit" << endl;
		average_r = 0;
	}
	else if (average_r > real_w)
	{
		//cout << "high limit" << endl;
		average_r = real_w;
	}
	t_n = Arpha * M * (pow(2, average_r));
	return t_n;
}

void LLSketch::print_sketch(unsigned t_lines)
{
	unsigned i = 0, j = 0;
	double R = 0;
	for (i = 0; i < M; i++)
	{
		R += sk[i];
	}
	cout << "R: " << R / M << endl;
	for (i = 0; i < t_lines; ++i)
	{
		cout << sk[i] << " ";
		cout << endl;
	}
	cout << endl;
}

void LLSketch::change_seed(unsigned t_seed)
{
	seed = t_seed;
}
