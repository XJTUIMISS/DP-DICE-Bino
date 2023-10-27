#include "batch_est.h"
#include "build/FMSConfig.h"
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

using namespace std;
using namespace chrono;

// calculate the RMSD under different M without noisy
void cal_err_M(unsigned long long n1, unsigned long long n2, unsigned long long n4, bool LL, bool FM)
{
	int begin = clock();
	double error = 0.0;
	double res[3];
	double union_inter[2];
	double fm_union_inter[2];
	FILE *fp;
	fp = fopen("fms_m.txt", "w");
	fflush(fp);
	if (LL == 1)
	{
		fprintf(fp, "LLsketch: \n");
		for (int i = 1; i <= 4; i++)
		{
			error = batch_ll_est(n1, n2, n4, 100, 1024 * pow(2, i - 1), union_inter);
			fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
			fflush(fp);
		}
		fprintf(fp, "\n");
	}
	if (FM == 1)
	{
		fprintf(fp, "FMsketch: \n");
		fflush(fp);
		for (int i = 1; i <= 4; i++)
		{
			error = batch_fm_est(n1, n2, n4, 100, 1024 * pow(2, i - 1), union_inter);
			fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
			fflush(fp);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "FMSsketch: \n");
	fflush(fp);
	for (int i = 1; i <= 4; i++)
	{
		error = batch_fms_est(n1, n2, n4, 100, 1024 * pow(2, i - 1), union_inter);
		fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
		fflush(fp);
	}
	fprintf(fp, "\n");
	fclose(fp);
	int end = clock();
	cout << "time : " << end - begin << "ms" << endl;
}

// calculate the RMSD under different M with noisy
// noisy_type: gaussian_noise(1) CDP(2) Binomial mechanism(3)
void cal_err_noisy_M(unsigned long long n1, unsigned long long n2, unsigned long long n4, bool LL, bool FM)
{
	int begin = clock();
	double error = 0.0;
	double res[3];
	double union_inter[2];
	double fm_union_inter[2];
	FILE *fp;
	fp = fopen("fms_noisy_m.txt", "w");
	fflush(fp);
	for (int noisy_type = 1; noisy_type <= 3; noisy_type++)
	{
		if (noisy_type == 1)
		{
			fprintf(fp, "noisy_type: Gaussian mechanism\n");
		}
		else if (noisy_type == 2)
		{
			fprintf(fp, "noisy_type: Distributed Discrete Gaussian mechanism\n");
		}
		else if (noisy_type == 3)
		{
			fprintf(fp, "noisy_type: MPCDP-CP\n");
		}
		if (LL == 1)
		{
			fprintf(fp, "LLsketch: \n");
			for (int i = 1; i <= 4; i++)
			{
				error = batch_noisy_ll_est(n1, n2, n4, 100, 1024 * pow(2, i - 1), noisy_type, 0.1, union_inter);
				fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
				fflush(fp);
			}
			fprintf(fp, "\n");
		}
		if (FM == 1)
		{
			fprintf(fp, "FMsketch: \n");
			fflush(fp);
			for (int i = 1; i <= 4; i++)
			{
				error = batch_noisy_fm_est(n1, n2, n4, 100, 1024 * pow(2, i - 1), noisy_type, 0.1, union_inter);
				fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
				fflush(fp);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "FMSsketch: \n");
		fflush(fp);
		for (int i = 1; i <= 4; i++)
		{
			error = batch_noisy_fms_est(n1, n2, n4, 100, 1024 * pow(2, i - 1), res, noisy_type, 0.1);
			fprintf(fp, "%lf, %lf, \n", res[2], res[0]);
			fflush(fp);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	int end = clock();
	cout << "time : " << end - begin << "ms" << endl;
}

// calculate the RMSD under different epsilon with noisy
// noisy_type: gaussian_noise(1) CDP(2) Binomial mechanism(3)
void cal_err_noisy_epsilon(unsigned long long n1, unsigned long long n2, unsigned long long n4, bool LL, bool FM)
{
	int begin = clock();
	double error = 0.0;
	double res[3];
	double union_inter[2];
	double fm_union_inter[2];
	FILE *fp;
	fp = fopen("fms_noisy_epsilon.txt", "w");
	fflush(fp);
	for (int noisy_type = 1; noisy_type <= 3; noisy_type++)
	{
		if (noisy_type == 1)
		{
			fprintf(fp, "noisy_type: Gaussian mechanism\n");
		}
		else if (noisy_type == 2)
		{
			fprintf(fp, "noisy_type: Distributed Discrete Gaussian mechanism\n");
		}
		else if (noisy_type == 3)
		{
			fprintf(fp, "noisy_type: MPCDP-CP\n");
		}
		if (LL == 1)
		{
			fprintf(fp, "LLsketch: \n");
			for (int i = 1; i <= 5; i++)
			{
				error = batch_noisy_ll_est(n1, n2, n4, 100, 4096, noisy_type, 0.1 * i, union_inter);
				fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
				fflush(fp);
			}
			fprintf(fp, "\n");
		}
		if (FM == 1)
		{
			fprintf(fp, "FMsketch: \n");
			fflush(fp);
			for (int i = 1; i <= 5; i++)
			{
				error = batch_noisy_fm_est(n1, n2, n4, 100, 4096, noisy_type, 0.1 * i, union_inter);
				fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
				fflush(fp);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "FMSsketch: \n");
		fflush(fp);
		for (int i = 1; i <= 5; i++)
		{
			error = batch_noisy_fms_est(n1, n2, n4, 100, 4096, res, noisy_type, 0.1 * i);
			fprintf(fp, "%lf, %lf, \n", res[2], res[0]);
			fflush(fp);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	int end = clock();
	cout << "time : " << end - begin << "ms" << endl;
}

int main(int argc, char *argv[])
{
	cout << " Version " << FMS_VERSION_MAJOR << "."
		 << FMS_VERSION_MINOR << std::endl;
	unsigned long long n[3];
	bool LL = 0, FM = 0;
	n[0] = 100'000;
	n[1] = 100'000;
	n[2] = 50'000;
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "c") == 0 || strcmp(argv[i], "C") == 0)
		{
			for (int j = 0; j < 3; j++)
			{
				i++;
				n[j] = atoi(argv[i]);
			}
		}
		if (strcmp(argv[i], "ll") == 0 || strcmp(argv[i], "LL") == 0)
		{
			LL = 1;
		}
		if (strcmp(argv[i], "fm") == 0 || strcmp(argv[i], "FM") == 0)
		{
			FM = 1;
		}
	}
	cal_err_M(n[0], n[1], n[2], LL, FM);
	cal_err_noisy_M(n[0], n[1], n[2], LL, FM);
	cal_err_noisy_epsilon(n[0], n[1], n[2], LL, FM);
	return 0;
}
