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
void cal_err_M(unsigned long long n1, unsigned long long n2, unsigned long long n4)
{
	int begin = clock();
	double error = 0.0;
	double res[3];
	double union_inter[2];
	double fm_union_inter[2];
	FILE *fp;
	fp = fopen("fms_m.txt", "w");
	fflush(fp);
	for (int i = 1; i <= 4; i++)
	{
		error = batch_ll_est(100'000, 100'000, 50'000, 100, 1024 * pow(2, i - 1));
		fprintf(fp, "%lf, \n", error);
		fflush(fp);
	}
	fprintf(fp, "\n");
	fprintf(fp, "FMSketch: \n");
	fflush(fp);
	for (int i = 1; i <= 4; i++)
	{
		error = batch_fm_est(100'000, 100'000, 50'000, 100, 1024 * pow(2, i - 1));
		fprintf(fp, "%lf, \n", error);
		fflush(fp);
	}
	fprintf(fp, "\n");
	fprintf(fp, "FMSSketch(3+1): \n");
	fflush(fp);
	for (int i = 1; i <= 4; i++)
	{
		error = batch_fms_intersection(100'000, 100'000, 50'000, 100, 1024 * pow(2, i - 1));
		fprintf(fp, "%lf, \n", error);
		fflush(fp);
	}
	fprintf(fp, "\n");
	fclose(fp);
	int end = clock();
	cout << "time : " << end - begin << "ms" << endl;
}

// calculate the RMSD under different M with noisy
// noisy_type: gaussian_noise(1) CDP(2) Binomial mechanism(3)
void cal_err_noisy_M(unsigned long long n1, unsigned long long n2, unsigned long long n4)
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
		fprintf(fp, "noisy_type: %d\n", noisy_type);
		fprintf(fp, "LLSketch: \n");
		for (int i = 1; i <= 4; i++)
		{
			error = batch_noisy_ll_est(100'000, 100'000, 50'000, 100, 1024 * pow(2, i - 1), noisy_type, 0.1, union_inter);
			fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
			fflush(fp);
		}
		fprintf(fp, "\n");
		fprintf(fp, "FMSketch: \n");
		fflush(fp);
		for (int i = 1; i <= 4; i++)
		{
			error = batch_noisy_fm_est(100'000, 100'000, 50'000, 100, 1024 * pow(2, i - 1), noisy_type, 0.1, union_inter);
			fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
			fflush(fp);
		}
		fprintf(fp, "\n");
		fprintf(fp, "FMSSketch(3+1): \n");
		fflush(fp);
		for (int i = 1; i <= 4; i++)
		{
			error = batch_noisy_fms_intersection(100'000, 100'000, 50'000, 100, 1024 * pow(2, i - 1), res, noisy_type, 0.1);
			fprintf(fp, "%lf, %lf,  %lf, \n", res[0], res[1], res[2]);
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
void cal_err_noisy_epsilon(unsigned long long n1, unsigned long long n2, unsigned long long n4)
{
	int begin = clock();
	double error = 0.0;
	double res[3];
	double union_inter[2];
	double fm_union_inter[2];
	FILE *fp;
	fp = fopen("fms_noisy_epsilon.txt", "w");
	fflush(fp);
	for (int noisy_type = 3; noisy_type <= 3; noisy_type++)
	{
		fprintf(fp, "noisy_type: %d\n", noisy_type);
		fprintf(fp, "LLSketch: \n");
		for (int i = 1; i <= 5; i++)
		{
			error = batch_noisy_ll_est(100'000, 100'000, 50'000, 100, 4096, noisy_type, 0.1 * i, union_inter);
			fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
			fflush(fp);
		}
		fprintf(fp, "\n");
		fprintf(fp, "FMSketch: \n");
		fflush(fp);
		for (int i = 1; i <= 5; i++)
		{
			error = batch_noisy_fm_est(100'000, 100'000, 50'000, 100, 4096, noisy_type, 0.1 * i, union_inter);
			fprintf(fp, "%lf, %lf, \n", union_inter[0], union_inter[1]);
			fflush(fp);
		}
		fprintf(fp, "\n");
		fprintf(fp, "FMSSketch(3+1): \n");
		fflush(fp);
		for (int i = 1; i <= 5; i++)
		{
			error = batch_noisy_fms_intersection(100'000, 100'000, 50'000, 100, 4096, res, noisy_type, 0.1 * i);
			fprintf(fp, "%lf, %lf,  %lf, \n", res[0], res[1], res[2]);
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
		 << FMS_VERSION_MINOR << ": noise" << std::endl;
	unsigned long long n[3];
	if (argc <= 1)
	{
		n[0] = 100'000;
		n[1] = 100'000;
		n[2] = 50'000;
	}
	else
	{
		for (int i = 1; i < argc; i++)
		{
			n[i - 1] = atoi(argv[i]);
		}
	}
	cal_err_M(n[0], n[1], n[2]);
	cal_err_noisy_M(n[0], n[1], n[2]);
	cal_err_noisy_epsilon(n[0], n[1], n[2]);
	return 0;
}
