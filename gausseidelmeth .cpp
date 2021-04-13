#include <stdio.h>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <iomanip>
#include <cstdlib>
#include <fstream>
using namespace std;

double** res;
int N, nstep, num_stream;
double eps, h;
double Sum2;
double Sum1;


typedef double(*func)(double x, double y);
double** SolveEquationDirikhle(func G, func F, int N, double eps);
double** SolveEquationDirikhlePar(func G, func F, int N, double eps);
double F(double x, double y);
double G(double x, double y);


int main(int argc, char** argv) {
	N = 200;
	eps = 0.0001;
	nstep = 10;
	num_stream =4;

	fstream file1, file2, file3, file4;
	//file1.open("file1.txt", ios::out);
	file2.open("file4.txt", ios::out);
	file3.open("diff.txt", ios::out);
	file4.open("N.txt", ios::out);

	N = 0;
	for (int k = 0; k < 8; k++) {
		N = N + 100;
		file4 << N << "\n";
		cout <<"N="<< N << endl;
		//Cons version program***************************************
		double start_posl = omp_get_wtime();
		res = SolveEquationDirikhle(&G, &F, N, eps);
		double finish_posl = omp_get_wtime();

		cout << "Time Posl = " << finish_posl - start_posl << endl;


		for (int i = 1; i < N + 1; i++) {
			for (int j = 1; j < N + 1; j++)
				Sum1 += res[i][j];
		}
		
		//parall version program***************************************

		double start_par = omp_get_wtime();
		res = SolveEquationDirikhlePar(&G, &F, N, eps);
		double finish_par = omp_get_wtime();

		for (int i = 1; i < N + 1; i++) {
			for (int j = 1; j < N + 1; j++)
				Sum2 += res[i][j];
		}
		cout << "Time Par = " << finish_par - start_par << endl;
		file2 << finish_par - start_par << "\n";
		cout << "Sum2-Sum1 = " << Sum2 - Sum1 << endl;
		file3 << Sum2 - Sum1 << "\n";
	}

	return 0;
}



double F(double x, double y)
{
	return sin(x) * sin(y);
}
double G(double x, double y)
{
	if (x == 0.) return 1. - 2. * y;
	if (x == 1.) return -1. + 2. * y;
	if (y == 0.) return 1. - 2. * x;
	if (y == 1.) return -1. + 2. * x;
}


double** SolveEquationDirikhle(func G, func F, int N, double eps)
{
	double const h = 1.0 / (N + 1);
	double** f = new double* [N];
	for (int i = 0; i < N; i++)
	{
		f[i] = new double[N];
		for (int j = 0; j < N; j++)
			f[i][j] = F((i + 1) * h, (j + 1) * h);
	}
	double** u = new double* [N + 2];
	double** u0 = new double* [N + 2];
	for (int i = 1; i < N + 1; i++)
	{
		u[i] = new double[N + 2];
		u0[i] = new double[N + 2];
		for (int j = 1; j < N + 1; j++) {
			u[i][j] = 0.;
			u0[i][j] = 0.;
		}
		u[i][0] = G(i * h, 0.);
		u[i][N + 1] = G(i * h, (N + 1) * h);
	}
	u[0] = new double[N + 2];
	u[N + 1] = new double[N + 2];
	for (int j = 0; j < N + 2; j++)
	{
		u[0][j] = G(0, j * h);
		u[N + 1][j] = G((N + 1) * h, j * h);
	}


	double max;
	double* mx = new double[N];
	int iter = 0;
	do
	{
		iter++;
		// нарастание волны (k - длина фронта волны)
		for (int k = 1; k < N + 1; k++)
		{
			mx[k] = 0;
			for (int i = 1; i < k + 1; i++)
			{
				int j = k + 1 - i;
				//cout << "i ="<<i<<"j =" << j << endl;
				double u0 = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j]
					+ u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
				double d = abs(u[i][j] - u0);
				if (d > mx[i])
					mx[i] = d;
			}
		}
		// затухание волны (k - длина фронта волны)
		for (int k = N - 1; k > 0; k--)
		{
			for (int i = N - k + 1; i < N + 1; i++)
			{
				int j = 2 * N - k - i +1;
				//cout << "i =" << i << "j =" << j << endl;
				double u0 = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j]
					+ u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
				double d = abs(u[i][j] - u0);
				if (d > mx[i])
					mx[i] = d;
			}
		}
		max = 0;
		for (int i = 1; i < N + 1; i++)
			if (mx[i] > max)
				max = mx[i];
	} while (max > eps);


	cout << "IterCnt = " << iter << endl;
	return u;
}


double** SolveEquationDirikhlePar(func G, func F, int N, double eps)
{
	h = 1.0 / (N + 1);
	double const h = 1.0 / (N + 1);
	double** f = new double* [N];
	for (int i = 0; i < N; i++)
	{
		f[i] = new double[N];
		for (int j = 0; j < N; j++)
			f[i][j] = F((i + 1) * h, (j + 1) * h);
	}
	double** u = new double* [N + 2];
	double** u0 = new double* [N + 2];
	double** u2 = new double* [N + 2];
	for (int i = 1; i < N + 1; i++)
	{
		u[i] = new double[N + 2];
		u0[i] = new double[N + 2];
		u2[i] = new double[N + 2];
		for (int j = 1; j < N + 1; j++) {
			u[i][j] = 0.;
			u0[i][j] = 0.;
			u2[i][j] = 0.;
		}

		u[i][0] = G(i * h, 0.);
		u[i][N + 1] = G(i * h, (N + 1) * h);
	}
	u[0] = new double[N + 2];
	u[N + 1] = new double[N + 2];
	for (int j = 0; j < N + 2; j++)
	{
		u[0][j] = G(0, j * h);
		u[N + 1][j] = G((N + 1) * h, j * h);
	
	}


	double max;
	double* mx = new double[N];
	int iter = 0;
	omp_set_num_threads(num_stream);
	do
	{
		iter++;
		// нарастание волны (k - длина фронта волны)
		for (int k = 1; k < N + 1; k++)
		{
			mx[k] = 0;
#pragma omp parallel for shared(u, f, mx)
			for (int i = 1; i < k + 1; i++)
			{
				int j = k + 1 - i;
				double u0 = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j]
					+ u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
				double d = abs(u[i][j] - u0);
				if (d > mx[i])
					mx[i] = d;
			}
		}
		// затухание волны (k - длина фронта волны)
		for (int k = N - 1; k > 0; k--)
		{
		#pragma omp parallel for shared(u, f, mx)
			for (int i = N - k + 1; i < N + 1; i++)
			{
				int j = 2 * N - k - i + 1;
				double u0 = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j]
					+ u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
				double d = abs(u[i][j] - u0);
				if (d > mx[i])
					mx[i] = d;
			}
		}
		max = 0;
		for (int i = 1; i < N + 1; i++)
			if (mx[i] > max) max = mx[i];
	} while (max > eps);

	cout << "IterCnt = " << iter << endl;
	return u;
}
