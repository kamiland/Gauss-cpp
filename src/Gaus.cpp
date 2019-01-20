#include <math.h>
#include "nrutil.h"
#include <iostream>
#include <iomanip>
using namespace std;

#define SWAP(a, b) {temp = (a); (a) = (b); (b) = temp; }

void gaussj(double **a,  int n,  double **b,  int m)
{
	int *indxc, *indxr, *ipiv;
	int i, icol, irow, j, k, l, ll;
	double big, dum, pivinv, temp;

	indxc = ivector(1, n);
	indxr = ivector(1, n);
	ipiv = ivector(1, n);

	for (j = 1; j <= n; j++)
		ipiv[j] = 0;

	for (i = 1; i <= n; i++)
	{
		big = 0.0;
		for (j = 1; j <= n; j++)
			if (ipiv[j] !=  1)
		for (k = 1; k <= n; k++)
		{
			if (ipiv[k] == 0)
			{
				if (fabs(a[j][k]) >=  big)
				{
					big = fabs(a[j][k]);
					irow = j;
					icol = k;
				}
			}
		}
		++(ipiv[icol]);
		 if (irow !=  icol)
		 {
			for (l = 1; l <= n; l++)
				SWAP(a[irow][l], a[icol][l])
			for (l = 1; l <= m; l++)
				SWAP(b[irow][l], b[icol][l])
		}

		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)
			nrerror("gaussj: Singular Matrix");

		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 1; l <= n; l++) a[icol][l] *=pivinv;
		for (l = 1; l <= m; l++) b[icol][l] *=pivinv;

		for (ll = 1; ll <= n; ll++)
			if (ll !=  icol)
			{
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l]*dum;
				for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l = n; l >= 1; l--)
	{
		if (indxr[l] !=  indxc[l])
			for (k = 1; k <= n; k++)
				SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	}

	free_ivector(ipiv, 1, n);
	free_ivector(indxr, 1, n);
	free_ivector(indxc, 1, n);

}

void printMatrix(double **A, int I, int J, int superPrecision = 0)
{
	for(int i = 1; i <= I; i++)
	{
		for(int j = 1; j <= J; j++)
		{
			if(abs(A[i][j]) < 0.001 && !superPrecision)
			{
				cout << setprecision(3) << "0\t";
			}
			else
			{
				cout << setprecision(3) << A[i][j] << "\t";
			}
			if(!(j % J))
				cout << endl;
		}
	}
	cout << endl;
}

void printMatrix(double *A, int I) //vector
{
	for(int i = 1; i <= I; i++)
	{
		cout << setprecision(3) << A[i] << endl;
	}
	cout << endl;
}

void printMatlab(double **A, int I, int J, int superPrecision = 0)
{
	cout << "[";
	for(int i = 1; i <= I; i++)
	{
		for(int j = 1; j <= J; j++)
		{
			if(abs(A[i][j]) < 0.000001 && !superPrecision)
			{
				cout << setprecision(5) << " 0";
			}
			else
			{
				cout << setprecision(5) << " " << A[i][j];
			}

			if(!(j % J) && i < I)
				cout << ";";
		}
	}
		cout << "]\n" << endl;
}


#define N 10
double **matrixMultiply(double **A, int nA, int mA, double **B, int nB, int mB)
{
	double **C;

	C = dmatrix(1, N, 1, N);
	C[nA][nB] = {0};

	for(int i = 1; i <= nA; i++)
	{
		for(int j = 1; j <= mB; j++)
		{
			for(int k = 1; k <= mA; k++)
			{
				C[i][j] += A[i][k] * B[k][j];

			}
		}
	}

	return C;
}


int main(void)
{
	double **A, **oldA, **check, **B, **MultipResult;

	A = dmatrix(1, N, 1, N);
	oldA = dmatrix(1, N, 1, N);
	B = dmatrix(1, N, 1, N);

	for(int i = 1; i <= N; i++)
	{
		for(int j = 1; j <= N; j++)
		{
			A[i][j] =(double) (i * j) / (i + j);
			oldA[i][j] = A[i][j];
			B[i][j] =(double) (i - j) / (i * j);
		}
	}
	printMatrix(A, N, N);
	printMatrix(B, N, N);
	printMatlab(A, N, N);
	printMatlab(B, N, N);

	gaussj(A, N, B, N);

	cout << "after gauss" << endl;
	printMatlab(A, N, N);
	printMatlab(B, N, N);

	MultipResult = matrixMultiply(A, N, N, B, N, N);

	check = matrixMultiply(A,N,N,oldA,N,N);

	printMatlab(check, N, N);
	printMatrix(check, N, N);


	return 0;
}












