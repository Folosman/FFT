#include <complex.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

#define M_PI 3.14159265358979323846
#define N 128
#define F0 25

typedef std::vector<std::complex<double>> compexVector ;

void FFT(compexVector &analysValue)
{
    int n = analysValue.size();
    if (n == 1) return;

    compexVector W(n);
    double alpha;

    for (int i = 0; i < n; i ++)
    {
            alpha = 2 * M_PI * i / n;
            W[i] = std::complex<double>(cos(alpha), sin(alpha));
    }

    compexVector A(n / 2), B(n / 2);
    for (int i = 0; i < n / 2; i++)
    {
        A[i] = analysValue[i * 2];
        B[i] = analysValue[i * 2 + 1];
    }
    FFT(A);
    FFT(B);
    for ( int i = 0; i < n; i++)
    {
        analysValue[i] = A[i % (n / 2)] + W[i] * B[i % (n / 2)];
    }
}

int main()
{
    double w = F0;

    compexVector harmonic;
    compexVector output;
    double tempVar;

    for (int i = 0; i < N; i++)
    {
        tempVar = sin(w * i);
        std::cout << tempVar << " ";
        harmonic.push_back(tempVar);
    }

    std::cout << "\n\n";

    FFT(harmonic);
    for (int i = 0; i < N; i++)
    {
        std::cout <<" " << harmonic[i] << "; ";
    }
    return 0;
}
