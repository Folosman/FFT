#include <complex.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>

#define M_PI 3.14159265358979323846

typedef std::vector<std::complex<double>> compexVector ;

compexVector FFT(const compexVector &analysValue)
{
    int n = analysValue.size();
    std::cout << "TESTrecurs";
    if (n == 1) return compexVector(1, analysValue[0]);

    compexVector W(n);
    double alpha;

    for (int i = 0; i < n; i ++)
    {
        alpha = 2 * M_PI * i / n;
        W[i] = std::complex<double>(cos(alpha), sin(alpha));
        std::cout << "TESTcomplex";
    }

    compexVector A(n / 2), B(n / 2);
    for (int i = 0; i < n / 2; i++)
    {
        A[i] = analysValue[i * 2];
        B[i] = analysValue[i * 2 + 1];
    }
    compexVector Arecursion = FFT(A);
    compexVector Brecursion = FFT(B);
    compexVector res(n);
    for ( int i = 0; i < n; i++)
    {
        res[i] = Arecursion[i % (n / 2)] + W[i] * Brecursion[i % (n / 2)];
        std::cout << res[i];
    }
    return res;
}

int main()
{
    std::cout << "TESTmain";
    int N = 128;
    int f0 = 100;
    double w = 2 * M_PI * f0;

    compexVector harmonic;

    for (int i = 0; i < N; i++)
    {
        harmonic[i] = sin(w * i);
        std::cout << "TESTharm";
    }

    FFT(harmonic);

    for (int i = 0; i < N; i++)
    {
        std::cout <<" " << harmonic[i] << " ";
    }

    return 0;
}
