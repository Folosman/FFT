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
    if (n == 1) return compexVector(1, analysValue[0]);

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
    compexVector Arecursion = FFT(A);
    compexVector Brecursion = FFT(B);
    compexVector res(n);
    for ( int i = 0; i < n; i++)
    {
        res[i] = Arecursion[i % (n / 2)] + W[i] * Brecursion[i % (n / 2)];
        std::cout << res[i] << " ";
    }
    std::cout << "\n\n\n";
    return res;
}

int main()
{
    int N = 128;
    int f0 = 100;
    double w = f0;

    compexVector harmonic;
    compexVector output;
    double tempVar;

    for (int i = 0; i < N; i++)
    {
        tempVar = sin(w * i);
        harmonic.push_back(tempVar);
    }
    std::cout << "\n\n";
    FFT(harmonic);
    return 0;
}
