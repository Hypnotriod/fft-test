
#include <complex>

#include "math.h"
#include "stdio.h"

#include "FFT.h"

#ifndef M_TWOPI
#define M_TWOPI (M_PI * 2)
#endif

#define SIGNAL_TABLE_HEIGHT 20
#define POW_TABLE_HEIGHT 20
#define TAPS_NUM 256

void printBaseFrequency(float * Pow, int N) {
    float frequency = 0;
    float fm = 0;
    float vm = 0;
    float va = 0;
    float vb = 0;
    float vc = 0;
    float vd = 0;
    float ncoef;
    int i;

    for (i = 1; i <= N; i++) {
        if (vm < Pow[i]) {
            vm = Pow[i];
            fm = i;
        }
    }
    
    i = (int)fm;
    if (i == 0) return;
    va = Pow[i - 1];
    vb = i > 1 ? Pow[i - 2] : 0.f;
    vc = Pow[i + 1];
    vd = Pow[i + 2];

    ncoef = vm / (vm + va + vb + vc + vd);
    frequency = fm
            - va / vm * ncoef
            - vb / vm * ncoef * 2.f
            + vc / vm * ncoef
            + vd / vm * ncoef * 2.f;

    printf("Base frequency?: %f\n", frequency);
}

void printFrequencies(float * Pow, int N, float threshold) {
    float frequency = 0;
    float fm = 0;
    float vm = 0;
    float va = 0;
    float vb = 0;
    float vc = 0;
    float vd = 0;
    float ncoef;
    int i;
    int n = 1;

    for (i = 1; i <= N; i++) {
        if (vm < Pow[i]) {
            vm = Pow[i];
            fm = i;
        }
        else if (vm > threshold) {
            i--;
            va = Pow[i - 1];
            vb = i > 1 ? Pow[i - 2] : 0.f;
            vc = Pow[i + 1];
            vd = Pow[i + 2];
            i += 3;
            
            ncoef = vm / (vm + va + vb + vc + vd);
            frequency = fm
                    - va / vm * ncoef
                    - vb / vm * ncoef * 2.f
                    + vc / vm * ncoef
                    + vd / vm * ncoef * 2.f;
            
            printf("Frequency #%i: %f\n", n++, frequency);
            vm = 0.f;
        }
    }
}

void printPeaks(float * Pow, int N, float threshold) {
    float vm = 0;
    int i;
    int n = 1;

    for (i = 1; i <= N; i++) {
        if (vm < Pow[i]) {
            vm = Pow[i];
        }
        else if (vm > threshold) {
            printf("Frequency #%i peak: %f\n", n++, vm);
            vm = 0.f;
        }
    }
}

void hannWindow(float * data, int N) {
    for (int i = 0; i < N; i++) {
        float multiplier = 0.5f * (1 - cos(2 * M_PI * i / (N - 1)));
        data[i] = multiplier * data[i];
    }
}

void generateSine(float * Re, float * Im, int N, float f, float amp) {
    int i;
    float p = f * M_TWOPI / N;
    for (i = 0; i < N; i++) {
        Re[i] += cos(p * i) * amp;
        Im[i] = 0.0;
    }
}

void generateSquare(float * Re, float * Im, int N, float f, float amp) {
    int i;
    float p = f * M_TWOPI / N;
    for (i = 0; i < N; i++) {
        Re[i] += cos(p * i) > 0 ? amp : -amp;
        Im[i] = 0.0;
    }
}

void fillPow(float * Re, float * Im, float * Pow, int N) {
    int i;
    for (i = 0; i < N; i++) {
        Pow[i] = (Re[i] * Re[i] + Im[i] * Im[i]);
    }
}

// https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
inline float sqrt7(float x)
 {
   unsigned int i = *(unsigned int*) &x;
   // adjust bias
   i  += 127 << 23;
   // approximation of square root
   i >>= 1;
   return *(float*) &i;
 }

void normalizePow(float * Pow, int N) {
    int i;
    float p = (4.f / N);
    for (i = 0; i < N; i++) {
        Pow[i] = sqrt7(Pow[i]) * p;
    }
}

void printPow(float * Pow, int N) {
    int i;
    int j;
    for (j = POW_TABLE_HEIGHT - 1; j >= 0; j--) {
        for (i = 1; i < N; i++) {
            if (Pow[i] * (float) POW_TABLE_HEIGHT >= j) {
                printf("#");
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }
}

void printSignal(float * Re, int N) {
    int i;
    int j;
    for (j = SIGNAL_TABLE_HEIGHT - 1; j >= 0; j--) {
        for (i = 0; i < N; i++) {
            if (j > SIGNAL_TABLE_HEIGHT / 2) {
                if (Re[i] * (float) (SIGNAL_TABLE_HEIGHT / 2) >= j - (SIGNAL_TABLE_HEIGHT / 2)) {
                    printf("#");
                } else {
                    printf(" ");
                }
            } else {
                if (Re[i] * (float) (SIGNAL_TABLE_HEIGHT / 2) <= j - (SIGNAL_TABLE_HEIGHT / 2)) {
                    printf("#");
                } else {
                    printf(" ");
                }
            }
        }
        printf("\n");
    }
}

void printCoefficients(float * Re, float * Im, float * Pow, int N) {
    printf("%10s | %10s | %10s\n", "Real", "Imaginary", "Power");
    int i;
    for (i = 0; i < N; i++) {
        printf("%10.6f | %10.6f | %10.6f\n", Re[i], Im[i], Pow[i]);
    }
}

void printDivider(int N) {
    int i;
    for (i = 0; i < N; i++) {
        printf("_");
    }
    printf("\n");
}

void test() {
    static float Re[TAPS_NUM];
    static float Im[TAPS_NUM];
    static float Pow[TAPS_NUM];

    generateSine(Re, Im, TAPS_NUM, 2.5f, 1.f);
    generateSine(Re, Im, TAPS_NUM, 35.21f, 0.8f);
    generateSine(Re, Im, TAPS_NUM, 78.f, 0.5f);
    generateSine(Re, Im, TAPS_NUM, 121.8f, 0.7f);

    printSignal(Re, TAPS_NUM / 2);
    printDivider(TAPS_NUM / 2);

    hannWindow(Re, TAPS_NUM);

    FFT(Re, Im, TAPS_NUM, log2(TAPS_NUM), FT_DIRECT);

    fillPow(Re, Im, Pow, TAPS_NUM);

//    printCoefficients(Re, Im, Pow, TAPS_NUM);

    // Getting more accurate frequencies detection before peaks normalization
    printBaseFrequency(Pow, (TAPS_NUM / 2));
    printFrequencies(Pow, (TAPS_NUM / 2), 1.f / (2.f / TAPS_NUM));
    
    printDivider(TAPS_NUM / 2);
    
    normalizePow(Pow, TAPS_NUM);
    
    printPow(Pow, TAPS_NUM / 2 + 1);
    
    printDivider(TAPS_NUM / 2);
    
    printPeaks(Pow, (TAPS_NUM / 2), 0.2f);

    fflush(stdout);
}

int main(int argc, char** argv) {
    test();
}