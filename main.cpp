
#include "math.h"
#include "stdio.h"

#include "FFT.h"

#define FR_TABLE_HEIGHT 20
#define TAPS_NUM 256

void determineBaseFrequency(float * Pow, int N) {
    float frequency = 0.f;
    float fm = 0;
    float vm = 0;
    float va = 0;
    float vb = 0;

    for (int i = 1; i <= N; i++) {
        if (vm < Pow[i]) {
            vm = Pow[i];
            fm = i;
            va = Pow[i - 1];
            vb = Pow[i + 1];
        }
    }

    frequency = fm - va / vm + vb / vm;

    printf("Base frequency?: %f\n", frequency);
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
    for (i = 0; i < TAPS_NUM; i++) {
        Re[i] += cos(p * i) * amp;
        Im[i] = 0.0;
    }
}

void generateSquare(float * Re, float * Im, int N, float f, float amp) {
    int i;
    float p = f * M_TWOPI / N;
    for (i = 0; i < TAPS_NUM; i++) {
        Re[i] += cos(p * i) > 0 ? amp : -amp;
        Im[i] = 0.0;
    }
}

void fillPow(float * Re, float * Im, float * Pow, int N) {
    int i;
    float d = 1.f;
    for (i = 0; i < N; i++) {
        Pow[i] = (Re[i] * Re[i] + Im[i] * Im[i]);
        if (d < Pow[i]) d = Pow[i];
    }
    for (i = 0; i < N; i++) {
        Pow[i] /= d;
    }
}

void test() {
    static float Re[TAPS_NUM];
    static float Im[TAPS_NUM];
    static float Pow[TAPS_NUM];
    int i;
    int j;

    generateSine(Re, Im, TAPS_NUM, 10.f, 0.7f);
    generateSquare(Re, Im, TAPS_NUM, 35.45f, 0.8f);
    generateSine(Re, Im, TAPS_NUM, 78.f, 0.5f);

    FFT(Re, Im, TAPS_NUM, log2(TAPS_NUM), FT_DIRECT);

    fillPow(Re, Im, Pow, TAPS_NUM);

    printf("%10s | %10s | %10s\n", "Real", "Imaginary", "Power");
    for (i = 0; i < TAPS_NUM; i++) {
        printf("%10.6f | %10.6f | %10.6f\n", Re[i], Im[i], Pow[i]);
    }

    for (j = FR_TABLE_HEIGHT - 1; j >= 0; j--) {
        for (i = 1; i < TAPS_NUM / 2 + 1; i++) {
            if (Pow[i] * (float) FR_TABLE_HEIGHT >= j) {
                printf("#");
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }

    determineBaseFrequency(Pow, (TAPS_NUM / 2));

    fflush(stdout);
}

int main(int argc, char** argv) {
    test();
}