#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.14159

float *splitEven(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 2, sizeof(float));
    for (unsigned i = 0; i < len / 2; i++)
        ptr[i] = y[i * 2];

    return ptr;
}

float *splitOdd(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 2, sizeof(float));
    for (unsigned i = 0; i < len / 2; i++)
        ptr[i] = y[i * 2 + 1];

    return ptr;
}

void descrete_fourier_transformation(unsigned no_of_samples, float *y, float *x_real, float *x_imag)
{
    float phi;
    for (unsigned i = 0; i < no_of_samples; i++)
    {

        x_real[i] = 0.0f;
        x_imag[i] = 0.0f;
        for (unsigned j = 0; j < no_of_samples; j++)
        {
            phi = (2 * pi * i * j) / (float)no_of_samples;
            x_real[i] += y[j] * cos(phi);
            x_imag[i] -= y[j] * sin(phi);
        }

        // freq_amp[i] = sqrt(pow(x_real, 2) + pow(x_imaginary, 2)) / (no_of_samples / 2);
        // freq_phase[i] = tan(x_imaginary / x_real) / (no_of_samples / 2) * (180 / pi) ;
    }
}

void redix2FFT(float *x_component_real, float *x_component_img, float *x_even_series_real, float *x_even_series_img, float *x_odd_series_real, float *x_odd_series_img, unsigned i, unsigned length)
{
    float cosphi, sinphi;
    unsigned half_len = length / 2;

    cosphi = cos(-(2 * pi * i) / length);
    sinphi = sin(-(2 * pi * i) / length);

    // printf("cos: %f, sin: %f j,\n", cosphi,sinphi);

    x_component_real[i] = x_even_series_real[i] + cosphi * x_odd_series_real[i] - sinphi * x_odd_series_img[i]; // first half real component
    x_component_img[i] = x_even_series_img[i] + cosphi * x_odd_series_img[i] + sinphi * x_odd_series_real[i];   // first imaginary component

    x_component_real[i + half_len] = x_even_series_real[i] - cosphi * x_odd_series_real[i] + sinphi * x_odd_series_img[i]; // second real component
    x_component_img[i + half_len] = x_even_series_img[i] - cosphi * x_odd_series_img[i] - sinphi * x_odd_series_real[i];   // second imaginary component
}

float **radix5FFT(float *y, unsigned i, unsigned length)
{
    float frac = -(pi / length);
    float cos2pi = cos(2 * frac);
    float cos4pi = cos(4 * frac);
    float sin2pi = sin(2 * frac);
    float sin4pi = sin(4 * frac);

    float **x_components = (float **)malloc(2 * sizeof(float *));
    x_components[0] = (float *)calloc(5, sizeof(float));
    x_components[1] = (float *)calloc(5, sizeof(float));

    x_components[0][0] = y[0] + y[1] + y[2] + y[3] + y[4];

    x_components[0][4] = x_components[0][1] = y[0] + y[1] * cos2pi + y[2] * cos4pi + y[3] * cos4pi + y[4] * cos2pi;
    x_components[1][1] = y[1] * sin2pi + y[2] * sin4pi + y[3] * -sin4pi + y[4] * -sin2pi;

    x_components[0][3] = x_components[0][2] = y[0] + y[1] * cos4pi + y[2] * cos2pi + y[3] * cos2pi + y[4] * cos4pi;
    x_components[1][2] = y[1] * sin4pi + y[2] * -sin2pi + y[3] * sin2pi + y[4] * -sin4pi;

    x_components[1][3] = -x_components[1][2];
    x_components[1][4] = -x_components[1][1];

    return x_components;
}

float **point10fft(float *y, unsigned length)
{
    float *y_even, *y_odd;

    float **x_even_series, **x_odd_series, **x_component;

    y_even = splitEven(y, length);
    y_odd = splitOdd(y, length);

    x_even_series = radix5FFT(y_even, 0, 5);
    x_odd_series = radix5FFT(y_odd, 0, 5);

    x_component = (float **)malloc(2 * sizeof(float));
    x_component[0] = (float *)calloc(length, sizeof(float)); // x_component[0] storing real component
    x_component[1] = (float *)calloc(length, sizeof(float)); // x_component[1] storing complex component

    for (unsigned i = 0; i < length / 2; i++)
    {
        redix2FFT(x_component[0], x_component[1], x_even_series[0], x_even_series[1], x_odd_series[0], x_odd_series[1], i, length);
    }
    return x_component;
}


float **cooleyTukey(float *y, unsigned length)
{
    float **x_even_series, **x_odd_series;
    float **x_component;

    unsigned size;
    if (length == 1)
    {
        x_component = (float **)malloc(2 * sizeof(float *));
        x_component[0] = (float *)calloc(1, sizeof(float));
        x_component[1] = (float *)calloc(1, sizeof(float));
        x_component[0][0] = *y;

        return x_component;
    }
    else
    {
        float *y_even, *y_odd;
        float cosphi, sinphi;
        y_even = splitEven(y, length);
        y_odd = splitOdd(y, length);

        x_even_series = cooleyTukey(y_even, length / 2);
        x_odd_series = cooleyTukey(y_odd, length / 2);

        x_component = (float **)malloc(2 * sizeof(float));
        x_component[0] = (float *)calloc(length, sizeof(float)); // x_component[0] storing real component
        x_component[1] = (float *)calloc(length, sizeof(float)); // x_component[1] storing complex component

        unsigned half_len = length / 2;

        // i is k in math eq. j is 0 - 4 inside radix 5 karnel
        for (unsigned i = 0; i < half_len; i++)
        {
            redix2FFT(x_component[0], x_component[1], x_even_series[0], x_even_series[1], x_odd_series[0], x_odd_series[1], i, length);
        }

        free(x_even_series[0]);
        free(x_even_series[1]);
        free(x_even_series);

        free(x_odd_series[0]);
        free(x_odd_series[1]);
        free(x_odd_series);

        free(y_even);
        free(y_odd);

        return x_component;
    }
}

void main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    float *y, *x_real, *x_imag;
    unsigned n;

    scanf("%d", &n);

    y = (float *)calloc(n, sizeof(float));

    for (unsigned i = 0; i < n; i++)
        scanf("%f", &y[i]);

    // float **radix_5_components = fivePointDescreteFourierTransform(y);

    // for (unsigned i = 0; i < 5; i++)
    //     printf("[%f, %f j, ]\n", radix_5_components[0][i], radix_5_components[1][i]);
    // printf("\n");

    float **fft_components = point10fft(y, 10);

    for (unsigned i = 0; i < 10; i++)
        printf("[%f, %f j, ]\n", fft_components[0][i], fft_components[1][i]);
    printf("\n");

    // printf("imagainary: ");
    // for (unsigned i = 0; i < 8; i++)
    //     printf("", fft_components[1][i]);

    // printf("\n");

    // descrete_fourier_transformation(8, y, x_real, x_imag);

    // printf("X_real: ");
    // for (unsigned i = 0; i < 8; i++)
    //     printf("%f, ", x_real[i]);

    // printf("\n");
    // printf("X_imaginary: ");
    // for (unsigned i = 0; i < 8; i++)
    //     printf("%f j, ", x_imag[i]);
    // printf("\n");
}
