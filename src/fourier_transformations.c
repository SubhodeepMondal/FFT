#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14159

void waveGenrator(unsigned no_of_freq, unsigned *freq, float *amplitude, unsigned no_of_samples, unsigned sampeling_freq, float *y)
{
    float time_step = 1.0f / sampeling_freq;
    for (unsigned i = 0; i < no_of_freq; i++)
    {
        for (unsigned j = 0; j < no_of_samples; j++)
        {
            y[j] += amplitude[i] * sin(2 * pi * (j * time_step) * freq[i]);
        }
    }
}

void descrete_fourier_transformation(unsigned no_of_samples, float *y, float *freq_amp, float *freq_phase)
{
    float x_real, x_imaginary, phi;

    for (unsigned i = 0; i < no_of_samples / 2; i++)
    {
        x_real = 0.0f;
        x_imaginary = 0.0f;
        for (unsigned j = 0; j < no_of_samples; j++)
        {
            phi = (2 * pi * i * j) / (float)no_of_samples;
            x_real += y[j] * cos(phi);
            x_imaginary += y[j] * sin(phi);
        }

        // freq_amp[i] = sqrt(pow(x_real, 2) + pow(x_imaginary, 2)) / (no_of_samples / 2);
        // freq_phase[i] = tan(x_imaginary / x_real) / (no_of_samples / 2) * (180 / pi) ;
    }
}

void fivePointDescreteFourierTransform(float* y, float* x_real, float* x_imag)
{
    float frac = pi / 5;
    x_real[0] = y[0] + y[1] + y[2] + y[3] + y[4];
    x_real[1] = y[0] + y[1] * cos(2 * frac) + y[2] * cos( 4 * frac) + y[3] * cos(6 * frac) + y[4] * cos(8 * frac);
    x_imag[1] = - (y[1] * sin(2 * frac) + y[2] * sin(4 * frac) + y[3] * sin(6 * frac) + y[4] * sin(8 * frac));
    
    x_real[2] = y[0] + y[1] * cos(2 * 2 * frac) + y[2] * cos( 2 * 4 * frac) + y[3] * cos(2 * 6 * frac) + y[4] * cos(2 * 8 * frac);
    x_imag[2] = - (y[1] * sin(2 * 2 * frac) + y[2] * sin(2 * 4 * frac) + y[3] * sin(2 * 6 * frac) + y[4] * sin(2 * 8 * frac));

    x_real[3] = x_real[2];
    x_imag[3] = x_real[3];
    
    x_real[4] = x_real[1];
    x_real[4] = x_real[1];
}

void main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    unsigned no_of_samples, no_of_freq, sampeling_frequency;
    unsigned *freq;
    float amplitude;
    float *y, *amp, *freq_amp, *freq_phase;

    // printf("Enter sampling frequency:\n");
    scanf("%d", &sampeling_frequency);
    // printf("Enter no of samples:\n");
    scanf("%d", &no_of_samples);

    // printf("Enter no of frequencies: \n");
    scanf("%d", &no_of_freq);

    freq = (unsigned *)malloc(no_of_freq * sizeof(unsigned));
    amp = (float *)malloc(no_of_freq * sizeof(float));

    // printf("Enter %d frequencies:\n", no_of_freq);
    for (unsigned i = 0; i < no_of_freq; i++)
        scanf("%d", &freq[i]);

    printf("Enter %d amplitudes:\n", no_of_freq);
    for (unsigned i = 0; i < no_of_freq; i++)
        scanf("%f", &amp[i]);

    printf("sampling frequency: %d \n No of samples:%d \n No of frequencies:%d\n", sampeling_frequency, no_of_samples, no_of_freq);
    for (unsigned i = 0; i < no_of_freq; i++)
        printf("freq[%d]: %d, amplitude[%d]:%f\n", i + 1, freq[i], i + 1, amp[i]);

    y = (float *)calloc(no_of_samples, sizeof(float));
    freq_amp = (float *)malloc(no_of_samples / 2 * sizeof(float));
    freq_phase = (float *)malloc(no_of_samples / 2 * sizeof(float));

    // printf("no_of samples: %d, amplitude: %f", no_of_samples, amplitude);

    waveGenrator(no_of_freq, freq, amp, no_of_samples, sampeling_frequency, y);

    descrete_fourier_transformation(no_of_samples, y, freq_amp, freq_phase);

    for (unsigned i = 0; i < no_of_samples / 2; i++)
    {
        printf("frequency: %d, frequency amp: %f, frequency phase(degree): %f\n", i, freq_amp[i], freq_phase[i]);
    }
    printf("\n");
}

// Enter sampling frequency:
// 1000
// Enter no of samples:
// 1000
// Enter no of frequencies
// 1
// Enter 1 frequencies80
// Enter 1 amplitudes1.8
// no_of samples: 1000, amplitude: 0.000000