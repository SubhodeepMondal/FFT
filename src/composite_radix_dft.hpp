#ifndef COMPOSITE_RADIX_DFT_HPP
#define COMPOSITE_RADIX_DFT_HPP
#include <cmath>

void radix_2_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j);
void radix_3_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j);
void radix_5_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j);
void radix_7_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j);
void mixed_radix_cooley_tukey(float *x, float *y, unsigned length,
                              unsigned stride, unsigned depth, unsigned *arr,
                              unsigned *stride_arr, unsigned k);
void radix_n_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j, unsigned length);
float **fft(float *y, unsigned length);
const double pi{3.14159265358979323846};
#endif // COMPOSITE_RADIX_DFT_HPP
