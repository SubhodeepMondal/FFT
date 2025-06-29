#ifndef COMPOSITE_RADIX_DFT_HPP
#define COMPOSITE_RADIX_DFT_HPP
#include <cmath>
#include <memory>

void radix_2_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j);
void radix_3_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j);
void radix_5_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j);
void radix_7_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j);
void radix_n_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j, unsigned length);
void mixed_radix_cooley_tukey(std::unique_ptr<float[]> &x, float *y,
                              unsigned length, unsigned stride, unsigned depth,
                              std::unique_ptr<unsigned[]> &arr,
                              std::unique_ptr<unsigned[]> &stride_arr,
                              unsigned k);
void fft(std::unique_ptr<std::unique_ptr<float[]>[]> &fft_components, float *y,
         unsigned length);
constexpr double pi{3.14159265358979323846};
#endif // COMPOSITE_RADIX_DFT_HPP
