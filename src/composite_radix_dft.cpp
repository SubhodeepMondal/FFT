#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <composite_radix_dft.hpp>
// This file implements the composite radix FFT algorithm using the Cooley-Tukey method.
// It supports radix-2, radix-3, radix-5, and radix-7 FFTs, as well as a general
// radix-n FFT for lengths that are not divisible by these radices.
// The algorithm recursively breaks down the FFT computation into smaller FFTs
// based on the input length and the available radices.
// The implementation uses a mixed-radix approach to handle various lengths efficiently.
// The FFT is computed in-place, and the results are stored in a 2D array where each
// row corresponds to a complex number (real and imaginary parts).

#define pi 3.1415926535

void radix_2_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j) {
  float x_temp[2], twiddle_fraction;
  unsigned i, idx;
  printf("Inside radix 2 kernel\n");
  twiddle_fraction = -(2.0 * pi * k) / 2.0;
  x_temp[0] = x_temp[1] = 0;
  for (i = 0; i < 2; i++) {
    idx = i * stride + offset;
    x_temp[0] += y[idx] * cos(twiddle_fraction * i); // real component
    x_temp[1] += y[idx] * sin(twiddle_fraction * i); // complex componnet
  }
  x[j * 2] = x_temp[0];
  x[j * 2 + 1] = x_temp[1];
}

void radix_3_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j) {
  float x_temp[2], twiddle_fraction;
  unsigned i, idx;
  printf("Inside radix 3 kernel\n");
  twiddle_fraction = -(2.0 * pi * k) / 3.0;
  x_temp[0] = x_temp[1] = 0;
  for (i = 0; i < 3; i++) {
    idx = i * stride + offset;
    x_temp[0] += y[idx] * cos(twiddle_fraction * i); // real component
    x_temp[1] += y[idx] * sin(twiddle_fraction * i); // complex componnet
  }
  x[j * 2] = x_temp[0];
  x[j * 2 + 1] = x_temp[1];
}

void radix_5_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j) {
  float x_temp[2], twiddle_fraction;
  unsigned i, idx;
  printf("Inside radix 5 kernel\n");
  twiddle_fraction = -(2.0 * pi * k) / 5.0;
  x_temp[0] = x_temp[1] = 0;
  for (i = 0; i < 5; i++) {
    idx = i * stride + offset;
    x_temp[0] += y[idx] * cos(twiddle_fraction * i); // real component
    x_temp[1] += y[idx] * sin(twiddle_fraction * i); // complex componnet
  }
  x[j * 2] = x_temp[0];
  x[j * 2 + 1] = x_temp[1];
}

void radix_7_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j) {
  float x_temp[2], twiddle_fraction;
  unsigned i, idx;
  printf("Inside radix 7 kernel\n");
  twiddle_fraction = -(2.0 * pi * k) / 5.0;
  x_temp[0] = x_temp[1] = 0;
  for (i = 0; i < 5; i++) {
    idx = i * stride + offset;
    x_temp[0] += y[idx] * cos(twiddle_fraction * i); // real component
    x_temp[1] += y[idx] * sin(twiddle_fraction * i); // complex componnet
  }
  x[j * 2] = x_temp[0];
  x[j * 2 + 1] = x_temp[1];
}

void radix_n_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j, unsigned length) {
  float x_temp[2], twiddle_fraction;
  unsigned i, idx;
  // printf("inside n kernel\n");
  twiddle_fraction = -(2.0 * pi * k) / (float)length;
  x_temp[0] = x_temp[1] = 0;
  for (i = 0; i < length; i++) {
    idx = i * stride + offset;
    x_temp[0] += y[idx] * cos(twiddle_fraction * i); // real component
    x_temp[1] += y[idx] * sin(twiddle_fraction * i); // complex componnet
  }
  x[j * 2] = x_temp[0];
  x[j * 2 + 1] = x_temp[1];
}

void mixed_radix_cooley_tukey(float *x, float *y, unsigned length,
                              unsigned stride, unsigned depth, unsigned *arr,
                              unsigned *stride_arr, unsigned k) {
  unsigned offset = 0;
  unsigned i, idx;
  switch (length) {
  case 7:
    for (i = 0; i < depth; i++)
      offset += arr[i] * stride_arr[i];
    idx = (depth != 0) ? arr[depth - 1] : 0;
    radix_7_FFT(x, y, stride, offset, k, idx);
    return;
  case 5:
    for (i = 0; i < depth; i++)
      offset += arr[i] * stride_arr[i];
    idx = (depth != 0) ? arr[depth - 1] : 0;
    radix_5_FFT(x, y, stride, offset, k, idx);
    return;
  case 3:
    for (i = 0; i < depth; i++)
      offset += arr[i] * stride_arr[i];
    idx = (depth != 0) ? arr[depth - 1] : 0;
    radix_3_FFT(x, y, stride, offset, k, idx);
    return;
  case 2:
    for (i = 0; i < depth; i++)
      offset += arr[i] * stride_arr[i];
    idx = (depth != 0) ? arr[depth - 1] : 0;
    radix_2_FFT(x, y, stride, offset, k, idx);
    return;
  default:
    unsigned radix = 0;
    if (length % 2 == 0) {
      radix = 2;
    } else if (length % 3 == 0) {
      radix = 3;
    } else if (length % 5 == 0) {
      radix = 5;
    } else if (length % 7 == 0) {
      radix = 7;
    }
    if (radix) {
      float x_comp[7][2] = {0};
      for (i = 0; i < radix; i++) {
        arr[depth] = i;
        stride_arr[depth] = stride;
        mixed_radix_cooley_tukey((float *)x_comp, y, length / radix,
                                 stride * radix, depth + 1, arr, stride_arr, k);

        float twiddle_factor = -(2.0 * pi * i * k) / (float)length;
        if (depth == 0) {
          x[0] += x_comp[i][0] * cos(twiddle_factor) -
                  x_comp[i][1] * sin(twiddle_factor); // real component
          x[1] += x_comp[i][1] * cos(twiddle_factor) +
                  x_comp[i][0] * sin(twiddle_factor); // complex component
        } else {
          x[arr[depth - 1] * 2] +=
              x_comp[i][0] * cos(twiddle_factor) -
              x_comp[i][1] * sin(twiddle_factor); // real component
          x[arr[depth - 1] * 2 + 1] +=
              x_comp[i][1] * cos(twiddle_factor) +
              x_comp[i][0] * sin(twiddle_factor); // complex component
        }
      }
    } else {
      for (i = 0; i < depth; i++)
        offset += arr[i] * stride_arr[i];
      idx = (depth != 0) ? arr[depth - 1] : 0;
      radix_n_FFT(x, y, stride, offset, k, idx, length);
      return;
    }
    break;
  }
}

float **fft(float *y, unsigned length) {

  float **fft_components = (float **)calloc(length, sizeof(float *));
  for (int i = 0; i < length; i++)
    fft_components[i] = (float *)calloc(2, sizeof(float));

  unsigned *arr = (unsigned *)calloc(length, sizeof(unsigned));
  unsigned *stride = (unsigned *)calloc(length, sizeof(unsigned));

  for (int i = 0; i < length; i++) {
    mixed_radix_cooley_tukey(fft_components[i], y, length, 1, 0, arr, stride,
                             i);
  }
  free(arr);
  free(stride);
  return fft_components;
}

// void main() {

//   freopen("input.txt", "r", stdin);
//   freopen("output.txt", "w", stdout);

//   unsigned i, length;
//   float *y;

//   scanf("%u", &length);

//   y = (float *)malloc(length * sizeof(float));

//   for (i = 0; i < length; i++)
//     y[i] = rand() % 10;
//   // y[i] = i + 1;
//   float **fft_component = fft(y, length);
//   printf("length: %d\n[ ", length);
//   for (i = 0; i < length; i++) {
//     printf("%0.0f, ", y[i]);
//   }
//   printf("]\n");

//   for (i = 0; i < length; i++) {
//     printf("%f,\t%fj\n", fft_component[i][0], fft_component[i][1]);
//   }
//   for (i = 0; i < length; i++)
//     free(fft_component[i]);
//   free(fft_component);
//   free(y);
// }
