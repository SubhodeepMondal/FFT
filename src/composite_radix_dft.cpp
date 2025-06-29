// This file implements the composite radix FFT algorithm using the Cooley-Tukey
// method. It supports radix-2, radix-3, radix-5, and radix-7 FFTs, as well as a
// general radix-n FFT for lengths that are not divisible by these radices. The
// algorithm recursively breaks down the FFT computation into smaller FFTs based
// on the input length and the available radices. The implementation uses a
// mixed-radix approach to handle various lengths efficiently. The FFT is
// computed in-place, and the results are stored in a 2D array where each row
// corresponds to a complex number (real and imaginary parts).
// #include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <composite_radix_dft.hpp>

void radix_2_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j) {
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

void radix_3_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j) {
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

void radix_5_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j) {
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

void radix_7_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j) {
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

void radix_n_FFT(std::unique_ptr<float[]> &x, float *y, unsigned stride,
                 unsigned offset, unsigned k, unsigned j, unsigned length) {
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

void mixed_radix_cooley_tukey(std::unique_ptr<float[]> &x, float *y,
                              unsigned length, unsigned stride, unsigned depth,
                              std::unique_ptr<unsigned[]> &arr,
                              std::unique_ptr<unsigned[]> &stride_arr,
                              unsigned k) {
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
      std::unique_ptr<std::unique_ptr<float[]>[]> x_comp =
          std::make_unique<std::unique_ptr<float[]>[]>(7);
      for (i = 0; i < 7; i++) {
        x_comp[i] = std::make_unique<float[]>(2);
      }
      for (i = 0; i < radix; i++) {
        arr[depth] = i;
        stride_arr[depth] = stride;
        mixed_radix_cooley_tukey(x_comp[0], y, length / radix, stride * radix,
                                 depth + 1, arr, stride_arr, k);

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

void fft(std::unique_ptr<std::unique_ptr<float[]>[]> &fft_components, float *y,
         unsigned length) {

  std::unique_ptr<unsigned[]> arr = std::make_unique<unsigned[]>(length);
  std::unique_ptr<unsigned[]> stride = std::make_unique<unsigned[]>(length);

  for (int i = 0; i < length; i++) {
    mixed_radix_cooley_tukey(fft_components[i], y, length, 1, 0, arr, stride,
                             i);
  }
  // free(arr);
  // free(stride);
}
