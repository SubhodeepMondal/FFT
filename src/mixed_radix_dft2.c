#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.1415926535

void radix_5_FFT(float *x, float *y, unsigned stride, unsigned offset,
                 unsigned k, unsigned j) {
  float x_temp[2], twiddle_fraction;
  unsigned i, idx;

  twiddle_fraction = -(2.0 * pi * k) / 5.0;
  x_temp[0] = x_temp[1] = 0;
  for (i = 0; i < 5; i++) {

    idx = i * stride + offset;
    x_temp[0] += y[idx] * cos(twiddle_fraction * i); // real component
    x_temp[1] += y[idx] * sin(twiddle_fraction * i); // complex componnet
    // printf("%d, ", idx);
  }
  // printf("\n");
  x[j * 2] = x_temp[0];
  x[j * 2 + 1] = x_temp[1];
}

void mixed_radix_cooley_tukey(float *x, float *y, unsigned length,
                              unsigned stride, unsigned depth, unsigned *arr,
                              unsigned *stride_arr, unsigned k) {
  if (length == 5) {
    unsigned offset = 0;
    for (int i = 0; i < depth; i++)
      offset += arr[i] * stride_arr[i];
    int idx = (depth != 0) ? arr[depth - 1] : 0;
    radix_5_FFT(x, y, stride, offset, k, idx);
  } else {
    float x_comp[5][2] = {0};
    for (int i = 0; i < 5; i++) {
      // printf("%d: ", i);
      arr[depth] = i;
      stride_arr[depth] = stride;
      mixed_radix_cooley_tukey((float *)x_comp, y, length / 5, stride * 5,
                               depth + 1, arr, stride_arr, k);

      float twiddle_factor = -(2.0 * pi * i * k) / (float)length;
      if (depth == 0) {
        x[0] += x_comp[i][0] * cos(twiddle_factor) -
                x_comp[i][1] * sin(twiddle_factor); // real component
        x[1] += x_comp[i][1] * cos(twiddle_factor) +
                x_comp[i][0] * sin(twiddle_factor); // complex component
      } else {
        x[arr[depth -1] * 2] += x_comp[i][0] * cos(twiddle_factor) -
                    x_comp[i][1] * sin(twiddle_factor); // real component
        x[arr[depth -1] * 2 + 1] += x_comp[i][1] * cos(twiddle_factor) +
                        x_comp[i][0] * sin(twiddle_factor); // complex component
      }
    }
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
  return fft_components;
}

void main() {

  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);

  unsigned i, length;
  float *y;
  unsigned *arr, *stride;

  scanf("%u", &length);

  y = (float *)malloc(length * sizeof(float));

  for (i = 0; i < length; i++)
    y[i] = rand() % 10;
  // y[i] = i + 1;

  arr = (unsigned *)malloc(length * sizeof(unsigned));
  stride = (unsigned *)malloc(length * sizeof(unsigned));

  float **fft_component = fft(y, length);
  printf("length: %d\n[ ", length);
  for (i = 0; i < length; i++) {
    printf("%0.0f, ", y[i]);
  }
  printf("]\n");

  for (i = 0; i < length; i++) {
    printf("%f,\t%fj\n", fft_component[i][0], fft_component[i][1]);
  }
}
