#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.14159

// void radix_5_FFT(float *y_0, float *y_1, float *y_2, float *y_3, float *y_4,
// float *x_0, float *x_1, float *x_2, float *x_3, float *x_4, unsigned i)
// {
//     float frac = -(i * pi / 5);
//     float cos2pi = cos(2 * frac);
//     float cos4pi = cos(4 * frac);
//     float sin2pi = sin(2 * frac);
//     float sin4pi = sin(4 * frac);
//     x_0[0] = y_0[0] + y_1[0] + y_2[0] + y_3[0] + y_4[0];
//     x_4[0] = x_1[0] = (y_0[0] + y_1[0] * cos2pi + y_2[0] * cos4pi + y_3[0] *
//     cos4pi + y_4[0] * cos2pi) - (y_1[1] * sin2pi + y_2[1] * sin4pi + y_3[1] *
//     -sin4pi + y_4[1] * -sin2pi); x_1[1] = (y_1[0] * sin2pi + y_2[0] * sin4pi
//     + y_3[0] * -sin4pi + y_4[0] * -sin2pi) + (y_1[1] * cos2pi + y_2[1] *
//     cos4pi + y_3[1] * cos4pi + y_4[1] * cos2pi); x_3[0] = x_2[0] = y_0[0] +
//     y_1[0] * cos4pi + y_2[0] * cos2pi + y_3[0] * cos2pi + y_4[0] * cos4pi -
//     (y_1[1] * sin4pi + y_2[1] * -sin2pi + y_3[1] * sin2pi + y_4[1] *
//     -sin4pi); x_2[1] = y_1[0] * sin4pi + y_2[0] * -sin2pi + y_3[0] * sin2pi +
//     y_4[0] * -sin4pi + (y_1[1] * cos4pi + y_2[1] * cos2pi + y_3[1] * cos2pi +
//     y_4[1] * cos4pi); x_3[1] = -x_2[1]; x_4[1] = -x_1[1];
// }

void radix_5_FFT(float **y, float **x, unsigned x_stride, unsigned y_stride,
                 unsigned k) {
  float x_temp[2], cosinetheta[5], sinetheta[5], twiddle_fraction;
  unsigned i, j, idx_x, idx_y;

  twiddle_fraction = -(2.0 * pi * k) / 5;
  for (i = 0; i < 5; i++) {
    cosinetheta[i] = cos(twiddle_fraction * i);
    sinetheta[i] = sin(twiddle_fraction * i);
  }

  x_temp[0] = x_temp[1] = 0;

  for (i = 0; i < 5; i++) {
    idx_x = i;// * x_stride;
    idx_y = i;// * y_stride;

    // printf("%f, %f\t", y[idx_y][0], y[idx_y][1]);
    x[idx_x][0] = y[idx_y][0] * cosinetheta[i] - y[idx_y][1] * sinetheta[i];
    x[idx_x][1] = y[idx_y][0] * sinetheta[i] + y[idx_y][1] * cosinetheta[i];
  }
}

void twiddleMultiplication(float **x, unsigned stride, unsigned length,
                           unsigned radix) {
  unsigned i, j, idx;
  float twiddle, costheta, sinetheta;
  float **x_temp;

  x_temp = (float **)malloc(length * sizeof(float *));
  for (i = 0; i < length; i++)
    x_temp[i] = (float *)calloc(2, sizeof(float));

  twiddle = -(2.0 * pi) / length;

  for (i = 0; i < length; i++) {
    for (j = 0; j < radix; j++) {
      costheta = cos(i * j * twiddle);
      sinetheta = sin(i * j * twiddle);

      idx = j * stride;

      x_temp[i][0] += x[idx][0] * costheta - x[idx][1] * sinetheta;
      x_temp[i][1] += x[idx][0] * sinetheta + x[idx][1] * costheta;
    }
  }
  for (i = 0; i < length; i++) {
    idx = i * stride;
    // printf("%u, ", idx);
    x[idx][0] = x_temp[i][0];
    x[idx][1] = x_temp[i][1];
  }
  for (i = 0; i < length; i++)
    free(x_temp[i]);
  free(x_temp);
}

void mixed_radix_cooley_tukey(float *y, float **x, unsigned length,
                              unsigned stride, unsigned k) {
  unsigned i, j, factor;
  float **y_con;

  if (length == 5) {
    for (i = 0; i < 5; i++) {
      x[i * stride][0] = y[i * stride];
      // printf("%f, ", y[i * stride]);
    }

    radix_5_FFT(x, x, stride, 1, 0);

    twiddleMultiplication(x, stride, length, 5);

    // printf("[ ");
    // for (i = 0; i < length; i++) {
    //   printf("%f,\t%fj\n", x[i*stride][0], x[i*stride][1]);
    // }

    // printf("]\n[");
  } else {
    if (!(length % 5))
      factor = 5;

    switch (factor) {
    case 5: {
      for (i = 0; i < 5; i++)
        mixed_radix_cooley_tukey((y + i), &x[i], length / 5, stride * 5, 0);

      for (i = 0; i < 5; i++)
        radix_5_FFT(&x[i], &x[i], stride, stride, i);

      twiddleMultiplication(x, stride, length, 5);
      break;
    }

    default:
      break;
    }
  }
}

void main() {

  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);

  unsigned i, length;
  float *y;

  scanf("%u", &length);

  y = (float *)malloc(length * sizeof(float));

  for (i = 0; i < length; i++)
    y[i] = rand() % 10;

  printf("[ ");
  for (i = 0; i < length; i++) {
    printf("%0.0f, ", y[i]);
  }
  printf("]\n[");

  float **fft_component = (float **)malloc(length * sizeof(float *));
  for (i = 0; i < length; i++)
    fft_component[i] = (float *)calloc(2, sizeof(float));
  mixed_radix_cooley_tukey(y, fft_component, length, 1, 1);

  for (i = 0; i < length; i++) {
    printf("%f,\t%fj\n", fft_component[i][0], fft_component[i][1]);
  }
}
