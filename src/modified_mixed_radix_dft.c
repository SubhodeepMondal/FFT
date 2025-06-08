#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14159

float *split5Mod0(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 5, sizeof(float));
    for (unsigned i = 0; i < len / 5; i++)
        ptr[i] = y[i * 5];

    return ptr;
}

float *split5Mod1(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 5, sizeof(float));
    for (unsigned i = 0; i < len / 5; i++)
        ptr[i] = y[i * 5 + 1];

    return ptr;
}

float *split5Mod2(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 5, sizeof(float));
    for (unsigned i = 0; i < len / 5; i++)
        ptr[i] = y[i * 5 + 2];

    return ptr;
}

float *split5Mod3(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 5, sizeof(float));
    for (unsigned i = 0; i < len / 5; i++)
        ptr[i] = y[i * 5 + 3];

    return ptr;
}

float *split5Mod4(float *y, unsigned len)
{
    float *ptr = (float *)calloc(len / 5, sizeof(float));
    for (unsigned i = 0; i < len / 5; i++)
        ptr[i] = y[i * 5 + 4];

    return ptr;
}

// radix 5 dft hard coded
void radix_5_FFT(float *y_0, float *y_1, float *y_2, float *y_3, float *y_4, float *x_0, float *x_1, float *x_2, float *x_3, float *x_4)
{
    float frac = -(pi / 5);
    float cos2pi = cos(2 * frac);
    float cos4pi = cos(4 * frac);
    float sin2pi = sin(2 * frac);
    float sin4pi = sin(4 * frac);

    x_0[0] = y_0[0] + y_1[0] + y_2[0] + y_3[0] + y_4[0];

    x_4[0] = x_1[0] = (y_0[0] + y_1[0] * cos2pi + y_2[0] * cos4pi + y_3[0] * cos4pi + y_4[0] * cos2pi) - (y_1[1] * sin2pi + y_2[1] * sin4pi + y_3[1] * -sin4pi + y_4[1] * -sin2pi);
    x_1[1] = (y_1[0] * sin2pi + y_2[0] * sin4pi + y_3[0] * -sin4pi + y_4[0] * -sin2pi) + (y_1[1] * cos2pi + y_2[1] * cos4pi + y_3[1] * cos4pi + y_4[1] * cos2pi);

    x_3[0] = x_2[0] = y_0[0] + y_1[0] * cos4pi + y_2[0] * cos2pi + y_3[0] * cos2pi + y_4[0] * cos4pi - (y_1[1] * sin4pi + y_2[1] * -sin2pi + y_3[1] * sin2pi + y_4[1] * -sin4pi);
    x_2[1] = y_1[0] * sin4pi + y_2[0] * -sin2pi + y_3[0] * sin2pi + y_4[0] * -sin4pi + (y_1[1] * cos4pi + y_2[1] * cos2pi + y_3[1] * cos2pi + y_4[1] * cos4pi);

    x_3[1] = -x_2[1];
    x_4[1] = -x_1[1];
}

// Twiddle multiplication for each point
void fivePointTwiddleMultiplication(float **r_inp, float *x, unsigned i, unsigned length)
{
    float frac;
    float cosphi[5], sinephi[5];
    unsigned j, k, one_5h_len;

    frac = -((2 * pi) / length);

    for (j = 0; j < 5; j++)
    {
        cosphi[j] = cos(frac * i * j);
        sinephi[j] = sin(frac * i * j);
    }

    // for (j = 0; j < 5; j++)
    //     for (k = 0; k < 5; k++)
    //         printf("cos[%d][%d]: %f, sin[%d][%d]: %f\n",j, k, cosphi[j][k], j, k, sinephi[j][k]);

    // real part calculation
    for (j = 0; j < 5; j++)
        x[0] += r_inp[j][0] * cosphi[j] - r_inp[j][1] * sinephi[j];

    // imaginary part calculation
    for (j = 0; j < 5; j++)
        x[1] += r_inp[j][1] * cosphi[j] + r_inp[j][0] * sinephi[j];
}

float **cooley_tokey(float *y, unsigned length)
{
    float **x_comp, **y_con;
    float *mod_0_ptr, *mod_1_ptr, *mod_2_ptr, *mod_3_ptr, *mod_4_ptr, *mod_5_ptr, *mod_6_ptr;
    float **mod_0, **mod_1, **mod_2, **mod_3, **mod_4, **mod_5, **mod_6;
    unsigned i, factor, mod_two, mod_three, mod_five, mod_seven, iter_len;

    mod_five = length % 5;

    if (!mod_five)
        factor = 5;

    switch (factor)
    {
    case 5:
    {
        // Base case if input length is 5, Radix 5 dft is invoked.
        if (length == 5)
        {
            y_con = (float **)malloc(5 * sizeof(float *));
            x_comp = (float **)malloc(5 * sizeof(float *));

            for (i = 0; i < 5; i++)
            {
                y_con[i] = (float *)calloc(2, sizeof(float));
                x_comp[i] = (float *)calloc(2, sizeof(float));

                y_con[i][0] = y[i];
            }
            radix_5_FFT(y_con[0], y_con[1], y_con[2], y_con[3], y_con[4], x_comp[0], x_comp[1], x_comp[2], x_comp[3], x_comp[4]);

            return x_comp;
        }
        else
        {
            float *r_inp[5];

            iter_len = length / 5;

            // spliting the dataset based on mod 5 indexes.
            mod_0_ptr = split5Mod0(y, length);
            mod_1_ptr = split5Mod1(y, length);
            mod_2_ptr = split5Mod2(y, length);
            mod_3_ptr = split5Mod3(y, length);
            mod_4_ptr = split5Mod4(y, length);

            // Recusively calling cooley_tokey for every mod of 5 i.e (0,...,4)
            mod_0 = cooley_tokey(mod_0_ptr, iter_len);
            mod_1 = cooley_tokey(mod_1_ptr, iter_len);
            mod_2 = cooley_tokey(mod_2_ptr, iter_len);
            mod_3 = cooley_tokey(mod_3_ptr, iter_len);
            mod_4 = cooley_tokey(mod_4_ptr, iter_len);

            // Allocating memory for x_component (output i.e real and complex components)
            x_comp = (float **)malloc(length * sizeof(float *));
            for (i = 0; i < length; i++)
                x_comp[i] = (float *)calloc(2, sizeof(float));

            // Twiddle Multiplication for each point.
            for (i = 0; i < length; i++)
            {

                r_inp[0] = mod_0[i % iter_len];
                r_inp[1] = mod_1[i % iter_len];
                r_inp[2] = mod_2[i % iter_len];
                r_inp[3] = mod_3[i % iter_len];
                r_inp[4] = mod_4[i % iter_len];

                fivePointTwiddleMultiplication(r_inp, x_comp[i], i, length);
            }

            // Free all allocated memory.
            free(mod_0_ptr);
            free(mod_1_ptr);
            free(mod_2_ptr);
            free(mod_3_ptr);
            free(mod_4_ptr);
            
            free(mod_0);
            free(mod_1);
            free(mod_2);
            free(mod_3);
            free(mod_4);

            return x_comp;
        }

        break;
    }
    }
}

void main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    unsigned i, length;
    float *y;

    scanf("%u", &length);

    y = (float *)malloc(length * sizeof(float));

    for (i = 0; i < length; i++)
        y[i] = rand() % 10;

    float **fft_component = cooley_tokey(y, length);

    for (i = 0; i < length; i++)
    {
        printf("%0.0f, ", y[i]);
    }

    printf("\n");

    for (i = 0; i < length; i++)
    {
        printf("%f,\t%fj\n", fft_component[i][0], fft_component[i][1]);
    }
}