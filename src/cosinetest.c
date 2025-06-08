#include <stdio.h>
#include <math.h>

#define pi 3.14159

void main()
{
    float cosine_values[5], sine_values[5];
    unsigned i, j, length, radix;
    length = 25;
    radix = 5;

    freopen("output.txt", "w", stdout);

    for (i = 0; i < radix; i++)
    {
        for (j = 0; j < radix; j++)
        {
            cosine_values[j] = cos(-(2 * pi * i * j) / (float)radix);
            sine_values[j] = sin(-(2 * pi * i * j) / (float)radix);
        }

        printf("\n[");
        for (j = 0; j < radix; j++)
        {
            printf("cos[%d]: %f, ", i*j, cosine_values[j]);
        }
        printf("]\n[");

        for (j = 0; j < radix; j++)
        {
            printf("sin[%d]: %f, ", i*j, sine_values[j]);
        }
        printf("]\n");
    }
}