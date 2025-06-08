#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <immintrin.h>

int main()
{
    unsigned no_of_elements = 1 << 26;
    unsigned shift = 17;
    shift <<= 28;

    srand(time(0));

    std::cout << no_of_elements << std::endl;

    float *a, *b, *d_lin, *d_vec;

    a = (float *)malloc(no_of_elements * sizeof(float));
    b = (float *)malloc(no_of_elements * sizeof(float));
    d_lin = (float *)malloc(no_of_elements * sizeof(float));
    d_vec = (float *)malloc(no_of_elements * sizeof(float));

    for (unsigned i = 0; i < no_of_elements; i++)
    {
        a[i] = (float)rand() / shift;
        b[i] = (float)rand() / shift;
    }


    std::chrono::steady_clock::time_point linear_start_time = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < no_of_elements; i++)
    {
        d_lin[i] = a[i] + b[i];
    }
    std::chrono::steady_clock::time_point linear_finish_time = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point vec_start_time = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < no_of_elements; i += 8)
    {
    
        __m256 a_vec = _mm256_loadu_ps(&a[i]);
        __m256 b_vec = _mm256_loadu_ps(&b[i]);

        __m256 c_vec = _mm256_add_ps(a_vec, b_vec);

        _mm256_storeu_ps(&d_vec[i], c_vec);
    }
    std::chrono::steady_clock::time_point vec_finish_time = std::chrono::steady_clock::now();

    // for (unsigned i = 0; i < no_of_elements; i++)
    // {
    //     if(!(i%8))
    //         std::cout << std::endl;
    //     std::cout << d_vec[i]  << ", ";
    // }
    std::cout << std::endl;
    auto escaped_lin = std::chrono::duration_cast<std::chrono::microseconds>(linear_finish_time-linear_start_time);
    auto escaped_vec = std::chrono::duration_cast<std::chrono::microseconds>(vec_finish_time-vec_start_time);
    std::cout << "Time(microseconds) spend of linear sum: " << escaped_lin.count() << std::endl;
    std::cout << "Time(microseconds) spend of vectorized sum: " << escaped_vec.count() << std::endl;

    return 0;
}