#include <composite_radix_dft.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

int main() {

  unsigned i, length;
  float *y;
  std::fstream file;
  file.open("../../tst/testfiles/input.txt", std::ios::in);
  //   std::fstream file2("../../tst/testfiles/output.txt", std::ios::out);

  if (file.is_open()) {
    // File operations here
    file >> length;
    std::cout << "Input length: " << length << std::endl;
    if (length <= 0) {
      std::cerr << "Invalid length: " << length << std::endl;
      return 1;
    }
    std::unique_ptr<std::unique_ptr<float []>[]> fft_component;

    y = new float[length];
    fft_component = std::make_unique< std::unique_ptr<float []>[]>(length);

    for (i = 0; i < length; i++) {
      fft_component[i] = std::make_unique<float[]>(2);
      y[i] = rand() % 10;
    }

    std::cout << "Input length: " << length << std::endl;
    std::cout << "Input values: [ ";
    for (i = 0; i < length; i++) {
      std::cout << y[i] << (i < length - 1 ? ", " : " ");
    }
    std::cout << "]" << std::endl;
    fft(fft_component, y, length);
    printf("length: %d\n[ ", length);
    for (i = 0; i < length; i++) {
      printf("%0.0f, ", y[i]);
    }
    printf("]\n");

    for (i = 0; i < length; i++) {
      printf("%f,\t%fj\n", fft_component[i][0], fft_component[i][1]);
    }
    // for (i = 0; i < length; i++)
    //   free(fft_component[i]);
    // free(fft_component);
    // free(y);
  } else {
    std::cerr << "Error opening file." << std::endl;
  }

  return 0;
}
