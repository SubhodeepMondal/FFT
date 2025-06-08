#include <iostream>
extern "C" int foo(void);

int main() {
	std::cout<<"Foo returns "<< foo() <<"\n";
	return 0;
}