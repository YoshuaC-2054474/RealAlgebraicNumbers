#include <iostream>
#define GINAC_DLL
#include <ginac/ginac.h>

using namespace GiNaC;

int main() {
	symbol x("x"), y("y");
	ex poly = x * y + pow(x, 2);
	std::cout << poly << std::endl;
	return 0;
}
