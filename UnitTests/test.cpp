#include <random>
#include <iostream>
#include <vector>
#include "../Lattice.hpp"
#include "../RandomNumberGenerator.hpp"
int main()
{
    std::vector<int> a={3};
	int b;
	std::vector<int> c;
    Site::eval_position_at(123, b, c);
	std::cout<< a[0]<<std::endl;
	
}