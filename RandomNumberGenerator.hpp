#ifndef RandomNumberGenerator_hpp
#define RandomNumberGenerator_hpp

#include <random>
#include <iostream>
#include <sstream>
#include <string>

class RandomNumberGenerator {
private:
	const int seed_init;
	std::mt19937_64 RndEngine;
	std::uniform_real_distribution<double> DistUnif01;
	std::uniform_int_distribution<> DistUnifSite;
	
public:
	RandomNumberGenerator(const int seed_spec, const int n_spins_spec = 0)
	: seed_init(seed_spec), RndEngine(seed_init), DistUnif01(0.0, 1.0), DistUnifSite(0, n_spins_spec - 1) {};
	~RandomNumberGenerator() {};
	
	/* uniform [0, 1) */
	double gen_rand01() { return DistUnif01(RndEngine); };
	
	int _seed_init(){return seed_init;}
	/* uniform [0, n_spins) */
	int gen_rand_site() { return DistUnifSite(RndEngine) ;}
	
	std::string save() const {
		std::stringstream save_data;
		save_data << RndEngine;
		return save_data.str();
	};

	void load(std::string save_data) {
		std::stringstream restore(save_data);
		restore >> RndEngine;
	};
};
#endif /* RandomNumberGenerator_hpp */
