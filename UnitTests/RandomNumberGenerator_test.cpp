
#include "catch2/catch_amalgamated.hpp"
#include "../RandomNumberGenerator.hpp"

/*
TEST_CASE("Random Number (64-bit Mersenne Twister in C++ library)") {
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed);


	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.455093, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.137216, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.967235, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.771364, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.244447, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.565477, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.871100, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.441108, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.999349, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.340951, 1e-1));
};

TEST_CASE("Save/Load Random Number (64-bit Mersenne Twister in C++ library)") {
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed);

	
	for (int n = 0; n < 99; n++) mtwist.gen_rand01();
	std::string save_data = mtwist.save();
	const double rn_sample = mtwist.gen_rand01();

	for (int n = 0; n < 999; n++) mtwist.gen_rand01();
	mtwist.load(save_data);
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(rn_sample, 1e-10));

};
*/