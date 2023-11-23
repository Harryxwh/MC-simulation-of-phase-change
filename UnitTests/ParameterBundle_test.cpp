#include "catch2/catch_amalgamated.hpp"
#include "../IsingSystem.hpp"
#include "../ParameterBundle.hpp"
/*
TEST_CASE("XML parameter paser test") {

	CubicLatticeParameterBundle Parameters("../InputParam.xml", true);

	REQUIRE(Parameters._L0() == 2);
	REQUIRE(Parameters._L1() == 2);
	REQUIRE(Parameters._L2() == 2);
	REQUIRE(Parameters._N_TemperaturePoints() == 56);
	REQUIRE_THAT(Parameters._Tmax(), Catch::Matchers::WithinRel(6.5, 1e-10));
	REQUIRE_THAT(Parameters._Tmin(), Catch::Matchers::WithinRel(0.1, 1e-10));

}


TEST_CASE("CubicLatticeIsingSystem Observables Brute-Force Counting Using XML") {
	CubicLatticeParameterBundle Parameters("../InputParam.xml", true);
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = Parameters._L0();
	L[1] = Parameters._L1();
	L[2] = Parameters._L2();

	std::vector<double> beta{Parameters._Tmin(), 0.5*(Parameters._Tmin()+Parameters._Tmax()), Parameters._Tmax()};

	const int n_spins = CubicLatticeIsingSystem::eval_n_spins(L);
	CubicLatticeDataBundle data_bundle(n_spins, beta);
	CubicLatticeIsingSystem test_system(L,data_bundle);

	test_system.exact_count();

	const double Jabs = std::fabs(test_system._J()); // 1
	const double mz = 0.0;
	for ( unsigned int beta_idx = 0; beta_idx < beta.size(); beta_idx ++ ) {

		const double b = beta[beta_idx] * Jabs;

		const double z = 2.0 * std::exp(-24.0 * b) * std::pow(1.0 + std::exp(4.0 * b), 4.0)
						* std::pow(1.0 + std::exp(8.0 * b), 2.0)
		 				* (1.0 - 4.0 * std::exp(4.0 * b) + 8.0 * std::exp(8.0 * b) 
						- 4.0 * std::exp(12.0 * b) + std::exp(16.0 * b));
		
		const double energy = - 24.0 * (-1.0 + 3.0 * std::exp(4.0 * b) - 5.0 * std::exp(8.0 * b) +
									 3.0 * std::exp(12.0 * b) - 3.0 * std::exp(16.0 * b) 
									 + 5.0 * std::exp(20.0 * b) - 3.0 * std::exp(24.0 * b)
									 + std::exp(28.0 * b)) /((1.0 + std::exp(4.0 * b)) * (1.0 + std::exp(8.0 * b))
									 * (1.0 - 4.0 * std::exp(4.0 * b) + 8.0 * std::exp(8.0 * b) 
										- 4.0 * std::exp(12.0 * b) + std::exp(16.0 * b))) ;

		const double mz_sq = 64.0 * std::exp(-12.0 * b)* (1.0 + 3.0 * std::exp(8.0 * b)+ 8.0 * std::exp(12.0 * b) 
									+ 3.0 * exp(16.0 * b) + 6.0 * std::exp(20.0 * b) + 9.0 * std::exp(24.0 * b) 
									+ 2.0 * std::exp(36.0 * b)) / z;
		
		const double ene_sq = 384.0 * std::exp(-24.0 * b) * (3.0 + 6.0 * std::exp(12.0 * b) + 5.0 * std::exp(16.0 * b)
											+ 2.0 * std::exp(20.0 * b) + 2.0 * std::exp(28.0 * b) 
											+ 5.0 * std::exp(32.0 * b)+ 6.0 * std::exp(36.0 * b)
											+ 3.0 * std::exp(48.0 * b)) / z;
		
		REQUIRE_THAT(data_bundle._get_exact_magz (beta_idx), Catch::Matchers::WithinRel(mz, 1e-12));
		REQUIRE_THAT(data_bundle._get_exact_magz_sq(beta_idx), Catch::Matchers::WithinRel(mz_sq, 1e-12));
		REQUIRE_THAT(data_bundle._get_exact_energy(beta_idx), Catch::Matchers::WithinRel(energy, 1e-12));
		REQUIRE_THAT(data_bundle._get_exact_C(beta_idx), Catch::Matchers::WithinAbs(beta[beta_idx] * beta[beta_idx] *(ene_sq - energy * energy), 1e-10));
	}
}
*/