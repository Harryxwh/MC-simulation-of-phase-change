#include "catch2/catch_amalgamated.hpp"
#include "../IsingSystem.hpp"


/*
TEST_CASE("CubicLatticeIsingSystem") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 3;
	L[1] = 3;
	L[2] = 3;
	CubicLatticeIsingSystem test_system(L);
	REQUIRE(test_system._n_spins() == 27);
	REQUIRE(test_system.eval_energy() == -81.0); // energy is a double
	REQUIRE(test_system.eval_mz() == 27.0); 
}


TEST_CASE("CubicLatticeIsingSystem 2") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 3;
	L[1] = 2;
	L[2] = 3;
	CubicLatticeIsingSystem test_system(L);
	test_system.set_spin(2, -1);
	test_system.set_spin(6, -1);
	test_system.set_spin(7, -1);
	test_system.set_spin(9, -1);
	test_system.set_spin(11, -1);
	
	REQUIRE(test_system._n_spins() == 18);
	REQUIRE(test_system.eval_energy() == -10.0);
	REQUIRE(test_system.eval_mz() == 8.0);
}


TEST_CASE("CubicIsingSystem Int. Rep. of States (1)") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 3;
	L[1] = 2;
	L[2] = 3;
	CubicLatticeIsingSystem test_system(L);
	test_system.set_state_by_code(259387);
	REQUIRE(test_system._sz(0) == 1);
	REQUIRE(test_system._sz(1) == 1);
	REQUIRE(test_system._sz(2) == -1);
	REQUIRE(test_system._sz(3) == 1);
	REQUIRE(test_system._sz(4) == 1);
	REQUIRE(test_system._sz(5) == 1);
	REQUIRE(test_system._sz(6) == -1);
	REQUIRE(test_system._sz(7) == -1);
	REQUIRE(test_system._sz(8) == 1);
	REQUIRE(test_system._sz(9) == -1);
	REQUIRE(test_system._sz(10) == 1);
	REQUIRE(test_system._sz(11) == -1);
	REQUIRE(test_system._sz(12) == 1);
	REQUIRE(test_system._sz(13) == 1);
	REQUIRE(test_system._sz(14) == 1);
	REQUIRE(test_system._sz(15) == 1);
	REQUIRE(test_system._sz(16) == 1);
	REQUIRE(test_system._sz(17) == 1);
}

TEST_CASE("CubicLatticeIsingSystem Integer Rep. of States (2)") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 2;
	L[1] = 2;
	L[2] = 2;
	CubicLatticeIsingSystem test_system(L);
	REQUIRE(test_system._maxrep_state() == 255);
	
	test_system.set_state_by_code(0);//all downward
	REQUIRE(test_system.eval_mz() == -8.0);
	REQUIRE(test_system.eval_energy() == -24.0);
	
	test_system.set_state_by_code(1);
	REQUIRE(test_system.eval_mz() == -6.0);
	REQUIRE(test_system.eval_energy() == -12.0);

	test_system.set_state_by_code(2);
	REQUIRE(test_system.eval_mz() == -6.0);
	REQUIRE(test_system.eval_energy() == -12.0);

	test_system.set_state_by_code(3);
	REQUIRE(test_system.eval_mz() == -4.0);
	REQUIRE(test_system.eval_energy() == -8.0);

	test_system.set_state_by_code(18);
	REQUIRE(test_system.eval_mz() == -4.0);
	REQUIRE(test_system.eval_energy() == 0.0);
}

TEST_CASE("CubicLatticeIsingSystem Observables Brute-Force Counting") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 2;
	L[1] = 2;
	L[2] = 2;

	// nan for too large beta
	//const std::vector<double> beta{0, 0.3, 1.1, 10.1, 100.1 };
	//const std::vector<double> beta{0, 0.3, 1.1, 10.1};
	
	const std::vector<double> beta{0, 0.3, 1.1, 5.1};
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
		REQUIRE_THAT(data_bundle._get_exact_C(beta_idx), Catch::Matchers::WithinRel(beta[beta_idx] * beta[beta_idx] *(ene_sq - energy * energy), 1e-8));
	}
}
*/