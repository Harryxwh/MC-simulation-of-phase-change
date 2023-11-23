#include "catch2/catch_amalgamated.hpp"
#include "../IsingSystem.hpp"
#include "../ParameterBundle.hpp"

/*
TEST_CASE("CubicLatticeIsingSystem_InputParameterTest") {
	CubicLatticeParameterBundle Parameters("../InputParam.xml");
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = Parameters._L0();
	L[1] = Parameters._L1();
	L[2] = Parameters._L2();
	int MaxSweepStep = Parameters._MaxSweepStep();
	int initstep = Parameters._InitStep();
	int seed = Parameters._Seed();
	REQUIRE(L[0] == 3);
	REQUIRE(L[1] == 3);
	REQUIRE(L[2] == 3);
	REQUIRE(MaxSweepStep == 100);
	REQUIRE(initstep == 30);
}
*/
/*
TEST_CASE("FlipAttempt_UsingMetropolis_Test") {
	std::vector<int> L(3);
	L[0] = 2;
	L[1] = 2;
	L[2] = 2;
	CubicLatticeIsingSystem system(L);
	
	system.set_state_by_code(1); // E = -12
	double R = 0.99;
	//Always flip when energy decreases
	REQUIRE(system.FlipAttemptUsingMetropolis(0,0.25,R)==1);

	system.set_state_by_code(0); // E = -24
	//The exact number should be 0.04978
	R = 0.049;
	REQUIRE(system.FlipAttemptUsingMetropolis(0,0.25,R)==1);
	R = 0.050;
	REQUIRE(system.FlipAttemptUsingMetropolis(0,0.25,R)==0);
}


TEST_CASE("Initialization_Test") {
	std::vector<int> L(3);
	L[0] = 2;
	L[1] = 2;
	L[2] = 2;
	CubicLatticeIsingSystem system(L);
	const int seed = 214523645;
	system.CubicHotInitialize();
	REQUIRE(system._sz(0) == -1);
	REQUIRE(system._sz(1) == 1);
	REQUIRE(system._sz(2) == -1);
	REQUIRE(system._sz(3) == -1);
	REQUIRE(system._sz(4) == -1);
	REQUIRE(system._sz(5) == 1);
	REQUIRE(system._sz(6) == -1);
	REQUIRE(system._sz(7) == 1);
}

TEST_CASE("Sweep_Test") {
	std::vector<int> L(3);
	L[0] = 2;
	L[1] = 2;
	L[2] = 2;
	CubicLatticeIsingSystem system(L);
	const int seed = 214523645;
	REQUIRE(system._n_spins()==8);
	RandomNumberGenerator mtwist(seed,system._n_spins());
	int site_idx;
	site_idx = mtwist.gen_rand_site();
	REQUIRE(site_idx ==2);
	site_idx = mtwist.gen_rand_site();
	REQUIRE(site_idx ==1);
	site_idx = mtwist.gen_rand_site();
	REQUIRE(site_idx ==6);
	site_idx = mtwist.gen_rand_site();
	REQUIRE(site_idx ==4);
	site_idx = mtwist.gen_rand_site();
	REQUIRE(site_idx ==3);
	
}
*/
/*
TEST_CASE("CubicLatticeIsingSystem_MCTest") {
	double ene_ave = 0.0;
	double mz_ave = 0.0;
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 3;
	L[1] = 3;
	L[2] = 3;
	const int seed = 214523645;
	double beta =  0.25;
	int MaxSweepStep = 100;
	int initstep =30;
	CubicLatticeIsingSystem system(L);
	system.CubicHotInitialize();
	
	for(int i=0;i<MaxSweepStep;i++)
	{
		system.sweep(beta);
		double ene = system.eval_energy();
		double mz = system.eval_mz();
		if(i<initstep) continue;
		ene_ave += ene/system._n_spins();
	}
	ene_ave /= (MaxSweepStep-initstep);
	REQUIRE_THAT(ene_ave, Catch::Matchers::WithinRel(-1.89397, 1e-5));
}
*/
bool stderr_validation(const int L_spec, const double b_spec, const int n_trials, const double dev_tol = 3.0) {

	//Size
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = L_spec;
	L[1] = L_spec;
	L[2] = L_spec;
	const int n_spins = CubicLatticeIsingSystem::eval_n_spins(L);

	//Seed
	int seed = 214523645;

	//MC parameters
	const double CI = 0.95;
	const int n_bins = 20;
	const double t_inv_0975_19 = 2.093;
	const double stddev_binomial = sqrt(CI * (1.0 - CI) / n_trials);

	const int n_samples_per_bin = 4000;
	const int mcs_thermalization = 1000;
	const int mcs_interval_btwn_bins = 100;

	std::vector<double> beta(1,b_spec);
	CubicLatticeDataBundle data_bundle(n_spins, beta);
	CubicLatticeIsingSystem system(L, data_bundle, seed);

	int n_in_CI_energy = 0;
	int n_in_CI_C = 0;
	int n_in_CI_magz_sq = 0;
	double exact_energy = 0;
	double exact_C = 0;
	double exact_magz_sq = 0;
	const int beta_idx = 0;
	
	for (int i_trial = 0; i_trial < n_trials; i_trial++) {
		
		if (i_trial==0) {
			system.exact_count();
			exact_energy = data_bundle._get_exact_energy(beta_idx)/n_spins;
			exact_C = data_bundle._get_exact_C(beta_idx)/n_spins;
			exact_magz_sq = data_bundle._get_exact_magz_sq(beta_idx)/n_spins/n_spins;
		}
		

		system.CubicHotInitialize();
		int i = 0;

		double ave_ene = 0.0;
		double ave_mz = 0.0;
		double ave_mz_q = 0.0;
		double Cv = 0.0;

		double stderr_ene = 0.0;
		double stderr_mz = 0.0;
		double stderr_mz_q =0.0;
		double stderr_Cv = 0.0;

		double ene_sq = 0.0;
		double mz_sq = 0.0;
		double mz_q_sq = 0.0;

		double ene_fluc_sq = 0.0;
		double ave_ene_fluc = 0.0;
		
		for(; i < mcs_thermalization; i++){
			system.sweep(beta[beta_idx]);
			double ene = system.eval_energy() / n_spins ;
			double mz = system.eval_mz() / n_spins;

		}
		for(int bin_idx=0; bin_idx<n_bins; bin_idx++)
		{
			double ene_m = 0.0;
			double ene_sq_m = 0.0;
			double mz_m = 0.0;
			double mz_sq_m = 0.0;
			double ene_fluc_m = 0.0;
			
			for(int j=0 ; j < n_samples_per_bin ; i++,j++)
			{
				system.sweep(beta[beta_idx]);
				double ene = system.eval_energy() / n_spins ;
				double mz = system.eval_mz() / n_spins;

				ene_m += ene;
				ene_sq_m += ene * ene;
				mz_m += mz;
				mz_sq_m += mz*mz;
			}
			
			ene_m /= (n_samples_per_bin);
			ene_sq_m /= (n_samples_per_bin);
			mz_m /= (n_samples_per_bin);
			mz_sq_m /= (n_samples_per_bin);
			

			ene_sq += ene_m * ene_m;
			mz_sq += mz_m * mz_m;
			mz_q_sq += mz_sq_m* mz_sq_m;

			ave_ene += ene_m;
			ave_mz += mz_m;
			ave_mz_q  += mz_sq_m;

			ene_fluc_m = ene_sq_m - ene_m * ene_m;
			ene_fluc_sq += ene_fluc_m * ene_fluc_m;
			ave_ene_fluc += ene_fluc_m;


			for(int j = 0 ; j < mcs_interval_btwn_bins; i++, j++ ){
				system.sweep(beta[beta_idx]);
				double ene = system.eval_energy() / n_spins ;
				double mz = system.eval_mz() / n_spins;

	
			}
		}
		ene_sq /= n_bins;
		mz_sq /= n_bins;
		ene_fluc_sq /= n_bins;

		ave_ene /= n_bins;
		ave_mz /= n_bins;
		ave_mz_q  /= n_bins;
		ave_ene_fluc /= n_bins;
		Cv = beta[beta_idx] *  beta[beta_idx] * ave_ene_fluc * n_spins;

		const double alpha = 1.96;
		
		stderr_ene = sqrt((ene_sq - ave_ene * ave_ene)/(n_bins-1));
		stderr_mz = sqrt((mz_sq - ave_mz * ave_mz)/(n_bins-1));
		stderr_mz_q = sqrt((mz_q_sq - ave_mz_q * ave_mz_q)/(n_bins-1));
		stderr_Cv = sqrt((ene_fluc_sq - ave_ene_fluc * ave_ene_fluc)/(n_bins-1));

		const double dlt_energy = std::fabs( ave_ene - exact_energy );
		if ( dlt_energy <= t_inv_0975_19 * stderr_ene ) n_in_CI_energy++;
		//std::cout<< ave_ene << "," << exact_energy <<","<<stderr_ene<<std::endl;

		const double dlt_C = std::fabs( Cv - exact_C );
		if ( dlt_C <= t_inv_0975_19 *  stderr_Cv ) n_in_CI_C++;
		//std::cout<< Cv << "," << exact_C <<","<<stderr_Cv<<std::endl;
		
		const double dlt_magz_sq = std::fabs( ave_mz_q - exact_magz_sq );
		if ( dlt_magz_sq <= t_inv_0975_19 * stderr_mz_q ) n_in_CI_magz_sq++;
		//std::cout<< ave_mz_q << "," << exact_magz_sq <<","<<stderr_mz_q<<std::endl;
		
	}
	double portion_in_CI_energy = static_cast<double>(n_in_CI_energy) / n_trials;
	std::cout << "T = " << 1.0 / b_spec << " : Energy : Portion within CI (95%) = " << portion_in_CI_energy << " (stddev for the binomial dist. = " << stddev_binomial << "; tol. parameter = " << dev_tol << ")";
	bool check_energy = (std::fabs(portion_in_CI_energy - CI) / stddev_binomial < dev_tol);
	if (check_energy) {
		std::cout << " ... OK" << std::endl;
	} else {
		std::cout << " (*)" << std::endl;
	}
	
	double portion_in_CI_C = static_cast<double>(n_in_CI_C) / n_trials;
	std::cout << "T = " << 1.0 / b_spec << " : C : Portion within CI (95%) = " << portion_in_CI_C << " (stddev for the binomial dist. = " << stddev_binomial << "; tol. parameter = " << dev_tol << ")";
	bool check_C = (std::fabs(portion_in_CI_C - CI) / stddev_binomial < dev_tol);
	if (check_C) {
		std::cout << " ... OK" << std::endl;
	} else {
		std::cout << " (*)" << std::endl;
	}
	

	double portion_in_CI_magz_sq = static_cast<double>(n_in_CI_magz_sq) / n_trials;
	std::cout << "T = " << 1.0 / b_spec << " : M^2 : Portion within CI (95%) = " << portion_in_CI_magz_sq << " (stddev for the binomial dist. = " << stddev_binomial << "; tol. parameter = " << dev_tol << ")";
	bool check_magz_sq = (std::fabs(portion_in_CI_magz_sq - CI) / stddev_binomial < dev_tol);
	if (check_magz_sq) {
		std::cout << " ... OK" << std::endl;
	} else {
		std::cout << " (*)" << std::endl;
	}

	return check_energy && check_magz_sq;
}

TEST_CASE("CubicLatticeIsingSystem_Stderr_validation") {
	const int n_trials = 100;

	REQUIRE(stderr_validation(2, 0.25, n_trials) == 1);
	
}