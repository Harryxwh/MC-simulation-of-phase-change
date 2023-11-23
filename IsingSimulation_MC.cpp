#include "Lattice.hpp"
#include "IsingSystem.hpp"
#include <vector>
#include <string>
#include <cmath>
#include "pugixml/pugixml.hpp"
#include "ParameterBundle.hpp"

using namespace std;
using namespace pugi;


void run_MC( std::vector<int> L_spec, std::vector<double> beta, const int N_TemperaturePoints, const int seed,
const int n_bins,const int mcs_thermalization, const int n_samples_per_bin,const int mcs_interval_btwn_bins,
bool output_data = true, string outf = "MC_data_out.csv",string outf_est = "MC_est_out.csv"){

	
	//Size
	const int dim = 3;
	std::vector<int> L = L_spec;
	const int n_spins = CubicLatticeIsingSystem::eval_n_spins(L);

	CubicLatticeDataBundle data_bundle(n_spins, beta);
	CubicLatticeIsingSystem system(L, data_bundle, seed);

	//const bool need_exact = (n_spins <= 32);
	const bool need_exact = 0;
	//Exact Count
	if ( need_exact ){
	system.exact_count();
	data_bundle.output_legends_exact();
	data_bundle.output_data_exact_all();
	std::cout<<"++++++Exact count finished++++++"<<endl<<endl;
	}

	//Monte Carlo

	std::ofstream ofs;
	std::ofstream ofs_est;
	ofs.open(outf, std::ios::out | std::ios::trunc); 
	ofs_est.open(outf_est, std::ios::out | std::ios::trunc);

	if(output_data) ofs<<"step,"<<"E_per_spin,"<<"Mz_per_spin,"<<"Mz_sq_perspin"<<endl;
	ofs_est<<"T,"<<"ave_ene,"<<"sigma_ene,"<<"ave_mz,"<<"sigma_mz,"<<"ave_mz_sq,"
				<<"sigma_mz_sq,"<<"Cv,"<<"sigma_Cv,"<<"ave_mz_qua,"<<"Binder_ratio,"<<endl;
	
	for(int beta_idx =0 ; beta_idx < N_TemperaturePoints; beta_idx++)
	//for(int beta_idx = N_TemperaturePoints-1 ; beta_idx >= 0; beta_idx--)
	{	
		system.CubicHotInitialize();
		int i = 0;

		double ave_ene = 0.0; // <E>
		double ave_mz = 0.0; // <M>
		double ave_mz_q = 0.0; // <M^2>
		double ave_mz_qua = 0.0;
		double Cv = 0.0;

		double ene_sq = 0.0; //<E(mean within a bin)^2>
		double mz_sq = 0.0;  //<M(mean within a bin)^2> 
		double mz_q_sq = 0.0; //<M(mean within a bin)^4>
		
		 
		// the two bunches of variables above are used to compute stderr
		double stderr_ene = 0.0;
		double stderr_mz = 0.0;
		double stderr_mz_q =0.0;
		double stderr_Cv = 0.0;

		double ene_fluc_sq = 0.0; // < (<E^2> - <E>^2)^2 >  (<E_fluc^2>)
		double ave_ene_fluc = 0.0;  // < <E^2> - <E>^2 >
		
		for(; i < mcs_thermalization; i++){
			system.sweep(beta[beta_idx]);
			double ene = system.eval_energy() / n_spins ;
			double mz = system.eval_mz() / n_spins;

			if(output_data)  ofs<<i<<","<<ene<<","<<mz<<","<<mz*mz<<endl;
		}
		for(int bin_idx=0; bin_idx<n_bins; bin_idx++)
		{
			//m means for one bin
			double ene_m = 0.0;
			double ene_sq_m = 0.0;
			double mz_qua_m = 0.0;
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
				mz_qua_m += mz * mz * mz * mz;
				
				if(output_data)  ofs<<i<<","<<ene<<","<<mz<<","<<mz*mz<<endl;
			}
			ene_m /= (n_samples_per_bin);
			ene_sq_m /= (n_samples_per_bin);
			mz_m /= (n_samples_per_bin);
			mz_sq_m /= (n_samples_per_bin);
			mz_qua_m /= (n_samples_per_bin);

			ene_sq += ene_m * ene_m;
			mz_sq += mz_m * mz_m;
			mz_q_sq += mz_sq_m * mz_sq_m;

			ave_ene += ene_m;
			ave_mz += mz_m;
			ave_mz_q  += mz_sq_m;
			ave_mz_qua += mz_qua_m;

			ene_fluc_m = ene_sq_m - ene_m * ene_m;
			ene_fluc_sq += ene_fluc_m * ene_fluc_m;
			ave_ene_fluc += ene_fluc_m;


			for(int j = 0 ; j < mcs_interval_btwn_bins; i++, j++ ){
				system.sweep(beta[beta_idx]);
				double ene = system.eval_energy() / n_spins ;
				double mz = system.eval_mz() / n_spins;

				if(output_data)  ofs<<i<<","<<ene<<","<<mz<<","<<mz*mz<<endl;
			}
		}
		ene_sq /= n_bins;
		mz_sq /= n_bins;
		ene_fluc_sq /= n_bins;

		ave_ene /= n_bins;
		ave_mz /= n_bins;
		ave_mz_q  /= n_bins;
		ave_mz_qua /= n_bins;
		ave_ene_fluc /= n_bins;
		Cv = beta[beta_idx] *  beta[beta_idx] * ave_ene_fluc * n_spins;

		const double alpha = 1.96;
		
		stderr_ene = sqrt((ene_sq - ave_ene * ave_ene)/(n_bins-1));
		stderr_mz = sqrt((mz_sq - ave_mz * ave_mz)/(n_bins-1));
		stderr_mz_q = sqrt((mz_q_sq - ave_mz_q * ave_mz_q)/(n_bins-1));
		stderr_Cv = sqrt((ene_fluc_sq - ave_ene_fluc * ave_ene_fluc)/(n_bins-1));

		double Binder_ratio = ave_mz_qua / (ave_mz_q * ave_mz_q);		
		ofs_est << 1/beta[beta_idx] <<"," << ave_ene << "," << alpha * stderr_ene << "," 
								<< ave_mz << "," << alpha  * stderr_mz << "," <<ave_mz_q << "," << alpha  * stderr_mz_q 
								<< "," << Cv << "," << alpha* stderr_Cv <<","<< ave_mz_qua <<"," << Binder_ratio << endl;
		std::cout<<" T = "<<1/beta[beta_idx]<<" finished"<<endl;
	}
	ofs.close();
}



int main(int argc,char *argv[]){
	/* Cubic lattice system (brute-force exact counting) */
	string input_XML = (argc>1) ? argv[1] : "../InputParam_MC.xml";
	
	CubicLatticeParameterBundle Parameters(input_XML);
	
	//Size
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = Parameters._L0();
	L[1] = Parameters._L1();
	L[2] = Parameters._L2();
	const int n_spins = CubicLatticeIsingSystem::eval_n_spins(L);

	//Temperatures
	double Tupb = Parameters._Tmax();
	double Tmin = Parameters._Tmin();
	int N_TemperaturePoints = Parameters._N_TemperaturePoints();
	const double Tdlt = (Tupb - Tmin) / N_TemperaturePoints;
	std::vector<double> beta(N_TemperaturePoints);
	for (int i = 0; i < N_TemperaturePoints; i++) beta[i] = 1.0 / (Tupb - i * Tdlt);

	//Seed
	int seed = Parameters._Seed();

	//MC parameters
	const int n_bins = Parameters._n_bins();
	const int mcs_thermalization = Parameters._mcs_thermalization();
	const int n_samples_per_bin = Parameters._n_samples_per_bin();
	const int mcs_interval_btwn_bins =  Parameters._mcs_interval_btwn_bins();

	const bool output_data = false;
	run_MC(L,beta,N_TemperaturePoints,seed,n_bins,mcs_thermalization,n_samples_per_bin,mcs_interval_btwn_bins,output_data);

	std::cout<<"++++++++MC finished++++++++++"<<endl;

	return 0;
}