#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>
#include "Lattice.hpp"
#include "Observable.hpp"
#include "RandomNumberGenerator.hpp"

class IsingSpin {
private:
	int sz; /* +/- 1 */

public:
	IsingSpin() {
		sz = 1;
	};
	~IsingSpin() {};

	int _sz() const {
		return sz;
		
	};

	void set_sz(int sz_spec) {
		assert(sz_spec == 1 || sz_spec == -1);
		sz = sz_spec;
	};

	void flip() {
		sz *= -1;
	};
};

class IsingSystem {
protected:
	const double J;
	const int n_spins;
	std::vector<IsingSpin> spin;
	const long long maxrep_state;
	
public:
	IsingSystem(const int n_spins_spec)
	: J(-1.0),
	n_spins(n_spins_spec),
	maxrep_state(static_cast<long long>(std::pow(2, n_spins)) - 1)
	{
		spin.resize(n_spins);//allocating space and creating instances 
	};

	virtual ~IsingSystem() {};

	double _J() const {
		return J;
	};
	
	int _n_spins() const {
		return n_spins;
	};
	
	int eval_n_spins_sq() const {
		return n_spins * n_spins;
	};

	int _sz(const int site_idx) const {
		return spin[site_idx]._sz();
	}
	
	void flip_spin(const int site_idx) {
		spin[site_idx].flip();
	};
	
	void set_spin(const int site_idx, int s_spec) {
		spin[site_idx].set_sz(s_spec);
	};
	
	double eval_mz() const {
		int mz_tmp = 0;
		for (int site_idx = 0; site_idx < n_spins; site_idx++) {
			mz_tmp += spin[site_idx]._sz();
		}
		return static_cast<double>(mz_tmp);
	};
	
	long long _maxrep_state() const { return maxrep_state; };
	
	void set_state_by_code(long long rep_state) {
		assert(rep_state >= 0 && rep_state <= maxrep_state);
		for (int site_idx = 0; site_idx < n_spins; site_idx++) {
			set_spin(site_idx, (rep_state % 2) ? 1 : -1);
			rep_state /= 2;
		}
	}
};


class CubicLatticeIsingSystem : public IsingSystem {
private:
	CubicLattice lattice;
	CubicLatticeDataBundle* ptrMCDB;
	const int beta_list_size;
	RandomNumberGenerator mtwist;


public:
	CubicLatticeIsingSystem(const std::vector<int>& L_spec, CubicLatticeDataBundle& MCDB_spec, const int seed)
	: IsingSystem(CubicLatticeIsingSystem::eval_n_spins(L_spec)),
	lattice(L_spec), ptrMCDB(&MCDB_spec), beta_list_size(static_cast<unsigned int>(MCDB_spec._beta().size())),
	mtwist(seed,lattice._N_sites()) {};
	
	CubicLatticeIsingSystem(const std::vector<int>& L_spec)
	: IsingSystem(CubicLatticeIsingSystem::eval_n_spins(L_spec)), lattice(L_spec), ptrMCDB(nullptr), beta_list_size(0), 
									mtwist(12345678,CubicLatticeIsingSystem::eval_n_spins(L_spec)) {};
	
	~CubicLatticeIsingSystem() { ptrMCDB = nullptr; };

	int _L0() const { return lattice._L0(); };
	int _L1() const { return lattice._L1(); };
	int _L2() const { return lattice._L2(); };
	std::string _system_size() const { return "# System size : " + std::to_string(_L0()) + " x " + std::to_string(_L1()); };
	

	double _beta(int beta_idx) const {
		assert(ptrMCDB != (DataBundle* )(NULL));
		assert(beta_idx >= 0 && beta_idx < beta_list_size);
		return ptrMCDB->_beta(beta_idx);
	};

	double eval_energy() const {
		double energy_tmp = 0;
		
		const int z_half =  lattice._z_common_half();
		for (int site_idx = 0; site_idx < lattice._N_sites(); site_idx++) {
			for (int bond_index = 0; bond_index < z_half; bond_index++) {
				const int j = lattice._NN_of_Site(site_idx, bond_index);
				energy_tmp += J * spin[site_idx]._sz() * spin[j]._sz();
			}
		}
		
		return energy_tmp;
	};

	void exact_count() {
		assert(ptrMCDB != (DataBundle* )(NULL));
		double ene, mz, weight;
		double ene_min = 0;
		#pragma omp parallel for shared(ptrMCDB, ene_min)
		for(long long rep_state = 0;rep_state <= maxrep_state;rep_state++) {
			//std::cout<<"num_threads = "<<omp_get_num_threads()<<" rep_state="<<rep_state<<std::endl;
			set_state_by_code(rep_state);
			ene = eval_energy();
			mz = eval_mz();

			if (ene < ene_min) { 
				const double dE = ene_min - ene;
				for ( int beta_idx = 0; beta_idx < beta_list_size; beta_idx++ ) {
					const double w_correction = std::exp(-ptrMCDB->_beta(beta_idx) * dE);
					ptrMCDB->exact_energy[beta_idx].weight_correction_by(w_correction);
					ptrMCDB->exact_magz[beta_idx].weight_correction_by(w_correction);
				}
				ene_min = ene;
			}

			for ( int beta_idx = 0; beta_idx < beta_list_size; beta_idx++ ) {
				weight = std::exp(-ptrMCDB->_beta(beta_idx) * (ene - ene_min));
				ptrMCDB->exact_energy[beta_idx].update_direct_by(ene, weight);
				ptrMCDB->exact_magz[beta_idx].update_direct_by(mz, weight);
			}

		}
		
		for ( int beta_idx = 0; beta_idx < beta_list_size; beta_idx++ ) {
			ptrMCDB->exact_energy[beta_idx].normalize_by_Z();
			ptrMCDB->exact_magz[beta_idx].normalize_by_Z();
		}
	};

	static int eval_n_spins(const std::vector<int>& L_spec) {
		int res = 1;
		for(int i = 0; i < static_cast<int>(L_spec.size()) ; i++)
		{
			res *= L_spec.at(i);
		}
		return res;
	};

	void CubicHotInitialize()
	{
		for (int site_idx = 0; site_idx < lattice._N_sites(); site_idx++) {
			double rd = mtwist.gen_rand01();
			if(rd>0.5) set_spin(site_idx, 1);
			else set_spin(site_idx, -1);
		}
	}

	int FlipAttemptUsingMetropolis(int site_idx, double beta,double R)
	{
		double energy_tmp = 0;
		for (int bond_index = 0; bond_index <lattice._z_common(); bond_index++) {
				const int j = lattice._NN_of_Site(site_idx, bond_index);
				energy_tmp += J * spin[site_idx]._sz() * spin[j]._sz();
			}
		double deltaE = -2 * energy_tmp;
		double p = (std::exp(-deltaE * beta)>1) ? 1 : std::exp(-deltaE * beta);
		if (p > R) return 1;
		return 0;
	}

	int FlipAttemptUsingHeatBath(int site_idx, double beta,double R)
	{
		double energy_tmp = 0;
		for (int bond_index = 0; bond_index <lattice._z_common(); bond_index++) {
				const int j = lattice._NN_of_Site(site_idx, bond_index);
				energy_tmp += J * spin[site_idx]._sz() * spin[j]._sz();
			}
		
		double p = std::exp(beta * energy_tmp)/(std::exp(-beta * energy_tmp)+std::exp(beta * energy_tmp));
		if (p > R) return 1;
		return 0;
	}

	void sweep(double beta)
	{	
		
		for (int i = 0 ;i < lattice._N_sites();i++)
		{
			int site_idx = mtwist.gen_rand_site();
			double R = mtwist.gen_rand01();
			if(FlipAttemptUsingMetropolis(site_idx,beta,R))
			{
				flip_spin(site_idx);
			}
		}
	}

	void output_observables_tofile(int MCstep,std::string outf){
	
		std::ofstream ofs;
		ofs.open(outf, std::ios::out | std::ios::app); 
		double ene = eval_energy()/n_spins;
		double mz = eval_mz()/n_spins;
		ofs << MCstep <<',';
		ofs << ene << ',';
		ofs << mz << ',';
		ofs << std::endl;
	} 

};
