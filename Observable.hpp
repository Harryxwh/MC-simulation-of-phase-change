#ifndef Observable_hpp
#define Observable_hpp

#include <cassert>
#include <algorithm>

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
	
class Observable {
protected:
	double q;
	double q_abs;
	double q_sq;
	double q_quar;
	double q_fluc;
	
private:
	double Z_thermal;
	std::string name;
	
public:
	Observable() : q(0.0), q_abs(0.0), q_sq(0.0), q_quar(0.0), q_fluc(0.0), Z_thermal(0.0), name("n/a") {};
	~Observable() {};
	
	double _q() const { return q; };
	double _q_abs() const { return q_abs; };
	double _q_sq() const { return q_sq; };
	double _q_quar() const { return q_quar; };
	double _q_fluc() const { return q_fluc;	};
	double eval_q_fluc() {
		q_fluc = q_sq - q * q;
		return q_fluc;
	};
	void set_q(const double val) { q = val; };
	void set_q_abs(const double val) { q_abs = val; };
	void set_q_sq(const double val) { q_sq = val; };
	void set_q_quar(const double val) { q_quar = val; };
	void set_name(std::string spec) { name = spec; };
	std::string _name() const { return name; };
	
	void update_sampling_by(const double& x) {
		q += x;
		q_abs += std::fabs(x);
		const double sq = x * x;
		q_sq += sq;
		q_quar += sq * sq;
	};

	void update_direct_by(const double& x, const double& weight) {
		q += x * weight;
		q_abs += std::fabs(x) * weight;
		const double sq = x * x;
		q_sq += sq * weight;
		q_quar += sq * sq * weight;
		Z_thermal += weight;
	};
	
	void weight_correction_by(const double& w) {
		q *= w;
		q_abs *= w;
		q_sq *= w;
		q_quar *= w;
		Z_thermal *= w;
	};

	void normalize_by_Z() {
		q /= Z_thermal;
		q_abs /= Z_thermal;
		q_sq /= Z_thermal;
		q_quar /= Z_thermal;
		eval_q_fluc();
	};
	
	void normalize_by_samples(int n_samples) {
		q /= n_samples;
		q_abs /= n_samples;
		q_sq /= n_samples;
		q_quar /= n_samples;
		eval_q_fluc();
	};
};


class MonteCarloObservable : public Observable {
private:
	const int n_bins;
	const int n_samples_per_bin;
	std::vector<Observable> MC_bin_data;
	std::vector<int> n_samples;
	int current_bin_index;
	
	double stderr_q;
	double stderr_q_abs;
	double stderr_q_sq;
	double stderr_q_quar;
	double stderr_q_fluc;
	
public:
	MonteCarloObservable(int n_bins_spec, int n_samples_per_bin_spec) : n_bins(n_bins_spec), 
						n_samples_per_bin(n_samples_per_bin_spec), current_bin_index(0), stderr_q(0), stderr_q_abs(0),
						 stderr_q_sq(0), stderr_q_quar(0), stderr_q_fluc(0.0) {
		MC_bin_data.resize(n_bins);
		n_samples.assign(n_bins, 0);
	};
	~MonteCarloObservable() {};

	void update_sampling_by(const double& x) {
		MC_bin_data[current_bin_index].update_sampling_by(x);
		++(n_samples[current_bin_index]);
	};
	
	bool check_if_bin_filled() {
		return (n_samples[current_bin_index] == n_samples_per_bin);
	};
	
	void switch_bin() {
		MC_normalize_bin();
		current_bin_index++;
	};
	
	void MC_normalize_bin() {
		MC_bin_data[current_bin_index].normalize_by_samples(n_samples_per_bin);
	};
	
	bool check_if_complete_sampling() const {
		bool check = true;
		for (int bin_idx = 0; bin_idx < n_bins; bin_idx++) {
			if (check) check = (n_samples[bin_idx] == n_samples_per_bin);
		}
		return check;
	};

	void MC_normalize() {
		assert(check_if_complete_sampling());
		for (int bin_idx = 0; bin_idx < n_bins; bin_idx++) {
			this->incorporate(MC_bin_data[bin_idx]);
		}
		this->normalize_by_bins();
	}
	
	void print(double b_spec) {
		for (int bin_idx = 0; bin_idx < n_bins; bin_idx++) {
			printf("%s for beta=%f; bin %d : %d samples \n", _name().c_str(), b_spec, bin_idx, n_samples[bin_idx]);
		}
		std::cout << std::endl;
	};

	double _stderr_q() const { return stderr_q; };
	double _stderr_q_abs() const { return stderr_q_abs; };
	double _stderr_q_sq() const { return stderr_q_sq; };
	double _stderr_q_quar() const { return stderr_q_quar; };
	double _stderr_q_fluc() const { return stderr_q_fluc; };

private:
	void incorporate(const Observable& RHS) {
		q += RHS._q();
		q_abs += RHS._q_abs();
		q_sq += RHS._q_sq();
		q_quar += RHS._q_quar();
		q_fluc += RHS._q_fluc();
		
		stderr_q += RHS._q() * RHS._q();
		stderr_q_abs += RHS._q_abs() * RHS._q_abs();
		stderr_q_sq += RHS._q_sq() * RHS._q_sq();
		stderr_q_quar += RHS._q_quar() * RHS._q_quar();
		stderr_q_fluc += RHS._q_fluc() * RHS._q_fluc();
	};
	
	void normalize_by_bins() {
		q /= n_bins;
		q_abs /= n_bins;
		q_sq /= n_bins;
		q_quar /= n_bins;
		q_fluc /= n_bins;

		stderr_q /= n_bins;
		stderr_q = sqrt((stderr_q - q * q) / (n_bins - 1));
		stderr_q_abs /= n_bins;
		stderr_q_abs = sqrt((stderr_q_abs - q_abs * q_abs) / (n_bins - 1));
		stderr_q_sq /= n_bins;
		stderr_q_sq = sqrt((stderr_q_sq - q_sq * q_sq) / (n_bins - 1));
		stderr_q_quar /= n_bins;
		stderr_q_quar = sqrt((stderr_q_quar - q_quar * q_quar) / (n_bins - 1));
		stderr_q_fluc /= n_bins;
		stderr_q_fluc = sqrt((stderr_q_fluc - q_fluc * q_fluc) / (n_bins - 1));
	};
};


class DataBundle {
protected:
	const int n_spins;
	std::vector<double> beta;
	const int n_bins;
	const int n_samples_per_bin;
	std::vector<MonteCarloObservable> energy;
	std::vector<Observable> exact_energy;
	std::vector<MonteCarloObservable> magz;
	std::vector<Observable> exact_magz;
	
public:
	DataBundle(int n_spins_spec, const std::vector<double>& beta_spec, int n_bins_spec = 0, int n_samples_per_bin_spec = 0) : n_spins(n_spins_spec), beta(DataBundle::eval_sorted_beta(beta_spec)), n_bins(n_bins_spec), n_samples_per_bin(n_samples_per_bin_spec),
	energy(beta_spec.size(), MonteCarloObservable(n_bins_spec, n_samples_per_bin_spec)), exact_energy(beta_spec.size()),
	magz(beta_spec.size(), MonteCarloObservable(n_bins_spec, n_samples_per_bin_spec)), exact_magz(beta_spec.size()) {
		for (auto& each: energy) each.set_name("energy");
		for (auto& each: exact_energy) each.set_name("energy (exact)");
		for (auto& each: magz) each.set_name("magnetization");
		for (auto& each: exact_magz) each.set_name("magnetization (exact)");
	};
	
	virtual ~DataBundle() {};
	
	std::vector<double> _beta() const { return beta; };
	
	double _beta(int beta_idx) const { return beta[beta_idx]; };
	
	double eval_T(int beta_idx) const { return 1.0 / beta[beta_idx]; };
	
	void update_sampling_energy(int beta_idx, const double& val) {
		energy[beta_idx].update_sampling_by(val);
	};
	
	void update_sampling_magz(int beta_idx, const double& val) {
		magz[beta_idx].update_sampling_by(val);
	};
	
	void MC_normalize(int beta_idx) {
		assert(check_if_complete_sampling_beta_index(beta_idx));
		energy[beta_idx].MC_normalize();
		magz[beta_idx].MC_normalize();
	};
	
	bool check_if_complete_sampling_beta_index(int beta_idx) const {
		bool check = true;
		if ( check ) check = energy[beta_idx].check_if_complete_sampling();
		if ( check ) check = magz[beta_idx].check_if_complete_sampling();
		return check;
	};
	
	bool check_if_bin_filled(int sampling_current_bin) {
		return sampling_current_bin == n_samples_per_bin;
	};

	void switch_bin(int beta_idx) {
		energy[beta_idx].switch_bin();
		magz[beta_idx].switch_bin();
	};
	
	void weight_correction(int beta_idx, const double& w_correction) {
		exact_energy[beta_idx].weight_correction_by(w_correction);
		exact_magz[beta_idx].weight_correction_by(w_correction);
	};

	void update_exact_energy(int beta_idx, const double& val, const double& weight) {
		exact_energy[beta_idx].update_direct_by(val, weight);
	};
	
	void update_exact_magz(int beta_idx, const double& mz, const double& weight) {
		exact_magz[beta_idx].update_direct_by(mz, weight);
	};

	void normalize_by_Z() {
		for ( unsigned int beta_idx = 0; beta_idx < beta.size(); beta_idx++ ) {
			exact_energy[beta_idx].normalize_by_Z();
			exact_magz[beta_idx].normalize_by_Z();
		}
	};

	double _get_exact_energy(int beta_idx) const {
		return exact_energy.at(beta_idx)._q();
	};
	double _get_exact_energy_per_spin(int beta_idx) const {
		return _get_exact_energy(beta_idx) / n_spins;
	};
	double _get_obs_energy(int beta_idx) const {
		return energy.at(beta_idx)._q();
	};
	double _get_obs_energy_per_spin(int beta_idx) const {
		return _get_obs_energy(beta_idx) / n_spins;
	};
	double _get_stderr_energy(int beta_idx) const {
		return energy.at(beta_idx)._stderr_q();
	};
	double _get_stderr_energy_per_spin(int beta_idx) const {
		return _get_stderr_energy(beta_idx) / n_spins;
	};
	
	double _get_exact_C(int beta_idx) const {
		return beta[beta_idx] * beta[beta_idx] * exact_energy.at(beta_idx)._q_fluc();
	};
	double _get_exact_C_per_spin(int beta_idx) const {
		return _get_exact_C(beta_idx) / n_spins;
	};
	double _get_obs_C(int beta_idx) const {
		return beta[beta_idx] * beta[beta_idx] * energy.at(beta_idx)._q_fluc();
	};
	double _get_obs_C_per_spin(int beta_idx) const {
		return _get_obs_C(beta_idx) / n_spins;
	};
	double _get_stderr_C(int beta_idx) const {
		return beta[beta_idx] * beta[beta_idx] * energy.at(beta_idx)._stderr_q_fluc();
	};
	double _get_stderr_C_per_spin(int beta_idx) const {
		return _get_stderr_C(beta_idx) / n_spins;
	};

	double _get_exact_magz(int beta_idx) const {
		return exact_magz.at(beta_idx)._q();
	};
	double _get_exact_magz_per_spin(int beta_idx) const {
		return _get_exact_magz(beta_idx) / n_spins;
	};
	double _get_obs_magz(int beta_idx) const {
		return magz.at(beta_idx)._q();
	};
	double _get_obs_magz_per_spin(int beta_idx) const {
		return _get_obs_magz(beta_idx) / n_spins;
	};
	double _get_stderr_magz(int beta_idx) const {
		return magz.at(beta_idx)._stderr_q();
	};
	double _get_stderr_magz_per_spin(int beta_idx) const {
		return _get_stderr_magz(beta_idx) / n_spins;
	};

	double _get_exact_magz_sq(int beta_idx) const {
		return exact_magz.at(beta_idx)._q_sq();
	};
	double _get_exact_magz_sq_per_spin(int beta_idx) const {
		return _get_exact_magz_sq(beta_idx) / eval_n_spins_sq();
	};
	double _get_obs_magz_sq(int beta_idx) const {
		return magz.at(beta_idx)._q_sq();
	};
	double _get_obs_magz_sq_per_spin(int beta_idx) const {
		return _get_obs_magz_sq(beta_idx) / eval_n_spins_sq();
	};
	double _get_stderr_magz_sq(int beta_idx) const {
		return magz.at(beta_idx)._stderr_q_sq();
	};
	double _get_stderr_magz_sq_per_spin(int beta_idx) const {
		return _get_stderr_magz_sq(beta_idx) / eval_n_spins_sq();
	};
	
	void output_legends_MC(const std::string system_size) {
		std::cout << system_size << std::endl;
		std::cout << "# T[1] E[2] err[3] C[4] err[5] M[6] err[7] M^2[8] err[9]";
	};
	
	void output_data_MC(int beta_idx) {
		std::cout << eval_T(beta_idx) << ' '
		<< _get_obs_energy_per_spin(beta_idx) << ' '
		<< _get_stderr_energy_per_spin(beta_idx) << ' '
		<< _get_obs_C_per_spin(beta_idx) << ' '
		<< _get_stderr_C_per_spin(beta_idx) << ' '
		<< _get_obs_magz_per_spin(beta_idx) << ' '
		<< _get_stderr_magz_per_spin(beta_idx) << ' '
		<< _get_obs_magz_sq_per_spin(beta_idx) << ' '
		<< _get_stderr_magz_sq_per_spin(beta_idx) << ' ';
	};
	
	void output_legends_exact(const std::string system_size) {
		std::cout << system_size << std::endl;
		std::cout << "# T[1] E[2] C[3] M[4] M^2[5]";
	};
	
	void output_data_exact(int beta_idx) {
		std::cout << eval_T(beta_idx) << ' '
		<< _get_exact_energy_per_spin(beta_idx) << ' '
		<< _get_exact_C_per_spin(beta_idx) << ' '
		<< _get_exact_magz_per_spin(beta_idx) << ' '
		<< _get_exact_magz_sq_per_spin(beta_idx) << ' ';
	};
	
	int _mcs_sampling() const { return n_samples_per_bin * n_bins; };
	
	int eval_n_spins_sq() const {	return n_spins * n_spins;	};

protected:
	static std::vector<double> eval_sorted_beta(std::vector<double> b) {
		std::sort(b.begin(), b.end());
		return b;
	};
};


class CubicLatticeDataBundle : public DataBundle {
	friend class CubicLatticeIsingSystem; //can access private members of ...

public:
	CubicLatticeDataBundle(int n_spins_spec, const std::vector<double>& beta_spec) : DataBundle(n_spins_spec, beta_spec) {};
	
	~CubicLatticeDataBundle() {};
	
	void output_legends_exact() {
		/*
		std::ofstream ofs;
		ofs.open(outf, std::ios::out | std::ios::trunc);  
		ofs << system_size << std::endl;
		ofs << "T[1],E[2],C[3],M[4],M^2[5]" << std::endl;
		*/
		//std::cout << system_size << std::endl;
		std::cout << "T[1],E[2],C[3],M[4],M^2[5]" << std::endl;
	};
	

	void output_data_exact_to_file(int beta_idx,std::string outf = "exact_out.csv") {
		std::ofstream ofs;
		ofs.open(outf, std::ios::out | std::ios::app); 
		ofs << eval_T(beta_idx) << ',';
		ofs << _get_exact_energy_per_spin(beta_idx) << ',';
		ofs << _get_exact_C_per_spin(beta_idx) << ',';
		ofs << _get_exact_magz_per_spin(beta_idx) << ',';
		ofs << _get_exact_magz_sq_per_spin(beta_idx) << ',';
		ofs << std::endl;
		ofs.close();
	};

	void output_data_exact(int beta_idx) {
		std::cout << eval_T(beta_idx) << ',';
		std::cout << _get_exact_energy_per_spin(beta_idx) << ',';
		std::cout << _get_exact_C_per_spin(beta_idx) << ',';
		std::cout << _get_exact_magz_per_spin(beta_idx) << ',';
		std::cout << _get_exact_magz_sq_per_spin(beta_idx) << ',';
		std::cout << std::endl;
	};
	
	void output_data_exact_all() {
		for ( unsigned int beta_idx = 0; beta_idx < beta.size(); beta_idx++ ) {
			output_data_exact(beta_idx);
		}
	};

	
};
#endif /* Observable_hpp */