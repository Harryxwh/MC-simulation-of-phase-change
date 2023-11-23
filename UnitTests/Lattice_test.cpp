#include "catch2/catch_amalgamated.hpp"
#include "../Lattice.hpp"


/*
int CubicLattice_N_sites(int L0, int L1, int L2) {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = L0;
	L[1] = L1;
	L[2] = L2;
	CubicLattice lattice(L);
	return lattice._N_sites();
}

int CubicLattice_NN(int L0, int L1, int L2, int site_idx, int bond_idx) {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = L0;
	L[1] = L1;
	L[2] = L2;
	CubicLattice lattice(L);
	return lattice._NN_of_Site(site_idx, bond_idx);
}
*/

/*
TEST_CASE("CubicLattice") {
	int N_SL_spec = 2;  
	int z_spec[2] = {3,3};
	const std::vector<int> L_spec = {3,3};
	Lattice Lat(N_SL_spec,z_spec, L_spec);
	
	REQUIRE(CubicLattice_N_sites(3,2,3) == 18);
	REQUIRE(CubicLattice_NN(3, 2, 3, 7, 0) == 8);
	REQUIRE(CubicLattice_NN(3, 2, 3, 7, 1) == 10);
	REQUIRE(CubicLattice_NN(3, 2, 3, 7, 2) == 13);
	REQUIRE(CubicLattice_NN(3, 2, 3, 7, 3) == 6);
	REQUIRE(CubicLattice_NN(3, 2, 3, 7, 4) == 10);
	REQUIRE(CubicLattice_NN(3, 2, 3, 7, 5) == 1);
	
}
*/
