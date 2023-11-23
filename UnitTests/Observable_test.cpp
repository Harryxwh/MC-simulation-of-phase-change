#include "catch2/catch_amalgamated.hpp"
#include "../Observable.hpp"

/*
TEST_CASE("Observable Initialization") {
	Observable test_obs;
	REQUIRE(test_obs._q() == 0);
	REQUIRE(test_obs._q_abs() == 0);
	REQUIRE(test_obs._q_sq() == 0);
	REQUIRE(test_obs._q_quar() == 0);
	
	test_obs.update_direct_by(10, 1.0);
	REQUIRE(test_obs._q() == 10);
	REQUIRE(test_obs._q_abs() == 10);
	REQUIRE(test_obs._q_sq() == 100);
	REQUIRE(test_obs._q_quar() == 10000);
	
	test_obs.update_direct_by(-10, 1.0);
	REQUIRE(test_obs._q() == 0);
	REQUIRE(test_obs._q_abs() == 20);
	REQUIRE(test_obs._q_sq() == 200);
	REQUIRE(test_obs._q_quar() == 20000);
	
	test_obs.normalize_by_Z();
	REQUIRE(test_obs._q() == 0);
	REQUIRE(test_obs._q_abs() == 10);
	REQUIRE(test_obs._q_sq() == 100);
	REQUIRE(test_obs._q_quar() == 10000);
}
*/