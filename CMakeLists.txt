cmake_minimum_required(VERSION 3.14)
project(ProjectIsing CXX)
set(CMAKE_CXX_STANDARD 14)

### Main target 
file( GLOB SRCS *.cpp *.hpp )
###add_executable(IsingSimulation ${SRCS})
add_executable(IsingSimulation_MC_Prototype IsingSimulation_MC_Prototype.cpp 
IsingSystem.hpp IsingSystem.cpp 
Lattice.hpp Lattice.cpp 
Observable.hpp Observable.cpp 

pugixml/pugixml.cpp)
target_compile_options(IsingSimulation_MC_Prototype PRIVATE -Wall -Wextra)

file( GLOB SRCS *.cpp *.hpp )
###add_executable(IsingSimulation ${SRCS})
add_executable(IsingSimulation_MC IsingSimulation_MC.cpp 
IsingSystem.hpp IsingSystem.cpp 
Lattice.hpp Lattice.cpp 
Observable.hpp Observable.cpp 

pugixml/pugixml.cpp)
target_compile_options(IsingSimulation_MC PRIVATE -Wall -Wextra)

### Unit test
enable_testing()
file( GLOB SRCS_Catch2 UnitTests/catch2/*.hpp UnitTests/catch2/*.cpp)
add_executable(test_IsingSimulation IsingSystem.hpp IsingSystem.cpp Lattice.hpp Lattice.cpp Observable.hpp Observable.cpp 
UnitTests/ParameterBundle_test.cpp 
UnitTests/MCtest.cpp 
UnitTests/RandomNumberGenerator_test.cpp 
pugixml/pugixml.cpp  ${SRCS_Catch2})
target_compile_options(test_IsingSimulation PRIVATE -Wall -Wextra)
add_test(test_IsingSimulation test_IsingSimulation)

### Load Catch2 from remote (commented out)
#include(FetchContent)
#FetchContent_Declare(catch2
#GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#GIT_TAG v3.1.0
#)
#FetchContent_MakeAvailable(Catch2)
#target_link_libraries(test_IsingSimulation Catch2::Catch2WithMain)

