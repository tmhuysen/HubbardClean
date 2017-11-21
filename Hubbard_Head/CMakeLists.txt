cmake_minimum_required(VERSION 3.8)
project(Hubbard_Head)

set(CMAKE_CXX_STANDARD 11)

#add_library(Hubbard_Head HubbardClass/include/bitset.h test.cpp HubbardClass/src/HubbardSolver.cpp HubbardClass/src/io.cpp HubbardClass/include/io.h HubbardClass/src/OneRDM.cpp HubbardClass/include/OneRDM.h Hubbard_Head/HubbardClass/src/Hubbard.cpp Hubbard_Head/HubbardClass/include/Hubbard.h unit_tests/basic_tests/test.cpp Hubbard_Head/HubbardClass/src/Lattice.cpp Hubbard_Head/HubbardClass/include/Lattice.h Hubbard_Head/HubbardClass/src/bitset.cpp Hubbard_Head/HubbardClass/src/AddressingMatrix.cpp Hubbard_Head/HubbardClass/include/AddressingMatrix.h Hubbard_Head/HubbardClass/src/SpinSector.cpp)







file(GLOB SOURCES "HubbardClass/src/*.cpp" "HubbardClass/include/*.h")


add_library(Hubbard_Head SHARED ${SOURCES} test.cpp test.h)

find_package(eigen3)
if (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIRS})
    target_link_libraries (Hubbard_Head_run ${EIGEN3_LIBRARIES})
endif ()

find_package(Boost)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries (Hubbard_Head_run ${Boost_LIBRARIES})
endif()

find_package(Armadillo)
if (Armadillo_FOUND)
    include_directories(${Armadillo_INCLUDE_DIR})
    target_link_libraries (Hubbard_Head_run ${Armadillo_LIBRARIES})
endif()