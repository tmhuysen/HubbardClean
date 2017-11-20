cmake_minimum_required(VERSION 3.8)
project(Hubbard_Head)

set(CMAKE_CXX_STANDARD 11)

add_library(Hubbard_Head HubbardClass/include/bitset.h test.cpp HubbardClass/src/HubbardSolver.cpp HubbardClass/src/io.cpp HubbardClass/include/io.h HubbardClass/src/OneRDM.cpp HubbardClass/include/OneRDM.h)
#file(GLOB helloworld_SRC
#        "*.h"
#        "*.cpp"
#        )


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