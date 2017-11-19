cmake_minimum_required(VERSION 3.8)
project(HubbardClean)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp Hubbard_Head/HubbardClass/src/Hubbard.cpp Hubbard_Head/HubbardClass/include/Hubbard.h unit_tests/basic_tests/test.cpp Hubbard_Head/HubbardClass/src/Lattice.cpp Hubbard_Head/HubbardClass/include/Lattice.h Hubbard_Head/HubbardClass/src/bitset.cpp Hubbard_Head/HubbardClass/src/AddressingMatrix.cpp Hubbard_Head/HubbardClass/include/AddressingMatrix.h Hubbard_Head/HubbardClass/src/SpinSector.cpp)
add_executable(Hubbard_Head_run ${SOURCE_FILES})

include_directories(Hubbard_Head)

add_subdirectory(Hubbard_Head)
add_subdirectory(unit_tests)



target_link_libraries(Hubbard_Head_run Hubbard_Head)


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