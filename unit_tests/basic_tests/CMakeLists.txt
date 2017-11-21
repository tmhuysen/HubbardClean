include_directories(${googletest_SOURCE_DIR}/include ${googletest_SOURCE_DIR})

add_executable(basic_tests
      test.cpp OneRDM_test.cpp phaseCheck_test.cpp Hubbard_io_test.cpp AddressingMatrix_test.cpp output_solver_hubbard_rdm_test.cpp)

target_link_libraries(basic_tests gtest gtest_main)
target_link_libraries(basic_tests Hubbard_Head)

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