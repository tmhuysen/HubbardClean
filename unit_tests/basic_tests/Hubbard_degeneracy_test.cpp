//
// Created by Wulfix on 20/11/2017.
//


#include "gtest/gtest.h"
#include "include/io.h"
#include "include/Hubbard.h"

TEST(io_test_3, test){
    // For four sites, an example hopping matrix can be found in the following file
    const std::string filename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/unit_tests/basic_tests/input_data/hopping_matrix4.data";
    arma::mat H_read = read_hopping_matrix_from_file(filename, ',');


    // Check if the wanted result is equal to the calculated one
    std::cout<<std::endl<<H_read<<std::endl;
    Lattice lat = Lattice(H_read);
    Hubbard hubby = Hubbard(4,lat);
    for(int i = 0;i<5;i++){
        Eigen::MatrixXd H = hubby.getSector(i);
        Eigen::EigenSolver<Eigen::MatrixXd> solver(H);
        Eigen::VectorXd eigenvalues = solver.eigenvalues().real().cast<double>();
        Eigen::MatrixXd eigenvectors = solver.eigenvectors().real().cast<double>();
        std::cout<<std::endl<<eigenvalues;

    }






}


