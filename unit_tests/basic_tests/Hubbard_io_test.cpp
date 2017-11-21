//
// Created by Wulfix on 20/11/2017.
//


#include "gtest/gtest.h"
#include "HubbardClass/include/io.h"
#include "HubbardClass/include/Hubbard.h"

TEST(io_test, test){
    // For four sites, an example hopping matrix can be found in the following file
    const std::string filename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/unit_tests/basic_tests/input_data/hopping_matrix2.data";
    arma::mat H_read = read_hopping_matrix_from_file(filename, ' ');

    // Specify the wanted result

    arma::mat H =       {{ 0,  -1,    0,  -1},
                        {-1,   0,   -1,   0},
                        { 0,  -1,    0,  -1},
                        {-1,   0,   -1,   0}};


    // Check if the wanted result is equal to the calculated one
    std::cout<<std::endl<<H_read<<std::endl;
    Lattice lat = Lattice(H_read);
    Hubbard hubby = Hubbard(4,lat);
    hubby.print();





}


