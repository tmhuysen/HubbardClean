//
// Created by Wulfix on 21/11/2017.
//

#include "gtest/gtest.h"
#include "include/Hubbard.h"
#include "include/OneRDM.h"

TEST(setLattice, test_visual_cout){

    const std::string filename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/unit_tests/basic_tests/input_data/hopping_matrix3.data";
    arma::mat H_read = read_hopping_matrix_from_file(filename, ' ');

    // Specify the wanted result

    arma::mat H =       {{ 0.123,  1,    0,  1},
                         {1,   0.5,   1,   0},
                         { 0,  1,    99.5,  1},
                         {1,   0,   1,   0.123}};

    int N = 4;
    int L = 4;
    Lattice lat = Lattice(H);
    Lattice newLat = Lattice(H_read);
    Hubbard hubby = Hubbard(4,lat);
    hubby.setLattice(newLat);
    double eigenTest = hubby.getGroundstates().at(0).eigenValue;
    double eigenExpected = -3.64898672;
    EXPECT_DOUBLE_EQ(eigenTest,eigenExpected);
    State x = hubby.getGroundstates().at(0);
    int sz = x.S_z;

    std::cout<<std::endl<<sz<<std::endl;
    Eigen::VectorXd v1;
    v1 = x.eigenVector;
    std::vector<double> v2;
    v2.resize(v1.size());
    Eigen::VectorXd::Map(&v2[0], v1.size()) = v1;
    OneRDM calc_test_rdm = OneRDM(v2,(N+sz)/2,(N-sz)/2,L);



    std::vector<double> testvec = {-0.14892,0.173082,3.69901e-16,2.94018e-16,-0.173082,0.189731,0.173082,8.50186e-16,-0.173082,-0.173082,0.379462,-0.173082,-1.7869e-16,-0.173082,0.14892,0.189731,-0.173082,1.87672e-16,-1.81443e-16,-0.173082,0.189731,0.14892,-0.173082,1.91932e-16,-0.173082,0.379462,-0.173082,-0.173082,7.76026e-16,0.173082,0.189731,-0.173082,-3.0344e-16,-2.85689e-16,0.173082,-0.14892};
    OneRDM calc_expect_rdm = OneRDM(testvec,2,2,4);
    arma::mat test_calc = calc_test_rdm.getRdmAlpha();
    arma::mat test_expect = calc_expect_rdm.getRdmAlpha();


    for(int i=0; i<test_calc.n_cols;i++){
        for(int j=0; j<test_calc.n_cols;j++){
            EXPECT_DOUBLE_EQ(test_calc(i,j),test_expect(i,j));
        }
    }

}