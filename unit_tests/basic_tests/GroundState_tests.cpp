//
// Created by Wulfix on 09/12/2017.
//

#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <include/OneRDM.h>
#include <include/Hubbard.h>
#include <include/io.h>

TEST(io_test_2, test){
    arma::arma_version ver;
    std::cout << "ARMA version: "<< ver.as_string() << std::endl;

    std::string filename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/unit_tests/basic_tests/input_data/test_input_23.txt";
    std::string fileout = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/unit_tests/basic_tests/input_data/test_rdm_23.txt";
    std::string fileout2 = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/unit_tests/basic_tests/input_data/test_ene_23.txt";
    std::ifstream is (filename);
    int counter = 0;
    std::ofstream outfile (fileout);
    std::ofstream outfile2 (fileout2);
    outfile<<std::setprecision(12);
    outfile2<<std::setprecision(12);
    arma::uword dim = 4;
    unsigned long elec = 4;
    unsigned long dimalt = 4;
    Lattice def = Lattice(dimalt);
    Hubbard instant = Hubbard(elec,def);
    char delimiter_char = ',';
    if (is.is_open()) {
        std::string line;
        while(std::getline(is, line)){
            try {
                counter++;

                arma::mat M = arma::zeros(dim, dim);

                // Read in the upper triangular part of the hopping matrix.
                // It is located on one line, so we can read in that line, and then break it down.
                std::stringstream linestream(
                        line); // convert the read line into a stringstream to perform std::getline() on it


                // The actual loop to read the upper triangular part of the hopping matrix
                for (arma::uword i = 0; i < dim; i++) {
                    for (arma::uword j = i; j < dim; j++) {
                        std::string value_as_string;
                        std::getline(linestream, value_as_string, delimiter_char);
                        M(i, j) = std::stod(value_as_string);
                    }
                }
                M = arma::symmatu(M);
                Lattice thisIn = Lattice(M);
                instant.setLattice(thisIn);
                std::vector<State> ground = instant.getGroundstates();
                outfile << ground.at(0).eigenValue << std::endl;
                outfile2 << ground.at(0).eigenValue << std::endl;
                for (State state : ground) {
                    int sz = state.S_z;
                    Eigen::VectorXd v1;
                    v1 = state.eigenVector;
                    std::vector<double> v2;
                    v2.resize(v1.size());
                    Eigen::VectorXd::Map(&v2[0], v1.size()) = v1;
                    OneRDM onesy = OneRDM(v2, (elec + sz) / 2, (elec - sz) / 2, dimalt);
                    outfile2<<std::endl<<sz<<" en "<<state.eigenValue<<std::endl;
                    writeOutTriMat(outfile2, onesy.getRdmAlpha());
                    writeOutTriMat(outfile2, onesy.getRdmBeta());



                }
                outfile2<<std::endl;
            }catch(const std::invalid_argument& ia){
                std::cerr << "Invalid argument: " << ia.what() << '\n';
                std::cout<<std::endl<<counter<<" WTF";

            }








        }
        is.close();
        outfile.close();
        outfile2.close();
    } else std::cout << "Unable to open file";





}