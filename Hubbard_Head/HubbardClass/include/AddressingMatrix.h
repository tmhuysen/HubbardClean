//
// Created by Wulfix on 17/11/2017.
//

#ifndef HUBBARDMODEL_ADDRESSINGMATRIX_H
#define HUBBARDMODEL_ADDRESSINGMATRIX_H

#include "HubbardClass/include/bitset.h"

#include <Eigen/Dense>
class AddressingMatrix {
private:
    int n_sites;
    int n_electrons;
    Eigen::MatrixXi addressMatrix;

    void generateMatrix();

public:
    AddressingMatrix()= default;
    AddressingMatrix(unsigned int n_sites, unsigned int n_electrons);
    unsigned long fetchAddress(boost::dynamic_bitset<> bitVector);
    boost::dynamic_bitset<> generateBinaryVector(unsigned long address);

    void print();




};


#endif //HUBBARDMODEL_ADDRESSINGMATRIX_H
