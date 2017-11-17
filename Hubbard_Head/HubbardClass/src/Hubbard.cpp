//
// Created by Wulfix on 16/11/2017.
//

#include "HubbardClass/include/Hubbard.h"

Hubbard::Hubbard(unsigned N, const Lattice &lattice): lattice(lattice) {
    if (N > (2 * this->lattice.getLength())) {
        throw std::invalid_argument("The lattice is not compatible with the given number of electrons.");
    } else {
        this->N = N;
    }
    generateAddressingMatrix();
    calculate();
}

void Hubbard::generateAddressingMatrix() {
    for(unsigned electrons = 0; electrons<N;electrons++){
        addressingMatrixs.push_back(AddressingMatrix(lattice.getLength(),electrons));
    }

}

void Hubbard::calculate() {

}
