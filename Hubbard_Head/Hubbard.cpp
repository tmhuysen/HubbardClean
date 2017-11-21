//
// Created by Wulfix on 16/11/2017.
//

#include "Hubbard.h"

Hubbard::Hubbard(unsigned N, const Lattice &lattice): lattice(lattice) {
    if (N > (2 * this->lattice.getLength())) {
        throw std::invalid_argument("The lattice is not compatible with the given number of electrons.");
    } else {
        this->N = N;
    }
    generateAddressingMatrix();

    calculate();
    solvedHubbard.push_back(HubbardSolver(*this));
}

void Hubbard::generateAddressingMatrix() {
    for(size_t electrons = 0; electrons<N+1;electrons++){
        addressing_list.push_back(AddressingMatrix(lattice.getLength(),electrons));
    }

}

void Hubbard::calculate() {
    for(unsigned i= 0; i<N+1;i++){
        spinSectors.push_back(SpinSector(*this,-static_cast<int>(N) + 2 * i));
    }


}

const std::vector<Hubbard::SpinSector> &Hubbard::getSpinSectors() const {
    return spinSectors;
}

void Hubbard::print(){
    for(auto x: spinSectors){
        std::cout<<std::endl;
        x.print_hamiltonian();
        std::cout<<std::endl;
        std::cout<<x.getS_z();
    }





}

std::vector<State> Hubbard::getGroundstates(){
        return solvedHubbard.at(0).getGroundstates();
};
