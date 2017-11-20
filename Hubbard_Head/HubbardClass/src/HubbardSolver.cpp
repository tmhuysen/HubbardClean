//
// Created by Wulfix on 19/11/2017.
//

#include "HubbardClass/include/Hubbard.h"

Hubbard::HubbardSolver::HubbardSolver(Hubbard &hubbard) : hubbard(hubbard){
    // Create groundstate with higehest energy possible (avoid using magic numbers).
    groundstates = { State {std::numeric_limits<double>::max(),arma::vec()} };
    // Diagonalize the Hamiltonian for every spin sector, and assign the eigenvalues and eigenvectors to the corresponding SpinSector instance
    for (auto& spin_sector : this->hubbard.spinSectors) {      // use reference, otherwise we would make a copy...
        auto H = spin_sector.getHamiltonian();

        arma::vec eigenvalues;
        arma::mat eigenvectors;

        arma::eig_sym(eigenvalues, eigenvectors, H);
        //We store all eigen solutions
        spin_sector.setEigenvalues(eigenvalues);
        spin_sector.setEigenvectors(eigenvectors);
        //We extract only the groundstate
        for (int i = 0; i<eigenvalues.size(); i++) {
            groundStates(State {eigenvalues[i], eigenvectors.col(static_cast<const arma::uword>(i)), spin_sector.getS_z()});
        }
    }
}
/** Adds a state to the groundstate list
 *     if the state is lower than the current energy of all assumed "groundstates"
 *     the list is cleared and the new lowest state is added.
 */
void Hubbard::HubbardSolver::groundStates(State key) {

    if(areSame(key,this->groundstates.at(0))){
        groundstates.push_back(key);
    }
    else{
        if(compareState(key,groundstates.at(0))){
            groundstates = std::vector<State> {key};
        }
    }


}

Hubbard::HubbardSolver::HubbardSolver() : hubbard(Hubbard()){

}


bool compareState(const State &o1, const State &o2){
    return o1.eigenValue < o2.eigenValue;
}
bool areSame(const State &o1, const State &o2) {

    double ELIPSON = (o1.eigenValue > o2.eigenValue) ?  o2.eigenValue/100000 : o1.eigenValue/100000;

    return fabs(o1.eigenValue - o2.eigenValue) < fabs(ELIPSON);
}