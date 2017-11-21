//
// Created by Wulfix on 18/11/2017.
//

#include "HubbardClass/include/Hubbard.h"


Hubbard::SpinSector::SpinSector(Hubbard &p, int S_z): p(p)
{

    this->S_z = S_z;

    // Convert N and S_z to N_alpha and N_beta
    this->N_alpha = (this->p.N + this->S_z) / 2;
    this->N_beta = (this->p.N - this->S_z) / 2;

    // The number of basis functions in this sector is (L choose N_alpha) x (L choose N_beta).
    // We will store both dimensions of basis vectors separated by spin and combined.
    auto n_bf_alpha_ = boost::math::binomial_coefficient<double>(this->p.lattice.getLength(), this->N_alpha);
    if (n_bf_alpha_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into unsigned long. - ALPHA");
    }
    this->n_bf_alpha = static_cast<unsigned long>(n_bf_alpha_);

    auto n_bf_beta_ = boost::math::binomial_coefficient<double>(this->p.lattice.getLength(), this->N_beta);
    if (n_bf_beta_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into unsigned long. - BETA");
    }
    this->n_bf_beta = static_cast<unsigned long>(n_bf_beta_);


    auto nbf_ = n_bf_alpha_*n_bf_beta_;
    if (nbf_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the sector is too high to be cast into unsigned long.");
    }
    this->n_bf = static_cast<unsigned long>(nbf_);

    // Upon initialization, allocate memory for the Hamiltonian matrix. Make sure this initialized a zero matrix.
    this->hamiltonian = Eigen::MatrixXd::Zero(this->n_bf, this->n_bf);

    calculateSector(0.0,1.0);


}



void Hubbard::SpinSector::calculateSector(double start, double end) {
    AddressingMatrix ad_a = AddressingMatrix(p.lattice.getLength(),N_alpha);
    AddressingMatrix ad_b = AddressingMatrix(p.lattice.getLength(),N_beta);
    //AddressingMatrix ad_b = p.addressing_list.at(N_beta);
    //Iterate over the bitvectors
    //These are separated (not sure if i'll keep this way)


    boost::dynamic_bitset<> beta_set = ad_b.generateBinaryVector(start * n_bf_beta);

    for (size_t i = 0; i < n_bf_beta * end; i++) {
        boost::dynamic_bitset<> alpha_set = ad_a.generateBinaryVector(start * n_bf_alpha);
        for (size_t j = 0; j < n_bf_alpha * end; j++) {


            //Itterate over hopping elements

            for (size_t x = 0; x < p.lattice.getLength(); x++) {
                for (size_t y = x; y < p.lattice.getLength(); y++) {

                    double element = p.lattice.getElement(x, y);

                    //If only evaluate non zero elements
                    if (std::abs(element) > 1.0e-06) {
                        boost::dynamic_bitset<> beta_target = boost::dynamic_bitset<>(beta_set);
                        boost::dynamic_bitset<> alpha_target = boost::dynamic_bitset<>(alpha_set);
                        if (x == y) {

                            if (annihilation(beta_target, x) && annihilation(alpha_target, x) &&
                                creation(alpha_target, x) && creation(beta_target, x)) {
                                addToHamiltonian(element, i * n_bf_alpha + j, i * n_bf_alpha + j);
                            }
                        } else {
                            if (annihilation(alpha_target, x) && creation(alpha_target, y)) {
                                double phase_factor = phaseCheck(alpha_set, alpha_target);
                                size_t address = ad_a.fetchAddress(alpha_target);
                                addToHamiltonian(phase_factor * element, i * n_bf_alpha + j, i * n_bf_alpha + address);
                            }
                            if (annihilation(beta_target, x) && creation(beta_target, y)) {
                                double phase_factor = phaseCheck(beta_set, beta_target);
                                size_t address = ad_b.fetchAddress(beta_target);
                                addToHamiltonian(phase_factor * element, i * n_bf_alpha + j, address * n_bf_alpha + j);
                            }

                        }
                    }
                }
            }
            alpha_set = next_bitset_permutation(alpha_set);

        }
        beta_set = next_bitset_permutation(beta_set);
    }


    symmetryFill();


}

void Hubbard::SpinSector::addToHamiltonian(double value, size_t index1, size_t index2) {
    hamiltonian(index1, index2) = value;

}

void Hubbard::SpinSector::symmetryFill() {
    Symmatu(hamiltonian);

}

const Eigen::MatrixXd &Hubbard::SpinSector::getHamiltonian() const {
    return hamiltonian;
}

int Hubbard::SpinSector::getS_z() const {
    return S_z;
}

void Hubbard::SpinSector::setEigenvalues(const Eigen::VectorXd &eigenvalues) {
    SpinSector::eigenvalues = eigenvalues;
}

void Hubbard::SpinSector::setEigenvectors(const Eigen::MatrixXd &eigenvectors) {
    SpinSector::eigenvectors = eigenvectors;
}

void Hubbard::SpinSector::print_hamiltonian() const {
    std::cout<<hamiltonian;

}




