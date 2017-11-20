//
// Created by Wulfix on 19/11/2017.
//

#include <utility>

#include "HubbardClass/include/OneRDM.h"


OneRDM::OneRDM(std::vector<double> coefs, unsigned nup, unsigned ndown, unsigned sites) {

    this -> nup = nup;
    this -> ndown = ndown;
    this -> sites = sites;
    this -> coefs = std::move(coefs);
    this -> ad_a = AddressingMatrix(sites,nup);
    this -> ad_b = AddressingMatrix(sites,ndown);

    auto n_bf_alpha_ = boost::math::binomial_coefficient<double>(this->sites, nup);
    if (n_bf_alpha_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into unsigned long.");
    }
    this->n_bf_alpha = static_cast<unsigned long>(n_bf_alpha_);

    auto n_bf_beta_ = boost::math::binomial_coefficient<double>(sites, ndown);
    if (n_bf_beta_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the separated sector is too high to be cast into unsigned long.");
    }
    this->n_bf_beta = static_cast<unsigned long>(n_bf_beta_);
    oneRDMDim = (sites*(sites+1))/2;
    rdmAlpha[oneRDMDim] = {};
    rdmBeta[oneRDMDim] = {};
    rdmCalc();


}

void OneRDM::rdmCalc(){
    for (size_t i = 0; i<sites; i ++){
        for(int j = 0; j<i+1; j++){
            rdmAlpha[(((i)*(i+1))/2 +j)] = oneRDMAlpha(i,j);
            rdmBeta[(((i)*(i+1))/2 +j)] = oneRDMBeta(i,j);


        }
    }

}

double OneRDM::oneRDMAlpha(unsigned i, unsigned j) {
    double coefsum = 0;
    boost::dynamic_bitset<> alpha_set = ad_a.generateBinaryVector(0);
    for (int l = 0; l<n_bf_alpha ; l++ ){
        boost::dynamic_bitset<> alpha_target = alpha_set;
        double element_set = coefs.at(l);
        if (annihilation(alpha_target, i) && creation(alpha_target, j)) {
            double phase_factor = phaseCheck(alpha_set, alpha_target);
            size_t address = ad_a.fetchAddress(alpha_target);
            for(size_t k = 0; k<n_bf_beta ; k++){
                coefsum += coefs.at(address + k*n_bf_beta)*element_set*phase_factor;
            }
        }
        alpha_set = next_bitset_permutation(alpha_set);

    }
    return coefsum;
}




double OneRDM::oneRDMBeta(unsigned i, unsigned j){
    double coefsum = 0;
    boost::dynamic_bitset<> beta_set = ad_b.generateBinaryVector(0);
    for (int l = 0; l<n_bf_alpha ; l++ ){
        boost::dynamic_bitset<> beta_target = beta_set;
        double element_set = coefs.at(l);
        if (annihilation(beta_target, i) && creation(beta_target, j)) {
            double phase_factor = phaseCheck(beta_set, beta_target);
            size_t address = ad_b.fetchAddress(beta_target);
            for(size_t k = 0; k<n_bf_alpha ; k++){
                coefsum += coefs.at(address + k)*element_set*phase_factor;
            }
        }
        beta_set = next_bitset_permutation(beta_set);

    }
    return coefsum;
}






