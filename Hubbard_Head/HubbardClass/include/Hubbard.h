//
// Created by Wulfix on 16/11/2017.
//

#ifndef HUBBARDCLEAN_HUBBARD_H
#define HUBBARDCLEAN_HUBBARD_H


#include "Lattice.h"
#include "AddressingMatrix.h"
#include <boost/math/special_functions.hpp>

struct State {
    double eigenValue;          // The energy of the solution, a.k.a. the eigenvalue
    arma::vec eigenVector; // The coefficients of the solution with respect to the given basis, a.k.a. the eigenvector corresponding to the eigenvalue
    int S_z;            // Total projected spin !!!TIMES 2!!!
    //      For total projected spin = 1/2, S_z = 1
    //      For total projected spin = -1, S_z = -2
};



class Hubbard
{
private:

    /**
     * nesting to prevent circular inclusion
     */
    class SpinSector {
    private:
        Hubbard& p; //parent
        int S_z;
        arma::vec eigenvalues;
        arma::mat eigenvectors;
    public:
        void setEigenvectors(const arma::mat &eigenvectors);

    public:
        void setEigenvalues(const arma::vec &eigenvalues);

    public:
        int getS_z() const;

    private:
        // Total projected spin !!!TIMES 2!!!
        //      For total projected spin = 1/2, S_z = 1
        //      For total projected spin = -1, S_z = -2
        unsigned N_alpha;   // Number of alpha electrons
        unsigned N_beta;    // Number of beta electrons

        unsigned long n_bf; // number of total basis functions
        unsigned long n_bf_alpha;
        unsigned long n_bf_beta;

        arma::mat hamiltonian;
    public:
        const arma::mat &getHamiltonian() const;

    private:

        void calculateSector(double start, double end);

        void addToHamiltonian(double value, size_t index1, size_t index2);
        void symmetryFill();
        //void calculateOneSpin(unsigned long bf_m_el, unsigned n_m_el,unsigned long bf_c_el, unsigned n_c_el, unsigned long start, unsigned long end);



    public:

        /** Default constructor
         *
         */
        SpinSector();

        /** Create a spin sector with
         *      S_z:    integer representation of total projected spin (times 2)
         *                  total projected spin = 1/2 -> S_z = 1
         *                  total projected spin = -1 -> S_z = -2
         */
        SpinSector(Hubbard& p, int S_z);




        /** Helper function to print all the basis vectors belonging to the spin sector
         */
        void print_basis() const;

        /** Helper function to print the Hamiltonian belonging to the spin sector
        */
        void print_hamiltonian() const;

    };

    /**
     * Solver
     */
    class HubbardSolver {

    private:
        Hubbard &hubbard;   // REFERENCE to a Hubbard instance. In this case, a reference IS needed, as we're going to make changes to the Hubbard instance (modifying the eigenvalues and eigenvectors of all its SpinSectors).

        std::vector<State> groundstates;      // The ground states  of the Hubbard system (i.e. the state with the lowest energy).


        /** Adds a state to the groundstate list
         * if the state is lower than the current energy of all assumed "groundstates"
         * the list is cleared and the new lowest state is added.
        */
        void groundStates(State key);
    public:
        /** Constructor based on a given Hubbard instance
         */
        HubbardSolver(Hubbard& hubbard);
        //HubbardSolver();


    };
public:
    Hubbard(unsigned N, const Lattice& lattice);
     //Hubbard & Hubbard():lattice(Lattice()){};
    const std::vector<SpinSector> &getSpinSectors() const;

private:
    unsigned N; // The number of electrons to place in the lattice
    const Lattice& lattice; // REFERENCE to a lattice. A reference is used to avoid copying.

    // Hubbard will serve as SpinSector factory.
    // this->spin_sectors is a std::vector of SpinSectors.
    //      They are ordered with increasing total projected spin
    //          spin_sectors[0]     is the spin sector with the lowest number of alpha electrons (0)
    //          spin_sectors[N]   is the spin sector with the highest number of alpha electrons (N)
    std::vector<SpinSector> spinSectors;
    std::vector<HubbardSolver> solver;

    //HubbardSolver solvedHubbard;

    // For each distribution of electrons we will need a different addressingsMatrix these may be same for several spin sectors
    // So we storing them allows us to only generate them once.
    std::vector<AddressingMatrix> addressing_list;

    void generateAddressingMatrix();
    void calculate();




};

//Compare tools for the State struct.
bool compareState(const State &o1, const State &o2);
bool areSame(const State &o1, const State &o2);


#endif //HUBBARDCLEAN_HUBBARD_H
