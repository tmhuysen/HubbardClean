//
// Created by Wulfix on 16/11/2017.
//

#ifndef HUBBARDCLEAN_HUBBARD_H
#define HUBBARDCLEAN_HUBBARD_H


#include "Lattice.h"
#include "AddressingMatrix.h"

class Hubbard
{
public:
    Hubbard(unsigned N, const Lattice& lattice);

private:
    unsigned N; // The number of electrons to place in the lattice
    const Lattice& lattice; // REFERENCE to a lattice. A reference is used to avoid copying.

    // Hubbard will serve as SpinSector factory.
    // this->spin_sectors is a std::vector of SpinSectors.
    //      They are ordered with increasing total projected spin
    //          spin_sectors[0]     is the spin sector with the lowest number of alpha electrons (0)
    //          spin_sectors[N]   is the spin sector with the highest number of alpha electrons (N)
    std::vector<SpinSector> spinSectors;
    // For each distribution of electrons we will need a different addressingsMatrix
    std::vector<AddressingMatrix> addressingMatrixs;

    void generateAddressingMatrix();
    void calculate();


    /**
     * nesting to prevent circular inclusion
     */
    class SpinSector {
    private:
        Hubbard& p; //parent
        int S_z;            // Total projected spin !!!TIMES 2!!!
        //      For total projected spin = 1/2, S_z = 1
        //      For total projected spin = -1, S_z = -2
        unsigned N_alpha;   // Number of alpha electrons
        unsigned N_beta;    // Number of beta electrons

        unsigned long N_bf;

        arma::mat hamiltonian;
        arma::vec eigenvalues;
        arma::mat eigenvectors;

        /**
        * nesting to prevent circular inclusion
        */
        class SingleSpinEvaluation{
        private:
            unsigned N_electrons;
            SpinSector& p; //parent
        public:
            SingleSpinEvaluation();
        };


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





};


#endif //HUBBARDCLEAN_HUBBARD_H
