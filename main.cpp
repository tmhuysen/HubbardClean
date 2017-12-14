#include <iostream>
#include <fstream>
#include <include/OneRDM.h>
#include <include/Hubbard.h>
#include <include/io.h>

void printVec(Eigen::VectorXd x){
    std::cout<<"[";

    std::cout<<x[0];
    for(int i = 1; i<x.innerSize();i++){
        std::cout<<",";
        std::cout<<x[i];

    }
    std::cout<< "]" <<std::endl;


}
void printVec2(std::vector<double> x){
    std::cout<<"[";

    std::cout<<x[0];
    for(int i = 1; i<x.size();i++){
        std::cout<<",";
        std::cout<<x[i];

    }
    std::cout<< "]" <<std::endl;


}
int main() {
    std::string filename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/randomized_fourring.txt";
    std::string fileout = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/ra.txt";
    std::string fileout2 = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/Hubbards_Qmark/HubbardClean/randomized_fourring_rdms.txt";
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
