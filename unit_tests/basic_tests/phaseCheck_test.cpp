//
// Created by Wulfix on 20/11/2017.
//

#include "gtest/gtest.h"
#include "include/bitset.h"

TEST(phaseCheck_test, test_output){

    boost::dynamic_bitset<> x(3);
    x[0] = true;
    boost::dynamic_bitset<> y(3);
    y[0] = true;
    double test = phaseCheck(x,y);
    std::cout<<std::endl<<test<<std::endl;


}