//
// Created by Wulfix on 20/11/2017.
//

#include "gtest/gtest.h"
#include "HubbardClass/include/OneRDM.h"

TEST(OneRDM_test, test_output){
    std::vector<double> testDoubles = {0.25,-0.25,0,0,0.25,-0.25,-0.25,0.25,0,0,-0.25,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0.25,-0.25,0,0,0.25,-0.25,-0.25,0.25,0,0,-0.25,0.25};
    OneRDM OneRDM_test_test_output_Test = OneRDM(testDoubles,2,2,4);
    OneRDM_test_test_output_Test.print();

    std::vector<double> testDoubles2 = {-0.25,0.25,-0.25,0.25,0.25,-0.25,0.25,-0.25,-0.25,0.25,-0.25,0.25,0.25,-0.25,0.25,-0.25};
    OneRDM OneRDM_test_test_output_Test2 = OneRDM(testDoubles2,1,3,4);
    OneRDM_test_test_output_Test2.print();

    std::vector<double> testDoubles3 = {-2,2,-2,-2,2,-2,-2,2,-2};
    OneRDM OneRDM_test_test_output_Test3 = OneRDM(testDoubles3,1,2,3);
    OneRDM_test_test_output_Test3.print();
}


