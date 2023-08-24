//
// Created by Goran Vinterhalter on 16/08/2023.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <string>

#include <catch2/catch_all.hpp>
#include <cnpy.h>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "npy/npy.hpp"

using namespace std;

TEST_CASE("Test usage of libnpy and cnpy"){

    const std::string NPY_PATH = NUCLEOSOMES_DIR "14.npy";

    SECTION("Test libnpy load") {
        std::vector<unsigned long> shape {};
        bool fortran_order;
        std::vector<double> data;

        npy::LoadArrayFromNumpy(NPY_PATH, shape, fortran_order, data);
//        fmt::print("data: {}", data);
    }

    SECTION("Test cnpy load") {

        cnpy::NpyArray arr = cnpy::npy_load(NPY_PATH);

        unsigned long n = arr.shape[0]; // length
        auto * loaded_data = arr.data<double>();

        // convert to vector<double>
        std::vector<double> data(n);
        for (long i = 0; i<n; i++){
            data[i] = loaded_data[i];
        }
//        fmt::print("data: {}", data);
    }

}

