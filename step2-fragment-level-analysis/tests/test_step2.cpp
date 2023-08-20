//
// Created by Goran Vinterhalter on 16/08/2023.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <string>

#include <catch2/catch_all.hpp>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "npy/npy.hpp"
#include "modalib/sample.h"
#include "modalib/features.h"
#include "modalib/kernel.h"


using namespace std;


void load_nuc_map(unordered_map<uint8_t, vector<long>> nuc_pos, const string &dir_path) {

    for(int chr_id = 1; chr_id < 27; chr_id++){
        auto file_path = dir_path + std::to_string(chr_id) + ".npy";

        // TODO: check if the file exists

        // output from reading a file
        std::vector<unsigned long> shape;
        bool fortran_order;
        std::vector<double> raw_data;

        // load raw_data
        npy::LoadArrayFromNumpy(file_path, shape, fortran_order, raw_data);

        // convert vector of double to vector of long (TODO why is it originally double)
        std::vector<long> data(raw_data.size());
        for(long i=0; i<raw_data.size(); i++){
            data[i] = static_cast<long>(raw_data[i]);
        }

        nuc_pos[chr_id] = std::move(data);
    }

//    // Print the contents of the map
//    for (auto & [key, value] : nuc_pos) {
//        fmt::print("key: {}, : {}\n", key, value);
//    }
}

void run() {

    // create a vector containing one sample
    string sample_path = STEP_1_OUT_DIR "GC110680";
    vector<modalib::_Sample> samples;
    samples.emplace_back(sample_path);

    // Load reference nucleosome map
    unordered_map<uint8_t, vector<long>> nuc_pos;
    load_nuc_map(nuc_pos, NUCLEOSOMES_DIR);

    // extract features
//    auto f_matrices =  modalib::extract_all_features(samples, nuc_pos);

    // extract features
//    auto k_matrices =  modalib::compute_kernel_matrices(samples, nuc_pos);


}


TEST_CASE("Step2 Test") {

    SECTION("Run one sample") {
        run();
    }
}
