//
// Created by Alexander on 6/14/2025.
//

#ifndef BOLTZMANNLATTICE_UTILS_H
#define BOLTZMANNLATTICE_UTILS_H

#include <array>
#include <limits>
#include <Eigen/Dense>
#include <chrono>

template <typename T>
inline void IFUCKEDUP(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)  {
    printf("PRINT MAT: (%td %td)\n",mat.rows(), mat.cols());
}

constexpr int WIDTH{800};
constexpr int HEIGHT{800};
constexpr double FPS_INTERVAL{5.0};
constexpr const char* FRAG_SHADER_PATH = "/Users/alexanderliu/CLionProjects/Boltzmann-Lattice/src/Shaders/clm.frag";
constexpr const char* VERT_SHADER_PATH = "/Users/alexanderliu/CLionProjects/Boltzmann-Lattice/src/Shaders/main.vert";
constexpr bool CURL_FLAG = false;

#endif //BOLTZMANNLATTICE_UTILS_H
