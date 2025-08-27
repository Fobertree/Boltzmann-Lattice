//
// Created by Alexander on 6/14/2025.
//

#ifndef BOLTZMANNLATTICE_UTILS_H
#define BOLTZMANNLATTICE_UTILS_H

#include <array>
#include <limits>
#include <Eigen/Dense>
#include <chrono>

inline void IFUCKEDUP(const Eigen::MatrixXd& mat)  {
    printf("PRINT MAT: (%td %td)\n",mat.rows(), mat.cols());
}

constexpr int WIDTH{800};
constexpr int HEIGHT{800};
constexpr const char* FRAG_SHADER_PATH = "/Users/alexanderliu/CLionProjects/Boltzmann-Lattice/src/Shaders/clm.frag";
constexpr const char* VERT_SHADER_PATH = "/Users/alexanderliu/CLionProjects/Boltzmann-Lattice/src/Shaders/tmp.vert";
constexpr bool CURL_FLAG = false;

template <typename T, size_t N>
constexpr T sum(const std::array<T, N>& arr) {
    T result = 0;
    for (const T& element : arr) {
        result += element;
    }
    return result;
}

template <typename T>
constexpr bool approx_equal(const T n1, const T n2) {
    return (std::abs(n1 - n2) < std::numeric_limits<T>(n1 - n2));
}

#endif //BOLTZMANNLATTICE_UTILS_H
