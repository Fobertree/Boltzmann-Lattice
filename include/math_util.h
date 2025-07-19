//
// Created by Alexander on 6/14/2025.
//

#ifndef BOLTZMANNLATTICE_MATH_UTIL_H
#define BOLTZMANNLATTICE_MATH_UTIL_H

#include <vector>
#include <format>
#include <type_traits>
#include <array>
#include "utils.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;
using ThreeD = std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd>>;

// MACROS AND CONSTS
// can change
constexpr int Nx = 400;
constexpr int Ny = 100;
constexpr double tau = 0.53;
constexpr int Nt = 30000;
constexpr int NL = 9;
constexpr double RIGHT_VEL = 2.3;
constexpr double CYLINDER_RADIUS_SQUARED = 13 * 13;

// typedefs
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

// likely some debate over std::array vs c-style array
// DO NOT CHANGE
constexpr std::array cxs = {0, 0, 1, 1, 1, 0,-1,-1,-1};
constexpr std::array cys = {0, 1, 1, 0,-1,-1,-1, 0, 1};
constexpr std::array<double, 10> weights = {4./9., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36.};

//static_assert(approx_equal(sum(weights), 1.), "Sum of weights must be 1")

// function declarations
float squared_distance(float x1, float y1, float x2, float y2);

void roll(ThreeD& m, int shift, int axis);

// bndryF = F[cylinder,:]
MatrixXd apply_boundary(ThreeD& m, MatrixXb& boundary, int sz);

void matrix_reorder(MatrixXd& m, const std::array<int, 9>& order);

MatrixXd sum_axis_two(ThreeD& m);

// name is a bit of a misnomer, just built for this program
ThreeD element_prod(ThreeD& m, const std::array<int, 9>& coeff);

MatrixXd element_div(MatrixXd& m1, MatrixXd& m2);

void apply_bndry_to_F(ThreeD& F, const MatrixXb& bndry, const MatrixXd& bndryF);

void apply_boundary_to_vel(MatrixXd& vel, MatrixXb& bndry);

MatrixXd get_curl(MatrixXd& ux, MatrixXd& uy);

// end necessary functions


// above functions should just be declarative to avoid stack frames

#endif //BOLTZMANNLATTICE_MATH_UTIL_H
