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
constexpr int NL = 9;                 // lattice points/directions
constexpr double RIGHT_VEL = 2.3;
constexpr double CYLINDER_RADIUS_SQUARED = 13 * 13;

// typedefs
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

constexpr std::array cxs = {0, 0, 1, 1, 1, 0,-1,-1,-1};
constexpr std::array cys = {0, 1, 1, 0,-1,-1,-1, 0, 1};
constexpr std::array<double, 10> weights = {4./9., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36.};

// function declarations
float squared_distance(float x1, float y1, float x2, float y2);

void roll(ThreeD& m, int shift, int axis, int idx);

// bndryF = F[cylinder,:]
MatrixXd apply_boundary(ThreeD& m, MatrixXb& boundary, int sz);

MatrixXd sum_axis_NL(ThreeD& m);

void apply_bndry_to_F(ThreeD& F, const MatrixXb& bndry, const MatrixXd& bndryF);

void apply_boundary_to_vel(MatrixXd& vel, MatrixXb& bndry);

MatrixXd get_curl(MatrixXd& ux, MatrixXd& uy);

MatrixXd get_velo(const ThreeD& F, const std::array<int, 9>& coeffs, const MatrixXd& rho);

#endif //BOLTZMANNLATTICE_MATH_UTIL_H
