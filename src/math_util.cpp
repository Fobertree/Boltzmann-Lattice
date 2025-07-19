//
// Created by Alexander on 6/25/2025.
//
#include <math_util.h>

float squared_distance(float x1, float y1, float x2, float y2) {
    return (x2-x1) * (x2-x1) + (y2-y1) * (y2-y1);
}


void matrix_fill(MatrixXd& m, double num) {
    using ll = long long;
    ll r{m.rows()}, c{m.cols()};

    for (ll row = 0; row < r; row++) {
        for (ll col = 0; col < c; c++) {
            m(r, col) = num;
        }
    }
}

void roll(ThreeD& m, int shift, int axis) {
    if (axis == 0) {
        // by row
        for (MatrixXd& submat : m) {
            int rows = static_cast<int>(submat.rows());
            if (rows == 0) continue;
            shift = (shift & rows + rows) % rows;

            MatrixXd rolled = submat;
            for (int i = 0; i < rows; i++) {
                rolled.row((i+shift) % rows) = submat.row(i);
            }
            submat = rolled;
        }
    }
    else if (axis == 1) {
        // column
        for (MatrixXd& submat : m) {
            int cols = static_cast<int>(submat.cols());
            if (cols == 0) continue;
            shift = (shift % cols + cols) % cols;

            MatrixXd rolled = submat;
            for (int j = 0; j < cols; ++j) {
                rolled.col((j + shift) % cols) = submat.col(j);
            }
            submat = rolled;
        }
    }
}


MatrixXd apply_boundary(ThreeD& m, MatrixXb& boundary, int sz) {
    // bndryF = F[cylinder,:]
    constexpr static int col_sz = 9; // 9 point lattice
    printf("apply_boundary::bndry: (%td %td). SZ: %d\n", boundary.rows(), boundary.cols(), sz);

    sz = 0;

    for (int i = 0; i < boundary.rows(); i++) {
        for (int j = 0; j < boundary.cols(); j++) {
            if (boundary(i,j)) sz++;
        }
    }

    MatrixXd out(sz, col_sz);
    int out_idx{0};

    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[0].rows(); j++) {
            if (boundary(i,j)) {
                out.row(out_idx++) = m[i].row(j);
            }
        }
    }
    return out;
}

void three_d_swap(ThreeD& m, const std::array<int, 9>& order) {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(9);
    auto perm_indices = perm.indices();

    Eigen::VectorXd indices(9);
    for (int i = 0; i < 9; ++i) {
        perm_indices(i) = order[i];
    }

    for (auto & mat : m) {
        mat = mat * perm;
    }
}

void matrix_reorder(MatrixXd& m, const std::array<int, 9>& order) {
    // REORDER BY COLUMN
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(9);
    auto perm_indices = perm.indices();

//    Eigen::VectorXd indices(9);
    for (int i = 0; i < 9; ++i) {
        perm.indices()(i) = order[i];
    }

    printf("(%td, %td) x (%td, %td)\n", m.rows(), m.cols(), perm.rows(), perm.cols());
    // if want to permute rows, do perm * mat
    m = m * perm;
    printf("success!\n");
}

MatrixXd sum_axis_two(ThreeD &m) {
    // return matrixXd or ThreeD (vector of MatrixXd)
    // emulate np.sum axis = 2
    // np.sum flattens a dimension, so just return a MatrixXd
    if (m.empty()) return {};

    MatrixXd out(m.size(), m[0].rows());
    out.setZero();

    for (int i = 0; i < m.size(); i++) {
        MatrixXd& mat = m[i];
        // create new MatrixXd for sum axis 2
        for (int j = 0; j < mat.rows(); j++) {
            for (int k = 0; k < mat.cols(); k++) {
                out(i,j) += mat(j,k);
            }
        }
    }
    return out;
}

ThreeD element_prod(ThreeD &m, const std::array<int, 9> &coeff) {
    // element-wise numpy style multiplication
    // m * coeff
    ThreeD out(Ny, MatrixXd(Nx, NL));
    for (int i = 0; i < m.size(); i++) {
        out[i].setZero();
        auto & mat = m[i];
        for (int j = 0; j < mat.rows(); j++) {
            for (int k = 0; k < mat.cols(); k++) {
                out[i](j,k) += mat(j,k) * coeff[k];
            }
        }
    }
    return out;
}

MatrixXd element_div(MatrixXd &m1, MatrixXd &m2) {
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    MatrixXd out(m1.rows(), m1.cols());

    for (int i = 0; i < m1.rows(); i++) {
        for (int j = 0; j < m1.cols(); j++) {
            out(i, j) = m1(i, j) / m2(i, j);
        }
    }
    return out;
}

void apply_bndry_to_F(ThreeD &F, const MatrixXb &bndry, const MatrixXd &bndryF) {
    // F[cylinder,:] = bndryF
    for (int i = 0; i < F.size(); i++) {
        MatrixXd& m = F[i];
        for (int j = 0; j < F[i].rows(); j++) {
            if (bndry(i,j)) {
                m.row(j).setConstant(bndry(i,j));
            }
        }
    }
    printf("OK\n");
}

void apply_boundary_to_vel(MatrixXd &vel, MatrixXb &bndry) {
    printf("%td %td\n", vel.rows(), vel.cols());
    for (int i = 0; i < vel.rows(); i++) {
        for (int j = 0; j < vel.cols(); j++) {
            if (bndry(i,j)) {
                vel(i,j) = 0.;
            }
        }
    }
}

MatrixXd get_curl(MatrixXd &ux, MatrixXd &uy) {
    assert(ux.rows() == uy.rows() && ux.cols() == uy.cols());
    const static Eigen::Index rows = ux.rows(), cols = ux.cols();
    // dfydx = ux[2:,1:-1] - ux[0:-2, 1:-1]
    MatrixXd dfydx = ux(Eigen::seq(2,rows - 1), Eigen::seq(1,cols-2)) -
            ux(Eigen::seq(0, rows-3), Eigen::seq(1, cols-2));
    // dfxdy = uy[1:-1, 2:] - uy[1:-1, 0:-2]
    MatrixXd dfxdy = uy(Eigen::seq(1, rows-2), Eigen::seq(2, cols-1)) -
            uy(Eigen::seq(1, rows-2), Eigen::seq(0, cols-3));
    return dfydx - dfxdy;
}
