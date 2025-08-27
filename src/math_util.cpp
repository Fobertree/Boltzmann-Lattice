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

void roll(ThreeD& m, int shift, int axis, int idx) {
    // block ops to save time

    if (idx > m.size()) {
        printf("ROLL::IDX OUT OF BOUNDS");
        return;
    }
    MatrixXd& submat = m[idx];

    if (axis == 0) {
        // by row
        int rows = static_cast<int>(submat.rows());
        if (rows == 0) return;
        int s = (shift % rows + rows) % rows;

        MatrixXd tmp = submat;
        submat.topRows(s) = tmp.bottomRows(s);
        submat.bottomRows(rows - s) = tmp.topRows(rows - s);
    }
    else if (axis == 1) {
        // column
        int cols = static_cast<int>(submat.cols());
        if (cols == 0) return;
        int s = (shift % cols + cols) % cols;

        MatrixXd tmp = submat;
        submat.leftCols(s) = tmp.rightCols(s);
        submat.rightCols(cols - s) = tmp.leftCols(cols - s);
    }
}


MatrixXd apply_boundary(ThreeD& F, MatrixXb& boundary, int sz) {
    // bndryF = F[cylinder,:]
    // F
    // TODO: This doesn't seem right
    constexpr static int col_sz = 9; // 9 point lattice
    sz = 0;

    for (int i = 0; i < boundary.rows(); i++) {
        for (int j = 0; j < boundary.cols(); j++) {
            if (boundary(i,j)) sz++;
        }
    }

    MatrixXd out(sz, col_sz);
    out.setZero();
    int out_idx{0};

    for (int i = 0; i < F[0].rows(); i++) {             // Ny
        for (int j = 0; j < F[0].cols(); j++) {         // Nx
            if (boundary(i,j)) {
                // ThreeD[NL][Ny]
//                out.row(out_idx++) = F[i].row(j);
                Eigen::VectorXd newRow(col_sz);
                newRow <<
                        F[0](i,j),
                        F[5](i,j),
                        F[6](i,j),
                        F[7](i,j),
                        F[8](i,j),
                        F[1](i,j),
                        F[2](i,j),
                        F[3](i,j),
                        F[4](i,j);

                out.row(out_idx++) = newRow;
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
    // if want to permute rows, do perm * mat
    m = m * perm;
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
    ThreeD out(NL, MatrixXd(Ny, Nx));
    for (int i = 0; i < m.size(); i++) {
        out[i].setZero();
        auto & mat = m[i];
        out[i] += mat * coeff[i];
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
    // F[:,cylinder] = bndryF
    for (MatrixXd& submat : F) {
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < F[i].rows(); j++) {
                if (bndry(i, j)) {
                    submat(i, j) = bndryF(i, j);
                }
            }
        }
    }
}

void apply_boundary_to_vel(MatrixXd &vel, MatrixXb &bndry) {
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

MatrixXd sum_axis_NL(ThreeD &m) {
    MatrixXd out(Ny, Nx);
    out.setConstant(0);

    for (const MatrixXd& submat : m) {
        out += submat;
    }
    return out;
}

MatrixXd get_velo(const ThreeD &F, const std::array<int, 9> &coeffs, const MatrixXd& rho) {
    MatrixXd vel(Ny, Nx); // (Ny, Nx)
    vel.setZero();

    for (int i = 0; i < NL; i++)
        vel += (F[i].array() * static_cast<double>(coeffs[i])).matrix();

    vel = vel.cwiseQuotient(rho);
    return vel;
}

template<typename T>
void vector_reorder(std::vector<T> &vec, const std::vector<int> & order) {
    for (int s = 1, d; s < order.size(); s++) {
        for (d = order[s]; d < s; d = order[d]) {
            if (d == s) {
                while (d = order[d], d != s)
                    std::swap(vec[s], vec[d]);
            }
        }
    }
}
