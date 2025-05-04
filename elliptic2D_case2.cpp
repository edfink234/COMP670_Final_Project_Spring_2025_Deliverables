#include "mole.h"
#include <iostream>
#include <armadillo>
//#include <stdexcept>
#include <tuple>

using namespace arma;
constexpr double pi = 3.14159;

// Compute the one-dimensional first-derivative operator on a uniform nodal grid.
// k       : even order of accuracy (must be divisible by 2)
// m_nodes : number of grid nodes
// dx      : uniform spacing
sp_mat nodal(int k, int m, double dx) {
    m = m - 1;
    int n_rows = m + 1;
    int n_cols = n_rows;

    sp_mat N(n_rows, n_cols);

    int len = k + 1;
    vec neighbors(len);
    neighbors(0) = -k / 2.0;

    for (int i = 1; i < len; ++i) {
        neighbors(i) = neighbors(i-1) + 1.0;
    }

    // Create k x k Vandermonde matrix
    mat A = fliplr(pow(repmat(neighbors.t(), len, 1), repmat(regspace(0, len-1), 1, len)));

    // First-order derivative
    vec b = zeros<vec>(len);
    b(len-2) = 1.0;

    // Solve for coefficients
    vec coeffs = solve(A, b);

    int j = 0;
    for (int i = k/2; i < n_rows - k/2; ++i) {
        for (int l = 0; l < len; ++l) {
            N(i, j+l) = coeffs(l);
        }
        j++;
    }

    int p = k/2;
    int q = k+1;

    sp_mat A_boundary(p, q);
    for (int i = 0; i < p; ++i) {
        vec neighbors_boundary(q);
        neighbors_boundary(0) = 1 - (i+1);
        neighbors_boundary(1) = neighbors_boundary(0) + 1;

        for (int j = 2; j < q; ++j) {
            neighbors_boundary(j) = neighbors_boundary(j-1) + 1;
        }

        mat V = fliplr(pow(repmat(neighbors_boundary.t(), q, 1), repmat(regspace(0, q-1), 1, q)));
        vec b_boundary = zeros<vec>(q);
        b_boundary(q-2) = 1.0;

        vec coeffs_boundary = solve(V, b_boundary);
        for (int l = 0; l < q; ++l) {
            A_boundary(i, l) = coeffs_boundary(l);
        }
    }

    // Insert A_boundary into N (upper-left)
    for (int i = 0; i < p; ++i) {
        for (int l = 0; l < q; ++l) {
            N(i, l) = A_boundary(i, l);
        }
    }

    // Permutation matrices
    sp_mat Pp = fliplr(speye(p, p));
    sp_mat Pq = fliplr(speye(q, q));

    sp_mat A_bottom = -Pp * A_boundary * Pq;

    // Insert A_bottom into N (lower-right)
    for (int i = 0; i < p; ++i) {
        for (int l = 0; l < q; ++l) {
            N(n_rows - p + i, n_cols - q + l) = A_bottom(i, l);
        }
    }

    // Scale N
    N *= (1.0 / dx);

    return N;
}

/**
 * Assemble the 2D first-derivative operator on a uniform nodal grid.
 *
 * k   : even order of accuracy
 * m   : number of nodes along x-axis
 * dx  : grid spacing along x-axis
 * n   : number of nodes along y-axis
 * dy  : grid spacing along y-axis
 */
sp_mat nodal2D(int k, int m, double dx, int n, double dy) {
    // Build the 1D operators along each axis
    sp_mat Nx = nodal(k, m, dx);  // size m×m
    sp_mat Ny = nodal(k, n, dy);  // size n×n

    // Identities of matching sizes
    sp_mat Im = speye<sp_mat>(m, m);
    sp_mat In = speye<sp_mat>(n, n);

    // Kronecker products:
    //   Along x-derivative directions: repeat Nx down each block of In
    sp_mat Kx = kron(In, Nx);        // size (n*m)×(n*m)
    //   Along y-derivative directions: repeat Ny across each block of Im
    sp_mat Ky = kron(Ny, Im);        // size (n*m)×(n*m)

    // Stack them vertically so that the first n*m rows approximate ∂/∂x
    // and the next n*m rows approximate ∂/∂y
    sp_mat N2 = join_cols(Kx, Ky);   // size (2*n*m)×(n*m)

    return N2;
}

/**
 * Compute the Jacobian determinant and metric derivatives on a 2D nodal grid.
 *
 * @param k   Even order of accuracy
 * @param X   Physical x–coordinates, size n×m
 * @param Y   Physical y–coordinates, size n×m
 * @return    A tuple (J, Xe, Xn, Ye, Yn), each an vec of length m*n:
 *            - J:  Xe % Yn  -  Xn % Ye    (elementwise)
 *            - Xe: ∂x/∂e metric, from first m*n entries
 *            - Xn: ∂x/∂n metric, from last  m*n entries
 *            - Ye: ∂y/∂e metric, from first m*n entries
 *            - Yn: ∂y/∂n metric, from last  m*n entries
 *
 * Requires:
 *   sp_mat nodal2D(int k, int m, double dx, int n, double dy);
 */
void jacobian2D(int k, const mat& X_in, const mat& Y_in, vec& J, vec& Xe, vec& Xn, vec& Ye, vec& Yn) {
    int n = X_in.n_rows;
    int m = X_in.n_cols;

    // Reshape X and Y (flattened column-major order like MATLAB's reshape(X', [], 1))
    vec X = vectorise(X_in.t());
    vec Y = vectorise(Y_in.t());

    // Build the nodal derivative operator
    sp_mat N = nodal2D(k, m, 1.0, n, 1.0);

    // Apply operator
    vec Xd = N * X;
    vec Yd = N * Y;

    int mn = m * n;

    Xe = Xd.subvec(0, mn-1);
    Xn = Xd.subvec(mn, 2*mn-1);
    Ye = Yd.subvec(0, mn-1);
    Yn = Yd.subvec(mn, 2*mn-1);

    // Compute the Jacobian determinant
    J = Xe % Yn - Xn % Ye;
}

sp_mat DI2(int m, int n, const std::string& type) {
    if (type == "Dn") {
        int rows = (m+1)*n;
        int cols = (m+1)*n;

        // Prepare triplet lists
        std::vector<uword> row_indices;
        std::vector<uword> col_indices;
        std::vector<double> values;

        // Build bdry manually
        for (int i = 0; i < m; ++i) {
            row_indices.push_back(i);
            col_indices.push_back(i);
            values.push_back(-0.5);

            row_indices.push_back(i);
            col_indices.push_back(i+1);
            values.push_back(-0.5);

            row_indices.push_back(i);
            col_indices.push_back(i+(m+1));
            values.push_back(0.5);

            row_indices.push_back(i);
            col_indices.push_back(i+(m+2));
            values.push_back(0.5);
        }

        sp_mat bdry(rows, cols);
        for (size_t idx = 0; idx < row_indices.size(); ++idx) {
            bdry(row_indices[idx], col_indices[idx]) = values[idx];
        }

        // Build block
        sp_mat block(m+2, m+1);
        for (int i = 0; i < m; ++i) {
            block(i, i) += 0.25;
            block(i, i+1) += 0.25;
        }
        block(m, m) += 0.25;  // extra for size

        // Build pattern
        sp_mat pattern(n-2, n);
        for (int i = 0; i < n-2; ++i) {
            pattern(i, i) = -1.0;
            pattern(i, i+2) = 1.0;
        }

        sp_mat middle = kron(pattern, block);

        // Shift bdry
        sp_mat bdry_shifted = bdry;
        bdry_shifted = shift(bdry_shifted, (m+1)*(n-2), 1);  // shift columns

        // Assemble I
        sp_mat top(m+3, cols);
        sp_mat mid(2, cols);
        sp_mat bottom(m+3, cols);

        sp_mat I = join_cols(join_cols(join_cols(join_cols(top, bdry), mid), middle), join_cols(bdry_shifted, bottom));

        return I;
    }
    else { // "De" case
        int rows = (n+1)*m;
        int cols = (n+1)*m;

        // Build block
        sp_mat block(m+2, m);
        for (int i = 1; i < m-1; ++i) {
            block(i, i-1) += -0.25;
            block(i, i)   +=  0.25;
        }
        // Manually fix boundary
        block(0, 0) = -0.5;
        block(0, 1) =  0.5;
        block(m-1, m-2) = -0.5;
        block(m-1, m-1) =  0.5;

        // Build pattern
        sp_mat pattern(n, n+1);
        for (int i = 0; i < n; ++i) {
            pattern(i, i) = 1.0;
            pattern(i, i+1) = 1.0;
        }

        sp_mat middle = kron(pattern, block);

        sp_mat top(m+3, cols);
        sp_mat bottom(m+1, cols);

        sp_mat I = join_cols(join_cols(top, middle), bottom);

        return I;
    }
}

sp_mat make_diag_sparse(const mat& M) {
    uword N = M.n_elem;
    vec vals = vectorise(M.t());
    umat locations(2, N);
    locations.row(0) = regspace<urowvec>(0, N-1);
    locations.row(1) = regspace<urowvec>(0, N-1);
    return sp_mat(locations, vals, N, N);
}

sp_mat div2DCurv(int k, const mat& X_in, const mat& Y_in) {
    // Get the Jacobian determinant and metrics
    vec J_vec, Xe_vec, Xn_vec, Ye_vec, Yn_vec;
    jacobian2D(k, X_in, Y_in, J_vec, Xe_vec, Xn_vec, Ye_vec, Yn_vec);
    
    std::cout << "size(J_vec) = " << arma::size(J_vec) << '\n';
    std::cout << "size(Xe_vec) = " << arma::size(Xe_vec) << '\n';
    std::cout << "size(Xn_vec) = " << arma::size(Xn_vec) << '\n';
    std::cout << "size(Ye_vec) = " << arma::size(Ye_vec) << '\n';
    std::cout << "size(Yn_vec) = " << arma::size(Yn_vec) << '\n';

    int n = X_in.n_rows;
    int m = X_in.n_cols;

    // Reshape to (n, m) matrices
    mat J = reshape(J_vec, m, n).t();
    mat Xe = reshape(Xe_vec, m, n).t();
    mat Xn = reshape(Xn_vec, m, n).t();
    mat Ye = reshape(Ye_vec, m, n).t();
    mat Yn = reshape(Yn_vec, m, n).t();
    
    std::cout << "size(J) = " << arma::size(J) << '\n';
    std::cout << "size(Xe) = " << arma::size(Xe) << '\n';
    std::cout << "size(Xn) = " << arma::size(Xn) << '\n';
    std::cout << "size(Ye) = " << arma::size(Ye) << '\n';
    std::cout << "size(Yn) = " << arma::size(Yn) << '\n';

    // Logical grid
    vec Xl = linspace(1, m, m);
    vec Yl = linspace(1, n, n);

    // Staggered logical grid
    vec Xs = join_vert(vec({1.0}), regspace(1.5, 1.0, m-0.5), vec({(double)m}));
    vec Ys = join_vert(vec({1.0}), regspace(1.5, 1.0, n-0.5), vec({(double)n}));
    
    std::cout << "size(Xl) = " << arma::size(Xl) << '\n';
    std::cout << "size(Yl) = " << arma::size(Yl) << '\n';
    std::cout << "size(Xs) = " << arma::size(Xs) << '\n';
    std::cout << "size(Ys) = " << arma::size(Ys) << '\n';

    // Interpolate metrics onto staggered grid
    mat J_interp;
    interp2(Xl, Yl, J, Xs, Ys, J_interp, "linear", datum::nan);

    mat Xe_interp;
    interp2(Xl, Yl, Xe, Xs, Ys, Xe_interp, "linear", datum::nan);

    mat Xn_interp;
    interp2(Xl, Yl, Xn, Xs, Ys, Xn_interp, "linear", datum::nan);

    mat Ye_interp;
    interp2(Xl, Yl, Ye, Xs, Ys, Ye_interp, "linear", datum::nan);

    mat Yn_interp;
    interp2(Xl, Yl, Yn, Xs, Ys, Yn_interp, "linear", datum::nan);
    
    vec xgrid = linspace(1.0, (double)m, m);
    vec ygrid = linspace(1.0, (double)n, n);

    std::cout << "size(J_interp) = " << arma::size(J_interp) << '\n';
    std::cout << "size(Xe_interp) = " << arma::size(Xe_interp) << '\n';
    std::cout << "size(Xn_interp) = " << arma::size(Xn_interp) << '\n';
    std::cout << "size(Ye_interp) = " << arma::size(Ye_interp) << '\n';
    std::cout << "size(Yn_interp) = " << arma::size(Yn_interp) << '\n';
    
    sp_mat J_inv_diag = make_diag_sparse(1.0 / J_interp);
    sp_mat Xe_diag = make_diag_sparse(Xe_interp);
    sp_mat Xn_diag = make_diag_sparse(Xn_interp);
    sp_mat Ye_diag = make_diag_sparse(Ye_interp);
    sp_mat Yn_diag = make_diag_sparse(Yn_interp);
    
    std::cout << "size(J_inv_diag) = " << arma::size(J_inv_diag) << '\n';
    std::cout << "size(Xe_diag) = " << arma::size(Xe_diag) << '\n';
    std::cout << "size(Xn_diag) = " << arma::size(Xn_diag) << '\n';
    std::cout << "size(Ye_diag) = " << arma::size(Ye_diag) << '\n';
    std::cout << "size(Yn_diag) = " << arma::size(Yn_diag) << '\n';

    // Build 2D uniform mimetic divergence operator
    Divergence D(k, m-1, n-1, 1.0, 1.0);
    
    std::cout << "size(D) = " << arma::size(D) << '\n';

    sp_mat De = D.cols(0, m*(n-1)-1);
    sp_mat Dn = D.cols(m*(n-1), D.n_cols-1);
    
    std::cout << "size(De) = " << arma::size(De) << '\n';
    std::cout << "size(Dn) = " << arma::size(Dn) << '\n';

    // Build logical-to-physical divergence operators
    sp_mat Dx = J_inv_diag * (Yn_diag * De - Ye_diag * DI2(m-1, n-1, "Dn"));
    sp_mat Dy = J_inv_diag * (-Xn_diag * DI2(m-1, n-1, "De") + Xe_diag * Dn);

    // Final 2D curvilinear mimetic divergence operator
    sp_mat D_final = join_rows(Dx, Dy);

    return D_final;
}

int main()
{
    constexpr unsigned int e = 1; // Controls "rectangularity" of the grid, e = 0 -> completely rectangular
    constexpr unsigned int k = 2;
    unsigned int m = 15;
    unsigned int n = 20;
    constexpr double a = -pi;
    constexpr double b = 2*pi;
    constexpr double c = -pi;
    constexpr double d = pi;
    mat K = {{1, 0}, {0, 1}};
    
//    auto F = [&](const auto& X, const auto& Y)
//    {
//        return sin(X)*sin(Y);
//    };
    
    double dm = (b-a)/static_cast<double>(m);
    double dn = (d-c)/static_cast<double>(n);
    
    mat X, Y;
    
    auto f = [&](const auto& X, const auto& Y)
    {
//        auto fx = K(0, 0)*diff(F(X, Y), 1, 0);
//        auto fy = K(1, 1)*diff(F(X, Y), 1, 1);
//        return fx + fy;
        return -sin(X)*sin(Y)*(K(0, 0) + K(1, 1)); //the explicit symbolic result of the following:
        /*
         F = @(X, Y) sin(X).*sin(Y);
         syms X Y
         fx = K(1)*diff(F(X, Y), X);
         fy = K(4)*diff(F(X, Y), Y);
         f = matlabFunction(diff(fx, X)+diff(fy, Y));
         */
    };
    
    vec first = linspace(a, b, 1+(b-a)/dm);
    vec second = linspace(c, d, 1+(d-c)/dn);
//    std::cout << first << '\n';
//    std::cout << second << '\n';
    Utils obj;
    obj.meshgrid(first, second, X, Y);
//    std::cout << X << '\n';
//    std::cout << Y << '\n';
    Y += e*cos(X);
//    std::cout << Y << '\n';
    n = X.n_rows;
    m = X.n_cols;
    n = n-1;
    m = m-1;
//    printf("m = %d, n = %d\n", m, n);
    mat Ux = (X.rows(0, X.n_rows - 2) + X.rows(1, X.n_rows - 1)) / 2.0;
    mat Uy = (Y.rows(0, Y.n_rows - 2) + Y.rows(1, Y.n_rows - 1)) / 2.0;
//    std::cout << Ux.n_rows << ',' << Ux.n_cols << '\n';
//    std::cout << Uy.n_rows << ',' << Uy.n_cols << '\n';
    mat Z = zeros(n+1, m+1);
    
    std::ofstream plot_script("plot.gnu");
    if (!plot_script)
    {
        std::cerr << "Error: Failed to create GNUplot script.\n";
        return 1;
    }

//    plot_script << "set terminal qt size 800,800\n";
    plot_script << "set title 'Grid'\n";
    plot_script << "set size ratio 1\n";
    plot_script << "set title 'Grid' offset 0, -5\n";
//    plot_script << "unset key\n";
    plot_script << "set size 1.0, 1.0\n";
    plot_script << "unset xlabel\n";
    plot_script << "unset ylabel\n";
    plot_script << "set view 0, 90, 1, 1\n";  // Set a 3D view
    plot_script << "set xrange [" << Y.max() << ':' << Y.min() << "];\n";
    plot_script << "splot '-' using 1:2:3 title 'Nodal points' with linespoints pt 7 ps 0.35 lc rgb '#6ce4dc', "
                << "'-' using 1:2:3 title 'u' with points pt 1 ps 0.75 lc rgb 'black', "
                << "'-' using 1:2:3 title 'v' with points pt 2 ps 0.75 lc rgb 'black', "
                << "'-' using 1:2:3 title 'Centers' with points pt 5 ps 0.25 lc rgb 'red', "
                << "'-' using 1:2:3 title 'All centers' with points pt 6 ps 1 lc rgb 'red'\n";

//    splot '-' using 1:2:3 with linespoints pt 7 ps 0.5
    
    // Write the grid data points
    for (unsigned int j = 0; j <= m; ++j) {
        for (unsigned int i = 0; i <= n; ++i) {
            plot_script << Y(i, j) << " " << X(i, j) << " " << Z(i, j) << "\n";
        }
        plot_script << '\n';  // Blank line to separate different j values
    }
    plot_script << "e\n";  // End of first dataset

    // Write the points
    for (unsigned int j = 0; j < Uy.n_cols; ++j) {
        for (unsigned int i = 0; i < Uy.n_rows; ++i) {
            plot_script << Uy(i, j) << " " << Ux(i, j) << " " << 0.0 << "\n";
        }
        plot_script << '\n';  // Blank line to separate different j values
    }
    
    plot_script << "e\n";  // End of data
    
    mat Vx = (X.cols(0, X.n_cols - 2) + X.cols(1, X.n_cols - 1)) / 2.0;
    mat Vy = (Y.cols(0, Y.n_cols - 2) + Y.cols(1, Y.n_cols - 1)) / 2.0;
    
    // Write the points
    for (unsigned int j = 0; j < Vy.n_cols; ++j) {
        for (unsigned int i = 0; i < Vy.n_rows; ++i) {
            plot_script << Vy(i, j) << " " << Vx(i, j) << " " << 0.0 << "\n";
        }
        plot_script << '\n';  // Blank line to separate different j values
    }
    
    plot_script << "e\n";  // End of data
    
    mat Cx = (Vx.rows(0, Vx.n_rows - 2) + Vx.rows(1, Vx.n_rows - 1)) / 2.0;
    mat Cy = (Uy.cols(0, Uy.n_cols - 2) + Uy.cols(1, Uy.n_cols - 1)) / 2.0;
    
    // Write the points
    for (unsigned int j = 0; j < Cy.n_cols; ++j) {
        for (unsigned int i = 0; i < Cy.n_rows; ++i) {
            plot_script << Cy(i, j) << " " << Cx(i, j) << " " << 0.0 << "\n";
        }
        plot_script << '\n';  // Blank line to separate different j values
    }
    
    plot_script << "e\n";  // End of data
    
    // West-East sides
    Cx = join_rows(Ux.col(0), Cx);
    Cy = join_rows(Uy.col(0), Cy);
    Cx = join_rows(Cx, Ux.col(Ux.n_cols - 1));
    Cy = join_rows(Cy, Uy.col(Uy.n_cols - 1));
    
    // South-North sides
    mat new_row_top_Cx = join_rows(zeros<mat>(1, 1), Vx.row(0), zeros<mat>(1, 1));
    mat new_row_top_Cy = join_rows(zeros<mat>(1, 1), Vy.row(0), zeros<mat>(1, 1));
    Cx = join_cols(new_row_top_Cx, Cx);
    Cy = join_cols(new_row_top_Cy, Cy);
    mat new_row_bottom_Cx = join_rows(zeros<mat>(1, 1), Vx.row(Vx.n_rows - 1), zeros<mat>(1, 1));
    mat new_row_bottom_Cy = join_rows(zeros<mat>(1, 1), Vy.row(Vy.n_rows - 1), zeros<mat>(1, 1));
    Cx = join_cols(Cx, new_row_bottom_Cx);
    Cy = join_cols(Cy, new_row_bottom_Cy);
    
    // Corners
    Cx(0, 0) = X(0, 0);
    Cy(0, 0) = Y(0, 0);
    Cx(0, Cx.n_cols - 1) = X(0, X.n_cols - 1);
    Cy(0, Cy.n_cols - 1) = Y(0, Y.n_cols - 1);
    Cx(Cx.n_rows - 1, 0) = X(X.n_rows - 1, 0);
    Cy(Cy.n_rows - 1, 0) = Y(Y.n_rows - 1, 0);
    Cx(Cx.n_rows - 1, Cx.n_cols - 1) = X(X.n_rows - 1, X.n_cols - 1);
    Cy(Cy.n_rows - 1, Cy.n_cols - 1) = Y(Y.n_rows - 1, Y.n_cols - 1);
    
    // Write the points
    for (unsigned int j = 0; j < Cy.n_cols; ++j) {
        for (unsigned int i = 0; i < Cy.n_rows; ++i) {
            plot_script << Cy(i, j) << " " << Cx(i, j) << " " << 0.0 << "\n";
            //std::cout << Cy(i, j) << " " << Cx(i, j) << " " << 0.0 << "\n";
        }
        plot_script << '\n';  // Blank line to separate different j values
        puts("");
    }
    
    plot_script << "e\n";  // End of data
    
    plot_script.close();

    // Execute the Gnuplot script
    if (system("gnuplot -persist plot.gnu") != 0) {
        std::cerr << "Error : Failed to execute GNUplot.\n";
        return 1;
    }
    
    sp_mat D = div2DCurv(k, X, Y);
    std::cout << "D = " << D << '\n';
    
    return 0;
    

}

///opt/homebrew/Cellar/llvm/20.1.2/bin/clang++ -std=c++20 -O3 -fopenmp -DARMA_DONT_USE_WRAPPER -Wall -DARMA_SUPERLU_INCLUDE_DIR=\"/opt/homebrew/Cellar/superlu/7.0.0/include\"  -I.  -I/opt/homebrew/Cellar/armadillo/14.4.1/include -I../../src/cpp -I/opt/homebrew/Cellar/superlu/7.0.0/include -o elliptic2D_case2 elliptic2D_case2.cpp -L../../src/cpp -lmole -L/opt/homebrew/Cellar/openblas/0.3.29/lib -L/opt/homebrew/Cellar/armadillo/14.4.1/lib -L/opt/homebrew/Cellar/superlu/7.0.0/lib -lopenblas -lsuperlu -L/opt/homebrew/opt/superlu/lib -L/opt/homebrew/Cellar/armadillo/14.4.1 -Wl,-rpath,/opt/homebrew/Cellar/armadillo/14.4.1 -larmadillo -L/opt/homebrew/Cellar/superlu/7.0.0/lib -Wl,-rpath,/opt/homebrew/Cellar/superlu/7.0.0/lib -lsuperlu
