#include "catch.hpp"

#include "libirc/linalg.h"

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
#else
#error
#endif

using namespace irc;

TEST_CASE("Size", "[size]"){
    SECTION("Vector"){
        vec v = {1, 2, 3, 4, 5};
        
        CHECK(linalg::size(v) == 5);
    }
    
    SECTION("Matrix"){
        mat m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        
        CHECK(linalg::size(m) == 9);
    }
}

TEST_CASE("Number of rows", "[nrows]"){
    SECTION("Vector"){
        vec v = {1, 2, 3, 4, 5};
        
        CHECK(linalg::n_rows(v) == 5);
    }
    
    SECTION("Matrix"){
        mat m = {{1, 2, 3}, {4, 5, 6}};

        CHECK(linalg::n_rows(m) == 2);
    }
}

TEST_CASE("Number of columns", "[ncols]"){
    SECTION("Vector"){
        vec v = {1, 2, 3, 4, 5};
    
        CHECK(linalg::n_cols(v) == 1);
    }

    SECTION("Matrix"){
        mat m = {{1, 2, 3}, {4, 5, 6}};

        CHECK(linalg::n_cols(m) == 3);
    }
}

TEST_CASE("Norm and normalization", "[norm]"){
    SECTION("Vector"){
        vec v = {1, 2, 3, 4};
        
        Approx t(std::sqrt(1*1 + 2*2 + 3*3 + 4*4));
        CHECK(linalg::norm(v) == t);
        
        v = linalg::normalize(v);
        
        CHECK(linalg::norm(v) == Approx(1.0));
    }

    SECTION("Matrix"){
        mat m = {{1, 2, 3}, {4, 5, 6}};

        Approx t(std::sqrt(0.5 * (91 + std::sqrt(8065)))); // Obtained from WolframAlpha
        CHECK(linalg::norm(m) == t);
    }
}