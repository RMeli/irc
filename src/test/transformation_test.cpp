#include "catch.hpp"

#include "libirc/transformation.h"

#include "config.h"
#include "libirc/conversion.h"
#include "libirc/io.h"
#include "libirc/irc.h"
#include "libirc/molecule.h"
#include "libirc/wilson.h"

#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;

template<typename T>
using Mat = arma::Mat<T>;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;

template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
#else
#error
#endif

using namespace irc;

TEST_CASE("RMS") {
  using namespace transformation;

  vec v{1, 2, 3, 4};

  CHECK(rms(v) == Approx(std::sqrt(30 / 4.)));
}

TEST_CASE("Transformation") {

  bool verbose{true};

  SECTION("Cartesian to internal for ethanol") {
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;

    // Load molecule from file
    const auto molecule = load_xyz<vec3>(config::molecules_dir + "ethanol.xyz");

    // Compute interatomic distance for formaldehyde molecule
    mat dd{distances<vec3, mat>(molecule)};

    // Build graph based on the adjacency matrix
    UGraph adj{adjacency_matrix(dd, molecule)};

    // Compute distance matrix and predecessor matrix
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, molecule)};

    // Print bonds
    if (verbose) {
      print_bonds<vec3, vec>(to_cartesian<vec3, vec>(molecule), B);
    }

    // Compute angles
    std::vector<Angle> A{angles(dist, molecule)};

    // Print angles
    if (verbose) {
      print_angles<vec3, vec>(to_cartesian<vec3, vec>(molecule), A);
    }

    // Compute dihedral angles
    std::vector<Dihedral> D{dihedrals(dist, molecule)};

    // Print dihedral angles
    if (verbose) {
      print_dihedrals<vec3, vec>(to_cartesian<vec3, vec>(molecule), D);
    }

    // Compute linear angles
    std::vector<LinearAngle<vec3>> LA{linear_angles(dist, molecule)};

    // Print linear angles
    if (verbose) {
      print_linear_angles<vec3, vec>(to_cartesian<vec3, vec>(molecule), LA);
    }

    // Compute linear angles
    std::vector<OutOfPlaneBend> OOPB{out_of_plane_bends(dist, molecule)};

    // Print linear angles
    if (verbose) {
      print_out_of_plane_bends<vec3, vec>(to_cartesian<vec3, vec>(molecule),
                                          OOPB);
    }

    // Compute number of cartesian coordinates
    std::size_t n_c{3 * molecule.size()};

    // Allocate vector for cartesian positions
    vec x_c{linalg::zeros<vec>(n_c)};

    // Fill vector with cartesian positions
    for (std::size_t i{0}; i < molecule.size(); i++) {
      x_c(3 * i + 0) = molecule[i].position(0);
      x_c(3 * i + 1) = molecule[i].position(1);
      x_c(3 * i + 2) = molecule[i].position(2);
    }

    if (verbose) {
      // Print cartesian coordinates
      cout << "\nCartesian coordinates (a.u.):\n " << x_c << endl;
    }

    if (verbose) {
      // Compute and print internal redundant coordinates
      cout << "Internal redundant coordinates (a.u.):\n"
           << cartesian_to_irc<vec3, vec>(x_c, B, A, D, LA, OOPB) << endl;
    }
  }

  SECTION("Internal to cartesian for H2") {
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    using namespace wilson;

    // Load molecule from file
    Molecule<vec3> molecule{{{"H", {0., 0., 0.}}, {"H", {1., 0., 0.}}}};

    // Connectivity
    mat dd{distances<vec3, mat>(molecule)};
    UGraph adj{adjacency_matrix(dd, molecule)};
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, molecule)};

    // Print bonds
    if (verbose) {
      print_bonds<vec3, vec>(to_cartesian<vec3, vec>(molecule), B);
    }

    // Check number of bonds
    REQUIRE(B.size() == 1);

    // Wilson B matrix
    mat W = wilson_matrix<vec3, vec, mat>(
        molecule::to_cartesian<vec3, vec>(molecule), B);

    // Allocate vector for internal reaction coordinates
    vec q_irc{connectivity::cartesian_to_irc<vec3, vec>(
        molecule::to_cartesian<vec3, vec>(molecule), B, {}, {}, {}, {})};

    // Displacement in internal coordinates
    vec dq_irc{0.1};

    // Print displacement in internal coordinates
    if (verbose) {
      cout << "\nDisplacement in internal coordinates (a.u.):\n " << dq_irc
           << endl;
    }

    // Compute number of cartesian coordinates
    std::size_t n_c{3 * molecule.size()};

    // Allocate vector for cartesian positions
    vec x_c_old{linalg::zeros<vec>(n_c)};

    // Fill vector with cartesian positions
    for (std::size_t i{0}; i < molecule.size(); i++) {
      x_c_old(3 * i + 0) = molecule[i].position(0);
      x_c_old(3 * i + 1) = molecule[i].position(1);
      x_c_old(3 * i + 2) = molecule[i].position(2);
    }

    // Compute new cartesian coordinates
    const auto itc_result = irc_to_cartesian<vec3, vec, mat>(
        q_irc, dq_irc, x_c_old, B, {}, {}, {}, {});
    const auto x_c = itc_result.x_c;

    // Print cartesian coordinates
    if (verbose) {
      cout << "\nNew cartesian coordinates (a.u.):\n " << x_c << endl;
    }

    vec3 p1{x_c(0), x_c(1), x_c(2)};
    vec3 p2{x_c(3), x_c(4), x_c(5)};

    Approx target{q_irc(0) + dq_irc(0)};
    target.margin(1e-6);

    REQUIRE(distance(p1, p2) == target);
  }

  SECTION("Internal to cartesian for H2O") {
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    using namespace wilson;

    // Load molecule from file
    const auto molecule = load_xyz<vec3>(config::molecules_dir + "water.xyz");

    // Compute interatomic distance for formaldehyde molecule
    mat dd{distances<vec3, mat>(molecule)};

    // Build graph based on the adjacency matrix
    UGraph adj{adjacency_matrix(dd, molecule)};

    // Compute distance matrix and predecessor matrix
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, molecule)};

    // Print bonds
    if (verbose) {
      print_bonds<vec3, vec>(to_cartesian<vec3, vec>(molecule), B);
    }

    // Check number of bonds
    REQUIRE(B.size() == 2);

    // Compute angle
    std::vector<Angle> A{angles(dist, molecule)};

    // Print angles
    if (verbose) {
      print_angles<vec3, vec>(to_cartesian<vec3, vec>(molecule), A);
    }

    // Compute Wilson B matrix
    mat W = wilson_matrix<vec3, vec, mat>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A);

    // Allocate vector for internal reaction coordinates
    vec q_irc_old{connectivity::cartesian_to_irc<vec3, vec>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A, {}, {}, {})};

    // Displacement in internal coordinates
    vec dq_irc{0.0, 0.0, 1. / 180. * tools::constants::pi};

    // Compute new internal coordinates
    vec q_irc_new{q_irc_old + dq_irc};

    // Print new internal coordinates
    if (verbose) {
      cout << "\nNew internal coordinates:\n " << q_irc_new << endl;
    }

    // Allocate vector for cartesian positions
    vec x_c_old{to_cartesian<vec3, vec>(molecule)};

    // Compute new cartesian coordinates
    const auto itc_result = irc_to_cartesian<vec3, vec, mat>(
        q_irc_old, dq_irc, x_c_old, B, A, {}, {}, {});
    const auto x_c = itc_result.x_c;

    // Print cartesian coordinates
    if (verbose) {
      cout << "\nNew cartesian coordinates (a.u.):\n " << x_c << endl;
    }

    // Reconstruct atomic positions (points in 3D space)
    vec3 p1{x_c(0), x_c(1), x_c(2)};
    vec3 p2{x_c(3), x_c(4), x_c(5)};
    vec3 p3{x_c(6), x_c(7), x_c(8)};

    SECTION("Bond 1") {
      Approx target{q_irc_new(0)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p2) == target);
    }

    SECTION("Bond 2") {
      Approx target{q_irc_new(1)};
      target.margin(1e-4);
      REQUIRE(distance(p2, p3) == target);
    }

    SECTION("Angle") {
      Approx target{q_irc_new(2)};
      target.margin(1e-4);
      REQUIRE(angle(p1, p2, p3) == target);
    }
  }

  SECTION("Internal to cartesian for H2O2") {
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    using namespace wilson;

    // Load molecule from file
    const auto molecule =
        load_xyz<vec3>(config::molecules_dir + "hydrogen_peroxide.xyz");

    // Compute interatomic distance for formaldehyde molecule
    mat dd{distances<vec3, mat>(molecule)};

    // Build graph based on the adjacency matrix
    UGraph adj{adjacency_matrix(dd, molecule)};

    // Compute distance matrix and predecessor matrix
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, molecule)};

    // Print bonds
    if (verbose) {
      print_bonds<vec3, vec>(to_cartesian<vec3, vec>(molecule), B);
    }

    // Check number of bonds
    REQUIRE(B.size() == 3);

    // Compute angles
    std::vector<Angle> A{angles(dist, molecule)};

    // Print angles
    if (verbose) {
      print_angles<vec3, vec>(to_cartesian<vec3, vec>(molecule), A);
    }

    // Check number of angles
    REQUIRE(A.size() == 2);

    // Compute dihedral angles
    std::vector<Dihedral> D{dihedrals(dist, molecule)};

    // Print dihedrals
    if (verbose) {
      print_dihedrals<vec3, vec>(to_cartesian<vec3, vec>(molecule), D);
    }

    // Check number of dihedrals
    REQUIRE(D.size() == 1);

    // Compute Wilson B matrix
    mat W = wilson_matrix<vec3, vec, mat>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A, D);

    // Allocate vector for internal reaction coordinates
    vec q_irc_old{connectivity::cartesian_to_irc<vec3, vec>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A, D, {}, {})};

    // Displacement in internal coordinates
    vec dq_irc{0.0, 0.0, 0.0, 0.0, 0.0, 1. / 180 * tools::constants::pi};

    // Compute new internal coordinates
    vec q_irc_new{q_irc_old + dq_irc};

    // Print new internal coordinates
    if (verbose) {
      cout << "\nNew internal coordinates:\n " << q_irc_new << endl;
    }

    // Allocate vector for cartesian positions
    vec x_c_old{to_cartesian<vec3, vec>(molecule)};

    // Compute new cartesian coordinates
    const auto itc_result = irc_to_cartesian<vec3, vec, mat>(
        q_irc_old, dq_irc, x_c_old, B, A, D, {}, {});
    const auto x_c = itc_result.x_c;

    // Print cartesian coordinates
    if (verbose) {
      cout << "\nNew cartesian coordinates (a.u.):\n " << x_c << endl;
    }

    // Reconstruct atomic positions (points in 3D space)
    vec3 p1{x_c(0), x_c(1), x_c(2)};
    vec3 p2{x_c(3), x_c(4), x_c(5)};
    vec3 p3{x_c(6), x_c(7), x_c(8)};
    vec3 p4{x_c(9), x_c(10), x_c(11)};

    SECTION("Bond 1") {
      Approx target{q_irc_new(0)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p2) == target);
    }

    SECTION("Bond 2") {
      Approx target{q_irc_new(1)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p3) == target);
    }

    SECTION("Bond 3") {
      Approx target{q_irc_new(2)};
      target.margin(1e-4);
      REQUIRE(distance(p2, p4) == target);
    }

    SECTION("Angle 1") {
      Approx target{q_irc_new(3)};
      target.margin(1e-4);
      REQUIRE(angle(p2, p1, p3) == target);
    }

    SECTION("Angle 2") {
      Approx target{q_irc_new(4)};
      target.margin(1e-4);
      REQUIRE(angle(p1, p2, p4) == target);
    }

    SECTION("Dihedral") {
      Approx target{q_irc_new(5)};
      target.margin(1e-4);
      REQUIRE(dihedral(p4, p2, p1, p3) == target);
    }
  }

  SECTION("Internal to cartesian for CO2") {
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    using namespace wilson;

    // Load molecule from file
    const auto molecule =
        load_xyz<vec3>(config::molecules_dir + "carbon_dioxide.xyz");

    // Compute interatomic distance for molecule
    mat dd{distances<vec3, mat>(molecule)};

    // Build graph based on the adjacency matrix
    UGraph adj{adjacency_matrix(dd, molecule)};

    // Compute distance matrix and predecessor matrix
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, molecule)};

    // Print bonds
    if (verbose) {
      print_bonds<vec3, vec>(to_cartesian<vec3, vec>(molecule), B);
    }

    // Check number of bonds
    REQUIRE(B.size() == 2);

    // Compute angles
    std::vector<Angle> A{angles(dist, molecule)};

    // Print angles
    if (verbose) {
      print_angles<vec3, vec>(to_cartesian<vec3, vec>(molecule), A);
    }

    // Check number of angles
    REQUIRE(A.size() == 0);

    // Compute dihedral angles
    std::vector<Dihedral> D{dihedrals(dist, molecule)};

    // Print dihedrals
    if (verbose) {
      print_dihedrals<vec3, vec>(to_cartesian<vec3, vec>(molecule), D);
    }

    // Check number of dihedrals
    REQUIRE(D.size() == 0);

    // Compute linear angles
    std::vector<LinearAngle<vec3>> LA{linear_angles(dist, molecule)};

    // Print linear angles
    if (verbose) {
      print_linear_angles<vec3, vec>(to_cartesian<vec3, vec>(molecule), LA);
    }

    // Check number of linear angles
    REQUIRE(LA.size() == 2);

    // Compute Wilson B matrix
    mat W = wilson_matrix<vec3, vec, mat>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A, D, LA, {});

    // Allocate vector for internal reaction coordinates
    vec q_irc_old{connectivity::cartesian_to_irc<vec3, vec>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A, D, LA, {})};

    // Displacement in internal coordinates
    vec dq_irc{0.0, 0.0, 0.0, 1. / 180 * tools::constants::pi};

    // Compute new internal coordinates
    vec q_irc_new{q_irc_old + dq_irc};

    // Print new internal coordinates
    if (verbose) {
      cout << "\nNew internal coordinates:\n " << q_irc_new << endl;
    }

    // Allocate vector for cartesian positions
    vec x_c_old{to_cartesian<vec3, vec>(molecule)};

    // Compute new cartesian coordinates
    const auto itc_result = irc_to_cartesian<vec3, vec, mat>(
        q_irc_old, dq_irc, x_c_old, B, A, D, LA, {});
    const auto x_c = itc_result.x_c;

    // Print cartesian coordinates
    if (verbose) {
      cout << "\nNew cartesian coordinates (a.u.):\n " << x_c << endl;
    }

    // Reconstruct atomic positions (points in 3D space)
    vec3 p1{x_c(0), x_c(1), x_c(2)};
    vec3 p2{x_c(3), x_c(4), x_c(5)};
    vec3 p3{x_c(6), x_c(7), x_c(8)};

    SECTION("Bond CO{1}") {
      Approx target{q_irc_new(0)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p2) == target);
    }

    SECTION("Bond CO{2}") {
      Approx target{q_irc_new(1)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p3) == target);
    }

    SECTION("Angle OCO") {
      Approx target{q_irc_new(3)};
      target.margin(1e-4);
      REQUIRE(angle(LA[1], x_c) == target);
    }
  }
}

TEST_CASE("Issue 41") {
  SECTION("Big change in water") {
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    using namespace wilson;

    auto molecule = molecule::Molecule<vec3>{
        {{"O", {0., 0., 0.}}, {"H", {1., 0., 0.}}, {"H", {0., 1., 0.}}}};
    multiply_positions(molecule, tools::conversion::angstrom_to_bohr);

    // Compute interatomic distance for formaldehyde molecule
    mat dd{distances<vec3, mat>(molecule)};

    // Build graph based on the adjacency matrix
    UGraph adj{adjacency_matrix(dd, molecule)};

    // Compute distance matrix and predecessor matrix
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, molecule)};

    // Check number of bonds
    REQUIRE(B.size() == 2);

    // Compute angle
    std::vector<Angle> A{angles(dist, molecule)};

    // Compute Wilson B matrix
    mat W = wilson_matrix<vec3, vec, mat>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A);

    // Allocate vector for internal reaction coordinates
    vec q_irc_old{connectivity::cartesian_to_irc<vec3, vec>(
        molecule::to_cartesian<vec3, vec>(molecule), B, A, {}, {}, {})};

    // Displacement in internal coordinates
    vec dq_irc_75{0.5, 0.5, 75. / 180. * tools::constants::pi};
    vec dq_irc_89{0.5, 0.5, 89. / 180. * tools::constants::pi};

    // Compute new internal coordinates
    vec q_irc_new{q_irc_old + dq_irc_75};
    CAPTURE(q_irc_new);

    // Allocate vector for cartesian positions
    vec x_c_old{to_cartesian<vec3, vec>(molecule)};

    // Compute new cartesian coordinates
    const auto itc_result_single = irc_to_cartesian_single<vec3, vec, mat>(
        q_irc_old, dq_irc_75, x_c_old, B, A, {}, {}, {});
    CHECK(!itc_result_single.converged);

    const auto itc_result = irc_to_cartesian<vec3, vec, mat>(
        q_irc_old, dq_irc_75, x_c_old, B, A, {}, {}, {});
    CHECK(itc_result.converged);
    const auto x_c = itc_result.x_c;

    const auto itc_result89 = irc_to_cartesian<vec3, vec, mat>(
        q_irc_old, dq_irc_89, x_c_old, B, A, {}, {}, {});
    CHECK(itc_result89.converged);

    // Print cartesian coordinates
    CAPTURE(x_c);

    // Reconstruct atomic positions (points in 3D space)
    vec3 p1{x_c(0), x_c(1), x_c(2)};
    vec3 p2{x_c(3), x_c(4), x_c(5)};
    vec3 p3{x_c(6), x_c(7), x_c(8)};

    SECTION("Bond 1") {
      Approx target{q_irc_new(0)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p2) == target);
    }

    SECTION("Bond 2") {
      Approx target{q_irc_new(1)};
      target.margin(1e-4);
      REQUIRE(distance(p1, p3) == target);
    }

    SECTION("Angle") {
      Approx target{q_irc_new(2)};
      target.margin(1e-4);
      REQUIRE(angle(p2, p1, p3) == target);
    }
  }
}

TEST_CASE("Problems") {
  using namespace molecule;
  using namespace connectivity;
  using namespace transformation;
  using namespace tools::conversion;
  using namespace io;
  using namespace std;
  const std::string initial = "22                                                                               \n"
                              "                                                                                 \n"
                              "C      2.1180331035    0.6524228183    0.6375019190                              \n"
                              "C      1.6441147713    0.5728544661   -0.7089733020                              \n"
                              "C      0.3455255130    0.9969960630   -0.9865779900                              \n"
                              "C     -0.4928024308    1.5311634237    0.0779040585                              \n"
                              "C      0.0268787843    1.6422330060    1.4241830598                              \n"
                              "C      1.3383928814    1.1860121740    1.6726925298                              \n"
                              "C      3.4965764670    0.1257190961    0.9030552962                              \n"
                              "C      4.6959764678    0.8536394010    0.5212215417                              \n"
                              "C      5.4988467330    0.3563345564   -0.4663595922                              \n"
                              "C      5.0344184935   -0.9033654323   -1.1556157485                              \n"
                              "C      4.4134390629   -1.8701179927   -0.1932514017                              \n"
                              "C      3.6653961780   -1.3206891625    0.8105986697                              \n"
                              "O      4.7375270679   -0.8750792705   -2.3528797355                              \n"
                              "O     -1.6899719508    1.8941139777   -0.1913451212                              \n"
                              "H      4.5490076673   -2.9646774751   -0.2746417121                              \n"
                              "H      3.1446405179   -1.9810706859    1.5383502882                              \n"
                              "H      6.4506580916    0.8156713166   -0.7880847806                              \n"
                              "H      5.0034051403    1.7672286916    1.0725996322                              \n"
                              "H     -0.0769625058    0.9467302475   -2.0078701638                              \n"
                              "H      2.2878857021    0.1398504722   -1.5082722535                              \n"
                              "H      1.7375985743    1.2376935514    2.7051319417                              \n"
                              "H     -0.6059179497    2.0748344376    2.2229078545 ";

  auto initial_stream = std::istringstream(initial);
  const auto molecule = irc::io::load_xyz<vec3>(initial_stream);

  auto ircs = irc::IRC<vec3, vec, mat>(molecule);

  const auto dq = vec{-0.00292000000000003,-0.00216000000000038,-0.00181000000000031,-0.00205000000000011,-0.00191000000000008,-0.0026900000000003,-0.0025599999999999,-0.00316000000000027,-0.00038999999999989,-0.00287999999999977,-0.00274999999999981,-0.00337999999999994,-0.000129999999999963,0.000319999999999876,-0.00124000000000013,-0.00117000000000012,-0.00112000000000023,-0.0010699999999999,-0.00106000000000028,-0.000979999999999759,-0.00119999999999987,-0.000939999999999941,-0.00116000000000005,0.000770000000000159,0.000449999999999839,0.000480000000000036,-0.00143000000000004,-0.00109000000000004,0.000820000000000043,0.00105000000000022,5.00000000003276E-05,-0.000039999999999818,0.000789999999999846,-0.000039999999999818,0.000539999999999985,-0.00059000000000009,0.000969999999999693,-0.000609999999999999,0.000140000000000029,0.00104000000000015,0.000299999999999745,0.000909999999999744,0.000510000000000232,-0.000500000000000167,0.000360000000000138,-0.000509999999999788,-1.00000000000655E-05,-0.000119999999999898,0.000159999999999716,-0.000680000000000014,-0.000129999999999963,-0.000300000000000189,-0.000140000000000029,-0.000599999999999934,-0.000169999999999781,-0.000430000000000152,-0.000040000000000262,-0.000379999999999825,-0.000450000000000284,-4.03000000000001E-05,-8.29999999999928E-06,-3.76999999999982E-05,-1.29999999999922E-06,7.03000000000023E-05,3.55999999999967E-05,9.00000000001455E-05,-0.000059999999999949,0.000190000000000135,0.000209999999999821,-0.00156999999999985,0.0003709,-0.0005313,0.00139,0.000878999999999963,-6.60000000000105E-05,-0.000220000000000109,-0.000199999999999978,-0.000690000000000079,0.000230000000000063,0.00158000000000014,-0.0016799999999999,-1.99999999996869E-05,5.00000000003276E-05,-0.000659999999999883,0.000359999999999694,-0.00154999999999994,0.000860000000000083,0.000350000000000072,-0.000020000000000131,-0.0001467,0.00055000000000005,-0.000239999999999796,0.00141000000000013,-0.00123999999999991,0.000020000000000131,-0.000359999999999694,0.000198500000000001,-3.00000000001965E-05,0,-0.0000222,4.99999999998835E-05,-3.00000000001965E-05,-7.20999999999968E-06,4.89000000000045E-05,0,-4.99999999998835E-05,-1.74000000000007E-05,1.00000000000655E-05,0.000020000000000131,3.97999999999996E-05,-3.16999999999991E-05,0.0000149,-9.13000000000025E-05,4.21999999999922E-06,-1.29900000000001E-05,-1.41999999999998E-05,3.79999999999998E-05,0.000276999999999999,0.00015068,0.000118190000000001,-0.0004575};

  const auto x_c = irc::molecule::to_cartesian<vec3, vec>(molecule);
  const auto q = ircs.cartesian_to_irc(x_c);

  const auto itc_result = transformation::irc_to_cartesian<vec3, vec, mat>(
      q,
      dq,
      x_c,
      ircs.get_bonds(),
      ircs.get_angles(),
      ircs.get_dihedrals(),
      ircs.get_linear_angles(),
      ircs.get_out_of_plane_bends()
      );
  CAPTURE(ircs.get_bonds().size());
  CAPTURE(ircs.get_angles().size());
  CAPTURE(ircs.get_dihedrals().size());
  CAPTURE(ircs.get_linear_angles().size());
  CAPTURE(ircs.get_out_of_plane_bends().size());
  CHECK(itc_result.converged);
}