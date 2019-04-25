#include "catch.hpp"

#include "libirc/connectivity.h"

#include "config.h"
#include "libirc/atom.h"
#include "libirc/conversion.h"
#include "libirc/io.h"
#include "libirc/molecule.h"
#include "libirc/transformation.h"

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

// TODO (Peter) Extend for linear angles
template<typename Vector3, typename Vector, typename Matrix>
std::tuple<std::vector<connectivity::Bond>,
           std::vector<connectivity::Angle>,
           std::vector<connectivity::Dihedral>,
           std::vector<connectivity::LinearAngle<Vector3>>>
badla_from_molecule(const molecule::Molecule<Vector3>& mol) {

  using namespace connectivity;

  // Compute interatomic distance for formaldehyde molecule
  Matrix dd{distances<Vector3, Matrix>(mol)};

  // Build graph based on the adjacency matrix
  UGraph adj{adjacency_matrix(dd, mol)};

  // Compute distance matrix and predecessor matrix
  Matrix dist{distance_matrix<Matrix>(adj)};

  // Compute bonds
  std::vector<Bond> B{bonds(dist, mol)};

  // Compute angles
  std::vector<Angle> A{angles(dist, mol)};

  // Compute dihedral angles
  std::vector<Dihedral> D{dihedrals(dist, mol)};

  // Compute linear angles
  std::vector<LinearAngle<Vector3>> LA{linear_angles(dist, mol)};

  // Return bonds, angles and dihedral angles
  return std::make_tuple(B, A, D, LA);
}

TEST_CASE("Bonds, angles and dihedral angles") {
  using namespace connectivity;

  CHECK_THROWS(Bond{0, 0});

  Bond b1{0, 1};
  Bond b2{1, 0};
  Bond b3{b1};
  Bond b4{0, 2};

  CHECK(b1 == b2);
  CHECK(b1 == b3);
  CHECK(b1 != b4);

  CHECK_THROWS(Angle{0, 0, 1});
  CHECK_THROWS(Angle{0, 1, 0});
  CHECK_THROWS(Angle{1, 0, 0});
  CHECK_THROWS(Angle{0, 0, 0});

  Angle a1{0, 1, 2};
  Angle a2{2, 1, 0};
  Angle a3{a1};
  Angle a4{1, 2, 3};

  CHECK(a1 == a2);
  CHECK(a1 == a3);
  CHECK(a1 != a4);

  CHECK_THROWS(Dihedral{0, 0, 1, 2});
  CHECK_THROWS(Dihedral{0, 1, 0, 2});
  CHECK_THROWS(Dihedral{0, 1, 2, 0});
  CHECK_THROWS(Dihedral{1, 0, 0, 2});
  CHECK_THROWS(Dihedral{1, 0, 2, 0});
  CHECK_THROWS(Dihedral{1, 2, 0, 0});

  CHECK_THROWS(Dihedral{1, 0, 0, 0});
  CHECK_THROWS(Dihedral{0, 1, 0, 0});
  CHECK_THROWS(Dihedral{0, 0, 1, 0});
  CHECK_THROWS(Dihedral{0, 0, 0, 1});

  CHECK_THROWS(Dihedral{0, 0, 1, 1});
  CHECK_THROWS(Dihedral{0, 1, 0, 1});
  CHECK_THROWS(Dihedral{0, 1, 1, 0});
  CHECK_THROWS(Dihedral{1, 0, 1, 0});
  CHECK_THROWS(Dihedral{1, 1, 0, 0});

  CHECK_THROWS(Dihedral{0, 0, 0, 0});

  Dihedral d1{0, 1, 2, 3};
  Dihedral d2{3, 2, 1, 0};
  Dihedral d3{d1};
  Dihedral d4{1, 2, 3, 4};

  CHECK(d1 == d2);
  CHECK(d1 == d3);
  CHECK(d1 != d4);
}

TEST_CASE("Distance, angle and dihedral angle") {
  using namespace connectivity;

  // Define four points in 3D space
  const auto p1 = vec3{0.00, 0.00, -0.25};
  const auto p2 = vec3{0.00, 0.00, 1.50};
  const auto p3 = vec3{0.00, 1.00, 1.50};
  const auto p4 = vec3{5.00, 1.00, 1.50};

  CHECK(distance(p1, p2) == Approx(1.75));

  CHECK(angle(p1, p2, p3) == Approx(tools::constants::pi / 2.));

  CHECK(dihedral(p1, p2, p3, p4) == Approx(-tools::constants::pi / 2.));
}

// Regular bond
TEST_CASE("Connectivity for compressed H2") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  double d{0.5};

  // Define compressed H2 molecule
  Molecule<vec3> molecule{{"H", {0.0, 0.0, d}}, {"H", {0.0, 0.0, 0.0}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 1);

  // Check number of angles
  CHECK(A.empty());

  // Check number of dihedral angles
  CHECK(D.empty());

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, {}, {}, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 1);

  // Check bond length
  CHECK(q(0) == Approx(d * angstrom_to_bohr));
}

// Interfragment bond
TEST_CASE("Connectivity for stretched H2") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  double d{2.5};

  // Define compressed H2 molecule
  Molecule<vec3> molecule{{"H", {0.0, 0.0, d}}, {"H", {0.0, 0.0, 0.0}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 1);

  // Check number of angles
  CHECK(A.empty());

  // Check number of dihedral angles
  CHECK(D.empty());

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 1);

  // Check bond length
  CHECK(q(0) == Approx(d * angstrom_to_bohr));
}

// Angle
TEST_CASE("Connectivity for compressed H2O") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  double d1{0.8};
  double d2{0.6};
  double angle{134};

  double angle_rad{angle * deg_to_rad};
  double a{(180. - angle) / 2.};
  double a_rad{a * deg_to_rad};

  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"O", {0, 0, 0}},
      {"H", {d1 * std::cos(a_rad), -d1 * std::sin(a_rad), 0}},
      {"H", {-d2 * std::cos(a_rad), -d2 * std::sin(a_rad), 0}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 2);

  // Check number of angles
  CHECK(A.size() == 1);

  // Check number of dihedral angles
  CHECK(D.empty());

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 3);

  // Check bond 1
  CHECK(q(0) == Approx(d1 * angstrom_to_bohr));

  // Check bond 2
  CHECK(q(1) == Approx(d2 * angstrom_to_bohr));

  // Check angle
  CHECK(q(2) == Approx(angle_rad));
}

// Angle from interfragment bonds
TEST_CASE("Connectivity for stretched H2O") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  double d1{1.3};
  double d2{1.8};
  double angle{85};

  double angle_rad{angle * deg_to_rad};
  double a{(180. - angle) / 2.};
  double a_rad{a * deg_to_rad};

  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"O", {0, 0, 0}},
      {"H", {d1 * std::cos(a_rad), -d1 * std::sin(a_rad), 0}},
      {"H", {-d2 * std::cos(a_rad), -d2 * std::sin(a_rad), 0}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 2);

  // Check number of angles
  CHECK(A.size() == 1);

  // Check number of dihedral angles
  CHECK(D.empty());

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 3);

  // Check bond 1
  CHECK(q(0) == Approx(d1 * angstrom_to_bohr));

  // Check bond 2
  CHECK(q(1) == Approx(d2 * angstrom_to_bohr));

  // Check angle
  CHECK(q(2) == Approx(angle_rad));
}

// Dihedral
TEST_CASE("Connectivity for compressed H2O2") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  // Define H2O2 molecule
  Molecule<vec3> molecule{{"O", {0.000000, 0.734058, -0.052750}},
                          {"O", {0.000000, -0.734058, -0.052750}},
                          {"H", {0.839547, 0.880752, 0.422001}},
                          {"H", {-0.839547, -0.880752, 0.422001}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 3);

  // Check number of angles
  CHECK(A.size() == 2);

  // Check number of dihedral angles
  CHECK(D.size() == 1);

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 6);

  // Atomic positions
  vec3 p1{molecule[0].position};
  vec3 p2{molecule[1].position};
  vec3 p3{molecule[2].position};
  vec3 p4{molecule[3].position};

  // Check bond 1
  CHECK(q(0) == Approx(distance(p1, p2)));

  // Check bond 2
  CHECK(q(1) == Approx(distance(p1, p3)));

  // Check bond 3
  CHECK(q(2) == Approx(distance(p2, p4)));

  // Check angle 1
  CHECK(q(3) == Approx(angle(p2, p1, p3)));

  // Check angle 2
  CHECK(q(4) == Approx(angle(p1, p2, p4)));

  // Check dihedral
  CHECK(q(5) == Approx(dihedral(p3, p1, p2, p4)));
}

// Dihedral from interfragment bonds
TEST_CASE("Connectivity for stretched H2O2") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  // Define H2O2 molecule
  Molecule<vec3> molecule{{"O", {0.000000, 0.734058, -0.052750}},
                          {"O", {0.000000, -0.734058, -0.052750}},
                          {"H", {1.839547, 1.880752, 0.422001}},
                          {"H", {-1.839547, -1.880752, 0.422001}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 3);

  // Check number of angles
  CHECK(A.size() == 2);

  // Check number of dihedral angles
  CHECK(D.size() == 1);

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 6);

  // Atomic positions
  vec3 p1{molecule[0].position};
  vec3 p2{molecule[1].position};
  vec3 p3{molecule[2].position};
  vec3 p4{molecule[3].position};

  // Check bond 1
  CHECK(q(0) == Approx(distance(p1, p2)));

  // Check bond 2
  CHECK(q(1) == Approx(distance(p1, p3)));

  // Check bond 3
  CHECK(q(2) == Approx(distance(p2, p4)));

  // Check angle 1
  CHECK(q(3) == Approx(angle(p2, p1, p3)));

  // Check angle 2
  CHECK(q(4) == Approx(angle(p1, p2, p4)));

  // Check dihedral
  CHECK(q(5) == Approx(dihedral(p3, p1, p2, p4)));
}

// Recursive interfragment bond
TEST_CASE("Connectivity for two stretched H2 molecules") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  // Define H2O2 molecule
  Molecule<vec3> molecule{{"H", {0.0, 0.0, 0.0}},
                          {"H", {2.0, 0.0, 0.0}},
                          {"H", {5.0, 0.0, 1.0}},
                          {"H", {7.0, 0.0, 1.0}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 3);

  // Check number of angles
  CHECK(A.size() == 2);

  // Check number of dihedral angles
  CHECK(D.size() == 1);

  // Compute IRC
  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 6);

  // Atomic positions
  vec3 p1{molecule[0].position};
  vec3 p2{molecule[1].position};
  vec3 p3{molecule[2].position};
  vec3 p4{molecule[3].position};

  // Check bond 1
  CHECK(q(0) == Approx(distance(p1, p2)));

  // Check bond 2
  CHECK(q(1) == Approx(distance(p2, p3)));

  // Check bond 3
  CHECK(q(2) == Approx(distance(p3, p4)));

  // Check angle 1
  CHECK(q(3) == Approx(angle(p1, p2, p3)));

  // Check angle 2
  CHECK(q(4) == Approx(angle(p2, p3, p4)));

  // Check dihedral
  CHECK(q(5) == Approx(dihedral(p1, p2, p3, p4)));
}

// Hydrogen bond (without quasi-linear angles)
TEST_CASE("Connectivity for bent water dimer") {
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  double d1{1.3};
  double d2{1.4};
  double a_deg{102.03};

  double a{(180. - a_deg) * deg_to_rad};

  double cos_a{std::cos(a)};
  double sin_a{std::sin(a)};

  // Define water dimer
  Molecule<vec3> molecule{{"O", {-1.464, 0.099, -0.300}},
                          {"H", {-1.956, 0.624, -0.340}},
                          {"H", {-1.797, -0.799, 0.206}},
                          {"O", {1.369, 0.146, -0.395}},
                          {"H", {1.894, 0.486, 0.335}},
                          {"H", {0.451, 0.165, -0.083}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 5);

  // Check number of angles
  CHECK(A.size() == 5);

  // TODO: Check wit other codes (where dihedral 2-1-3-6 is added).
  // Check number of dihedral angles
  CHECK(D.size() == 3);

  // Compute IRC
  vec q{cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, {})};

  // Check number of IRC
  CHECK(linalg::size<vec>(q) == 13);

  // Atomic positions
  vec3 p1{molecule[0].position};
  vec3 p2{molecule[1].position};
  vec3 p3{molecule[2].position};
  vec3 p4{molecule[3].position};
  vec3 p5{molecule[4].position};
  vec3 p6{molecule[5].position};

  // Check bond 1
  CHECK(q(0) == Approx(distance(p1, p2)));

  // Check bond 2
  CHECK(q(1) == Approx(distance(p1, p3)));

  // Check bond 3
  CHECK(q(2) == Approx(distance(p4, p5)));

  // Check bond 4
  CHECK(q(3) == Approx(distance(p1, p6)));

  // Check bond 5
  CHECK(q(4) == Approx(distance(p4, p6)));

  // Check angle 1
  CHECK(q(5) == Approx(angle(p2, p1, p3)));

  // Check angle 2
  CHECK(q(6) == Approx(angle(p1, p6, p4)));

  // Check angle 3
  CHECK(q(7) == Approx(angle(p2, p1, p6)));

  // Check angle 4
  CHECK(q(8) == Approx(angle(p3, p1, p6)));

  // Check angle 5
  CHECK(q(9) == Approx(angle(p5, p4, p6)));

  // Check dihedral 1
  CHECK(q(10) == Approx(dihedral(p2, p1, p6, p4)));

  // Check dihedral 2
  CHECK(q(11) == Approx(dihedral(p3, p1, p6, p4)));

  // Check dihedral 3
  CHECK(q(12) == Approx(dihedral(p1, p6, p4, p5)));
}

TEST_CASE("Linear angles") {
  using irc::connectivity::angle;

  SECTION("orthogonal to X and Y") {
    const vec3 p1 = {0, 0, 0};
    const vec3 p2 = {1, 1, 1};
    const vec3 p3 = {2, 2, 2};

    CHECK(angle(p1, p2, p3) == tools::constants::pi);
  };

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  // Define H2O2 molecule
  Molecule<vec3> molecule{{"C", {0.0, 1.0, 2.0}},
                          {"C", {2.0, 1.0, 2.0}},
                          {"C", {4.0, 1.0, 2.0}},
                          {"C", {6.0, 1.0, 2.0}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::vector<LinearAngle<vec3>> LA;
  std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(molecule);

  CHECK(B.size() == 3);
  CHECK(A.size() == 0);
  // TODO(Peter) set to == 0 when dihedral bug is fixed
  CHECK(D.size() == 1);
  CHECK(LA.size() == 2 * 2);

  vec q{connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(molecule), B, A, D, LA)};

  // TODO(Peter) set to == 7 when dihedral bug is fixed
  REQUIRE(linalg::size<vec>(q) == 8);

  // TODO(Peter) Check element q(3) rather than q(7) when dihedral bug is fixed
  CHECK(q(4) == Approx(tools::constants::pi));
  CHECK(q(5) == Approx(tools::constants::pi));
  CHECK(q(6) == Approx(tools::constants::pi));
  CHECK(q(7) == Approx(tools::constants::pi));
}

// Quasi-linear angles
TEST_CASE("Connectivity for water dimer") {
  // TODO
}

// Interfragment bonds
TEST_CASE("Connectivity for benzene dimer") {
  // TODO
}

TEST_CASE("Connectivity test for CH2O") {
  // TODO
  /*
  using namespace std;

  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  double d1{1.3};
  double d2{1.4};
  double a_deg{102.03};

  double a{(180. - a_deg) * deg_to_rad};

  double cos_a{std::cos(a)};
  double sin_a{std::sin(a)};

  // Define water dimer
  Molecule<vec3> molecule{{"H", { 1.0686, -0.1411,  1.0408}},
                          {"C", { 0.5979,  0.0151,  0.0688}},
                          {"H", { 1.2687,  0.2002, -0.7717}},
                          {"O", {-0.5960, -0.0151, -0.0686}}};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);

  // Get bonds, angles and dihedrals
  std::vector<Bond> B;
  std::vector<Angle> A;
  std::vector<Dihedral> D;
  std::tie(B, A, D) = bad_from_molecule<vec3, vec, mat>(molecule);

  // Check number of bonds
  CHECK(B.size() == 3);

  // Check number of angles
  CHECK(A.size() == 3);

  // Check number of dihedral angles
  CHECK(D.size() == 1);

  // Compute IRC
  vec q{
      cartesian_to_irc<vec3, vec>(to_cartesian<vec3, vec>(molecule), B, A, D)};

  // Check number of IRC
  REQUIRE(linalg::size<vec>(q) == 7);
   */
}

TEST_CASE("Connectivity test for 2-butyne") {
  // TODO
}

TEST_CASE("Connectivity of molecule database") {
  using namespace io;

  using namespace connectivity;
  using namespace molecule;
  using namespace tools;

  struct ConnectivityResult {
    std::string filename;
    std::size_t n_bonds;
    std::size_t n_angles;
    std::size_t n_dihedrals;
    std::size_t n_linear_angles;
  };

  auto results =
      std::vector<ConnectivityResult>{{"carbon_dioxide.xyz", 2, 0, 0, 2},
                                      {"hydrogen_peroxide.xyz", 3, 2, 1, 0},
                                      {"ethanol.xyz", 8, 13, 12, 0},
                                      {"glycerol.xyz", 13, 21, 27, 0},
                                      {"octane.xyz", 25, 48, 63, 0},
                                      {"phenol.xyz", 13, 19, 26, 0},
                                      {"indene.xyz", 18, 30, 44, 0},
                                      {"toluene.xyz", 15, 24, 30, 0},
                                      {"caffeine.xyz", 25, 43, 54, 0}};

  for (const auto& molecule_parameters : results) {
    CAPTURE(molecule_parameters.filename);

    const auto mol =
        load_xyz<vec3>(config::molecules_dir + molecule_parameters.filename);

    std::vector<Bond> B;
    std::vector<Angle> A;
    std::vector<Dihedral> D;
    std::vector<LinearAngle<vec3>> LA;
    std::tie(B, A, D, LA) = badla_from_molecule<vec3, vec, mat>(mol);

    CHECK(B.size() == molecule_parameters.n_bonds);
    CHECK(A.size() == molecule_parameters.n_angles);
    CHECK(D.size() == molecule_parameters.n_dihedrals);
    CHECK(LA.size() == molecule_parameters.n_linear_angles);
  }
}
