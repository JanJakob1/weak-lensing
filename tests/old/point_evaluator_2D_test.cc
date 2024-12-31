// Copyright (C) 2011-2021 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for
// more details.
//
// You should have received a copy of the European Union Public Licence (EUPL)
// v1.2 along with HiFlow3.  If not, see
// <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

/// \author Jonathan Schwegler

#define BOOST_TEST_MODULE point_evaluator2D

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

static const char *PARAM_FILENAME = "point_eval_test.xml";
#ifndef MESH_DATADIR
#define MESH_DATADIR "./"
#endif
static const char *DATADIR = MESH_DATADIR;

// System includes.
#include "hiflow.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

// All names are imported for simplicity.
using namespace hiflow;
using namespace hiflow::doffem;
using namespace hiflow::la;
using namespace hiflow::mesh;

// Shorten some datatypes with typedefs.
typedef LADescriptorCoupledD LAD;
typedef LAD::DataType Scalar;
typedef LAD::VectorType Vector;
typedef LAD::MatrixType Matrix;

// Rank of the master process.
const int MASTER_RANK = 0;

// Dimension of the problem.
const int DIMENSION = 2;

double exact_sol(std::vector< double > pt) {
  return std::accumulate(pt.begin(), pt.end(), 0.);
}

class LocalMassMatrixAssembler
  : private AssemblyAssistant< DIMENSION, double > {
public:
  void operator()(const Element< double, DIMENSION > &element,
                  const Quadrature< double > &quadrature, LocalMatrix &lm) {
    AssemblyAssistant< DIMENSION, double >::initialize_for_element(element,
        quadrature);
    // compute local matrix
    const int num_q = num_quadrature_points();
    for (int q = 0; q < num_q; ++q) {
      const double wq = w(q);
      const int n_dofs = num_dofs(0);
      for (int i = 0; i < n_dofs; ++i) {
        for (int j = 0; j < n_dofs; ++j) {
          lm(dof_index(i, 0), dof_index(j, 0)) +=
            wq * phi(j, q) * phi(i, q) * std::abs(detJ(q));
        }
      }
    }
  }

  void operator()(const Element< double, DIMENSION > &element,
                  const Quadrature< double > &quadrature, LocalVector &lv) {
    AssemblyAssistant< DIMENSION, double >::initialize_for_element(element,
        quadrature);
    // compute local vector
    const int num_q = num_quadrature_points();
    for (int q = 0; q < num_q; ++q) {
      const double wq = w(q);
      const int n_dofs = num_dofs(0);
      for (int i = 0; i < n_dofs; ++i) {
        lv[dof_index(i, 0)] += wq * f(x(q)) * phi(i, q) * std::abs(detJ(q));
      }
    }
  }

  double f(Vec< DIMENSION, double > pt) {
    double rhs_sol;
    std::vector< double > converted_point(DIMENSION);
    for (int i = 0; i < DIMENSION; ++i) {
      converted_point[i] = pt[i];
    }
    rhs_sol = exact_sol(converted_point);

    return rhs_sol;
  }
};

class PointEvalTest {
public:
  PointEvalTest(const std::string param_filename)
    : comm_(MPI_COMM_WORLD), rank_(-1), num_partitions_(-1),
      params_(DATADIR + param_filename, MASTER_RANK, MPI_COMM_WORLD), rhs_(0),
      sol_(0), matrix_(0) {
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &num_partitions_);
  }

  // Main algorithm

  void run() {
    // Construct / read in the initial mesh.
    build_initial_mesh();
    // Initialize space and linear algebra.
    prepare_system();
    // Compute the stiffness matrix and right-hand side.
    assemble_system();
    // Solve the linear system.
    solve_system();
    // Compute the error to the exact solution.
    compute_error();
  }

  ~PointEvalTest() {
    delete matrix_;
    delete sol_;
    delete rhs_;
  }

private:
  std::string mesh_filename_;
  void build_initial_mesh();
  // Setup space, linear algebra, and compute Dirichlet values.
  void prepare_system();
  // Compute the matrix and rhs.
  void assemble_system();
  // Compute solution x.
  void solve_system();
  // Compute errors compared to exact solution.
  void compute_error();
  // MPI communicator.
  MPI_Comm comm_;
  // Local process rank and number of processes.
  int rank_, num_partitions_;

  // Local mesh and mesh on master process.
  MeshPtr mesh_, master_mesh_;
  // Solution space.
  VectorSpace< double, DIMENSION > space_;
  // Parameter data read in from file.
  PropertyTree params_;
  // Linear algebra couplings helper object.
  LaCouplings couplings_;
  // Vectors for solution and load vector.
  CoupledVector< Scalar > *rhs_, *sol_;
  // System matrix.
  CoupledMatrix< Scalar > *matrix_;

  // Global assembler.
  StandardGlobalAssembler< double, DIMENSION > global_asm_;
};

BOOST_AUTO_TEST_CASE(point_evaluator2D, *utf::tolerance(1.0e-10)) {

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;

  MPI_Init(&argc, &argv);

  // set default parameter file
  std::string param_filename(PARAM_FILENAME);

  PointEvalTest tester(param_filename);
  tester.run();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}

void PointEvalTest::build_initial_mesh() {
  // Read in the mesh on the master process. The mesh is chosen according to the
  // dimension of the problem.
  if (rank_ == MASTER_RANK) {
    std::string mesh_name;
    switch (DIMENSION) {
    case (2): {
      mesh_name = params_["Mesh"]["Filename2"].get< std::string >(
                    "unit_square_tri_quad.inp");
    }
    break;
    case (3): {
      mesh_name = params_["Mesh"]["Filename3"].get< std::string >(
                    "unit_cube_tetras_3d.inp");
    }
    break;
    }

    std::string mesh_filename;
    mesh_filename = std::string(DATADIR) + mesh_name;

    master_mesh_ = read_mesh_from_file(mesh_filename, DIMENSION, DIMENSION, 0);

    // Refine the mesh until the initial refinement level is reached.
    const int initial_ref_lvl = params_["Mesh"]["RefLevel"].get< int >(2);
    for (int r = 0; r < initial_ref_lvl; ++r) {
      master_mesh_ = master_mesh_->refine();
    }
  }

  // Distribute mesh over all processes, and compute ghost cells
  int uniform_ref_steps;
  MeshPtr local_mesh = partition_and_distribute(master_mesh_, MASTER_RANK,
                       comm_, uniform_ref_steps);
  SharedVertexTable shared_verts;
  mesh_ = compute_ghost_cells(*local_mesh, comm_, shared_verts);
}

void PointEvalTest::prepare_system() {
  // Assign degrees to each element.
  const int fe_degree = params_["Mesh"]["FeDegree"].get< int >(1);
  std::vector< int > degrees(1, fe_degree);
  std::vector< bool > is_cg(1, true);
  std::vector< FEType > fe_ansatz (1, FE_TYPE_LAGRANGE);
  
  // Initialize the VectorSpace object.
  space_.Init(*mesh_, fe_ansatz, is_cg, degrees, CUTHILL_MCKEE);

  // Setup couplings object.
  couplings_.Init(comm_);

  std::vector< int > global_offsets;
  std::vector< int > ghost_dofs;
  std::vector< int > ghost_offsets;
  space_.get_la_couplings(global_offsets, ghost_dofs, ghost_offsets);

  couplings_.InitializeCouplings(global_offsets, ghost_dofs, ghost_offsets);

  // Compute the matrix graph.
  SparsityStructure sparsity;
  global_asm_.compute_sparsity_structure(space_, sparsity);

  // Setup linear algebra objects.
  delete matrix_;
  delete rhs_;
  delete sol_;

  CoupledMatrixFactory< Scalar > CoupMaFact;
  matrix_ = CoupMaFact
            .Get(params_["LinearAlgebra"]["NameMatrix"].get< std::string >(
                   "CoupledMatrix"))
            ->params(params_["LinearAlgebra"]);
  matrix_->Init(comm_, couplings_);
  CoupledVectorFactory< Scalar > CoupVecFact;
  rhs_ = CoupVecFact
         .Get(params_["LinearAlgebra"]["NameVector"].get< std::string >(
                "CoupledVector"))
         ->params(params_["LinearAlgebra"]);
  sol_ = CoupVecFact
         .Get(params_["LinearAlgebra"]["NameVector"].get< std::string >(
                "CoupledVector"))
         ->params(params_["LinearAlgebra"]);
  rhs_->Init(comm_, couplings_);
  sol_->Init(comm_, couplings_);

  // Initialize structure of LA objects.
  matrix_->InitStructure(
    vec2ptr(sparsity.diagonal_rows), vec2ptr(sparsity.diagonal_cols),
    sparsity.diagonal_rows.size(), vec2ptr(sparsity.off_diagonal_rows),
    vec2ptr(sparsity.off_diagonal_cols), sparsity.off_diagonal_rows.size());

  rhs_->InitStructure();
  sol_->InitStructure();

  // Zero all linear algebra objects.
  matrix_->Zeros();
  rhs_->Zeros();
  sol_->Zeros();
}

void PointEvalTest::assemble_system() {
  // Assemble matrix and right-hand-side vector.
  LocalMassMatrixAssembler local_asm;
  global_asm_.assemble_matrix(space_, local_asm, *matrix_);
  global_asm_.assemble_vector(space_, local_asm, *rhs_);
}

void PointEvalTest::solve_system() {
  LinearSolver< LAD > *solver_;
  LinearSolverFactory< LAD > SolFact;
  solver_ =
    SolFact.Get(params_["LinearSolver"]["Name"].get< std::string >("CG"))
    ->params(params_["LinearSolver"]);
  solver_->SetupOperator(*matrix_);
  solver_->Solve(*rhs_, sol_);

  delete solver_;
}

void PointEvalTest::compute_error() {
  sol_->Update();

  srand(time(NULL));
  int num_points = params_["Test"]["NumPts"].get< int >(1000);
  bool with_comm = params_["Test"]["WithCommunication"].get< bool >(false);

  // measure the initialization time of the PointEvaluator
  double timer_init = MPI_Wtime();

  // set up the PointEvaluator
<<<<<<< HEAD
  GridGeometricSearch< double, DIMENSION > geom_search(mesh_);
  PointEvaluator< double, DIMENSION > evaluator(space_);
  double timer_eval = MPI_Wtime();
  timer_init = timer_eval - timer_init;

  EvalFeFunction< LAD, DIMENSION > fun(space_, *sol_);
=======
  GridGeometricSearch geom_search(mesh_);
  PointEvaluator< double, double > evaluator(space_);
  double timer_eval = MPI_Wtime();
  timer_init = timer_eval - timer_init;

  FeEvalOnCellScalar< LAD > fun(space_, *sol_);
>>>>>>> master_19_01_22_2c0f

  double max_diff = 0;
  for (int i = 0; i < num_points; ++i) {
    // get some "random" test points.
    Vec<DIMENSION, double > pt;
    for (int j = 0; j < DIMENSION; ++j) {
      pt[j] = ((double)rand() / (RAND_MAX));
    }
    double value;

    bool has_point;
    if (with_comm) {
      has_point = evaluator.evaluate_fun_global(fun, &geom_search, pt, value, comm_);
    } else {
      has_point = evaluator.evaluate_fun(fun, &geom_search, pt, value);
    }

    if (has_point) {
      double diff = std::abs(exact_sol(pt) - value);
      if (max_diff < diff) {
        max_diff = diff;
      }
    }
  }

  timer_eval = MPI_Wtime() - timer_eval;

  // get the time needed for just finding the points for comparison
  double timer_finder = MPI_Wtime();
  for (int i = 0; i < num_points; ++i) {
    Vec<DIMENSION, double > pt;
    for (int j = 0; j < DIMENSION; ++j) {
      pt[j] = ((double)rand() / (RAND_MAX));
    }
    std::vector< double > ref_coords;
    std::vector< int > cells;
    //std::vector< std::vector< Coordinate > > ref_points;
    geom_search.find_cell(pt, cells/*, ref_points*/);
  }
  timer_finder = MPI_Wtime() - timer_finder;

  std::cout << "Time for initializing GridGeometricSearch: " << timer_init
            << std::endl;
  std::cout << "Time for evaluating " << num_points << " points: " << timer_eval
            << std::endl;
  std::cout << "Time for finding " << num_points << " points: " << timer_finder
            << std::endl;
  std::cout << "Mean Cells per GridCell: " << geom_search.mean_mesh_cells()
            << std::endl;
  std::cout << "Max Cells per GridCell: " << geom_search.max_mesh_cells()
            << std::endl;

  std::cout << "Max Diff (~L_inf) on process " << rank_ << ":  " << max_diff
            << std::endl;

  BOOST_TEST(max_diff == 0.0);
}
