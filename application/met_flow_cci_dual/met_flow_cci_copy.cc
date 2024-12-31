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

#include "met_flow_cci.h"

#include <fstream>

MetFlowCCIApp::MetFlowCCIApp(const std::string &param_filename)
    : MetFlowDualApp(param_filename) {
  MetFlowApp::local_asm_primal_ = &my_local_asm_;

  is_pod_active_ = true;
  is_dual_pod_active_ = true;
}

MetFlowCCIApp::~MetFlowCCIApp() {}

// read periodic or non-periodic mesh.

void MetFlowCCIApp::read() {
  MetFlowApp::prepare_periodicity();
MetFlowApp:
  read_periodic_mesh();
}

// prepare application

void MetFlowCCIApp::prepare() { MetFlowApp::prepare(); }

// prepare application

void MetFlowCCIApp::prepare_dual() {
  MetFlowDualApp::prepare_dual();

  prepare_goal_functional();
}

void MetFlowCCIApp::prepare_time() {
  MetFlowApp::prepare_time();

  bool success = false;
  success = my_local_asm_.set_timestepping_method(method_);
  interminable_assert(success);
  // my_local_asm_.set_time_discretization_to_crank_nicolson();

  my_local_asm_.set_dT(delta_t_);
}

void MetFlowCCIApp::prepare_goal_functional() {
  // choose which simulation parameters
  int ref_level = params_["Mesh"]["InitialRefLevel"].get< int >();
  double distance = params_["CCI"]["Distance"].get< double >();

  bool L4 = false;
  if (ref_level == 4)
    L4 = true;

  bool diverging = false;
  if (distance == 400.)
    diverging = true;

  if (goal_contribution_ == "Vorticity") {
    goal_functional_ = new GoalFunctionalVorticity< DIM >();
    goal_functional_->force_type_active() = false;
    goal_functional_->initial_type_active() = true;
    double max_val = 0.;

    if (L4 && diverging)
      max_val = 0.000745097;
    if (!L4 && diverging)
      max_val = 0.000723312;
    if (L4 && !diverging)
      max_val = 0.00125267;
    if (!L4 && !diverging)
      max_val = 0.000865425;

    reinterpret_cast< GoalFunctionalVorticity< DIM > * >(goal_functional_)
        ->lower_bound() = max_val * 0.5;
    // std::cout << "Goal-Functional: " << goal_functional_->get_description()
    // << std::endl;
  }

  if (goal_contribution_ == "KineticEnergy") {
    goal_functional_ = new GoalFunctionalKineticEnergy< DIM >();
    goal_functional_->force_type_active() = false;
    goal_functional_->initial_type_active() = true;
    double max_val = 0.;

    if (L4 && diverging)
      max_val = 0.0293852;
    if (!L4 && diverging)
      max_val = 0.030158;
    if (L4 && !diverging)
      max_val = 0.0369141;
    if (!L4 && !diverging)
      max_val = 0.0375191;

    reinterpret_cast< GoalFunctionalKineticEnergy< DIM > * >(goal_functional_)
        ->bound() = max_val * 0.9;
    // std::cout << "Goal-Functional: " << goal_functional_->get_description()
    // << std::endl;
  }

  my_local_asm_.set_goal_functional(goal_functional_);
}

// non-adaptive run

void MetFlowCCIApp::run() {

  this->prepare();

  int time_step = 0;

  std::cout << "Solving instationary Navier Stokes with " << method_
            << std::endl;
  std::cout << "--->  Visualizing solution at time index " << time_step
            << " with current time: " << delta_t_ * time_step << std::endl;
  std::cout << "==============================" << std::endl << std::endl;

  MetFlowApp::visualize_solution(0);
  sol_prev_.CloneFrom(sol_);

  // time-stepping loop
  while (time_step < num_time_steps_) {
    my_local_asm_.set_prev_solution(sol_prev_);

    solve_nlp();

    compute_vorticity();

    compute_storm_track();

    sol_prev_.CloneFrom(sol_);

    ++time_step;

    if (time_step % visual_step_ == 0) {
      std::cout << "--->  Visualizing solution at time index " << time_step
                << " with current time: " << delta_t_ * time_step << std::endl;
      MetFlowApp::visualize_solution(time_step);
      std::cout << "==============================" << std::endl << std::endl;
    }
  } // time-stepping loop
}

void MetFlowCCIApp::adaptive_run() {
  int adaption_counter = params_["Mesh"]["RestartAt"].get< int >();
  num_adaptions_ = params_["Mesh"]["NumAdaptionCycles"].get< int >();
  int num_cells = params_["Mesh"]["ConstNumCells"].get< int >();
  double accept_percent = params_["Mesh"]["MeshTolerance"].get< int >();
  int accept_cells = accept_percent / 100.;

  if (adaption_counter != 0)
    restart();

  while (adaption_counter < num_adaptions_) {
    std::cout << "Number of adaptions " << num_adaptions_ << "\n";
    std::cout << "ADAPTION STEP " << adaption_counter << "\n";

    // Prepare application
    this->prepare();

    // Prepare storm track
    compute_vorticity();
    prepare_storm_track(adaption_counter);
    compute_storm_track();
    write_storm_track(0, adaption_counter);

    // Prepare error estimators
    prepare_indicators();
    estimate_indicators();
    write_error(0, adaption_counter);

    // Now the main loop ...
    int time_step = 0;

    std::cout << "Solving instationary Navier Stokes with " << method_
              << std::endl;
    std::cout << "--->  Visualizing solution at time index " << time_step
              << " with current time: " << delta_t_ * time_step << std::endl;
    std::cout << "==============================" << std::endl << std::endl;

    visualize_solution(0, adaption_counter);
    visualize_mesh(adaption_counter);

    sol_prev_.CloneFrom(sol_);

    // Time Stepping Loop
    while (time_step < num_time_steps_) {
      my_local_asm_.set_prev_solution(sol_prev_);

      success = solve_nlp();

      if (!success) {
        break;
      }

      sol_prev_.CloneFrom(sol_);

      ++time_step;

      if ((time_step % visual_step_ == 0) || (time_step % track_step_ == 0) ||
          ((time_step % indicator_step_ == 0) &&
           (refinement_criterion_ == "vorticity")))
        compute_vorticity();

      if (time_step % visual_step_ == 0) {
        std::cout << "--->  Visualizing solution at time index " << time_step
                  << " with current time: " << delta_t_ * time_step
                  << std::endl;
        visualize_solution(time_step, adaption_counter);
        std::cout << "==============================" << std::endl << std::endl;
      }

      if (time_step % indicator_step_ == 0) {
        estimate_indicators();
        write_temp_indicators(adaption_counter, time_step);
        write_error(time_step, adaption_counter);
      }
      if (time_step % track_step_ == 0) {
        compute_storm_track();
        write_storm_track(time_step, adaption_counter);
      }
    }

    write_indicators(adaption_counter);
    visualize_indicators(adaption_counter);
    write_log_file(adaption_counter);

    ++adaption_counter;

    if (!success) {
      int restart_point;
      if (time_step > 2 * indicator_step_)
        restart_point =
            time_step - (time_step % indicator_step_) - indicator_step_;
      else
        restart_point = time_step - (time_step % indicator_step_);

      std::cout << " LOADING STATE OF KILLED JOB at " << restart_point
                << "ind count " << adaption_counter << "\n";
      load_temp_indicator(adaption_counter - 1, restart_point);
    }

    if (adaption_counter <= num_adaptions_) {
      // Adapt mesh
      if (refinement_criterion_ == "ns_aposteriori") {
        // adapt_mesh_const_num_strategy(adaption_counter, ns_apost);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, ns_apost);
      } else if (refinement_criterion_ == "aposteriori") {
        // adapt_mesh_const_num_strategy(adaption_counter, apost);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, apost);
      } else if (refinement_criterion_ == "energy") {
        // adapt_mesh_const_num_strategy(adaption_counter, energy);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, energy);
      } else if (refinement_criterion_ == "vorticity") {
        // adapt_mesh_const_num_strategy(adaption_counter, vort);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, vort);
      } else {
        std::cout << "Refinement Criterion unknown!\n";
        break;
      }
    }
  }
}

void MetFlowCCIApp::dwr_run() {

  int adaption_counter = params_["Mesh"]["RestartAt"].get< int >();
  num_adaptions_ = params_["Mesh"]["NumAdaptionCycles"].get< int >();
  int num_cells = params_["Mesh"]["ConstNumCells"].get< int >();
  double accept_percent = params_["Mesh"]["MeshTolerance"].get< int >();
  int accept_cells = accept_percent / 100.;

  //************************************************************
  // POD

  // PRIMAL PROBLEM
  // For HDF5
  std::string filename_primal =
      params_["WeightedErrorEstimator"]["PrimalSnapshots"].get< std::string >();
  std::string groupname_primal =
      params_["WeightedErrorEstimator"]["PrimalSnapshotsGroup"]
          .get< std::string >();
  std::string prefix_primal =
      params_["WeightedErrorEstimator"]["PrimalSnapshotsPrefix"]
          .get< std::string >();
  std::vector< std::string > datasetname_primal;

  // setup finite element ansatz (PRIMAL PROBLEM)
  const int velocity_deg =
      params_["FiniteElements"]["VelocityDegree"].get< int >();
  const int pressure_deg =
      params_["FiniteElements"]["PressureDegree"].get< int >();
  std::vector< int > degrees(DIM + 1, velocity_deg);
  degrees[DIM] = pressure_deg;

  std::string out_filename_primal =
      params_["WeightedErrorEstimator"]["PrimalBasis"].get< std::string >();
  std::string out_groupname_primal =
      params_["WeightedErrorEstimator"]["PrimalBasisGroup"]
          .get< std::string >();
  std::string out_prefix_primal =
      params_["WeightedErrorEstimator"]["PrimalBasisPrefix"]
          .get< std::string >();
  int basis_length_primal = 0;

  // DUAL PROBLEM
  // For HDF5
  std::string filename_dual =
      params_["WeightedErrorEstimator"]["DualSnapshots"].get< std::string >();
  std::string groupname_dual =
      params_["WeightedErrorEstimator"]["DualSnapshotsGroup"]
          .get< std::string >();
  std::string prefix_dual =
      params_["WeightedErrorEstimator"]["DualSnapshotsPrefix"]
          .get< std::string >();
  std::vector< std::string > datasetname_dual;

  // setup finite element ansatz (DUAL PROBLEM)
  const int velocity_deg_dual =
      params_["FiniteElementsDual"]["VelocityDegree"].get< int >();
  const int pressure_deg_dual =
      params_["FiniteElementsDual"]["PressureDegree"].get< int >();
  std::vector< int > degrees_dual(DIM + 1, velocity_deg_dual);
  degrees_dual[DIM] = pressure_deg_dual;

  std::string out_filename_dual =
      params_["WeightedErrorEstimator"]["DualBasis"].get< std::string >();
  std::string out_groupname_dual =
      params_["WeightedErrorEstimator"]["DualBasisGroup"].get< std::string >();
  std::string out_prefix_dual =
      params_["WeightedErrorEstimator"]["DualBasisPrefix"].get< std::string >();
  int basis_length_dual = 0;
  //************************************************************

  //   if (adaption_counter != 0)
  //     restart();

  while (adaption_counter < num_adaptions_) {
    std::cout << "Number of adaptions " << num_adaptions_ << "\n";
    std::cout << "ADAPTION STEP " << adaption_counter << "\n";

    // Prepare application
    prepare();
    prepare_dual();

    int time_step = 0;

    visualize_mesh(adaption_counter);

    // Log data
    std::map< std::string, double > log_data;

    // ************************************************************
    // PRIMAL PROBLEM
    // ************************************************************

    std::cout << "***************************************************"
              << std::endl;
    std::cout << "NOW THE PRIMAL PROBLEM (AC=" << adaption_counter << ")"
              << std::endl;
    std::cout << "***************************************************"
              << std::endl;

    Timer timer;
    timer.reset();
    timer.start();

    std::string primal_mode =
        params_["WeightedErrorEstimator"]["PrimalMode"].get< std::string >();

    if (primal_mode == "Calculate") {
      std::cout << "-> PrimalMode = " << primal_mode << std::endl;

      time_step = 1;

      // initialized or resume simulation
      std::cout << "-> Initializing Simulation" << std::endl;

      if (params_["WeightedErrorEstimator"]["PrimalOffset"].get< int >()) {
        // PrimalOffset is the index of the last existing file
        time_step =
            1 + params_["WeightedErrorEstimator"]["PrimalOffset"].get< int >();
        std::cout << "-> Resuming Simulation at step " << time_step
                  << std::endl;

        // Take care for HDF5 data set names
        for (int i = 0; i < time_step; ++i) {
          std::stringstream sp;
          sp << prefix_primal << i;
          datasetname_primal.push_back(sp.str());
        }

        // read last primal solution from file
        // std::stringstream ss7;
        // ss7 << "out/primal." << time_step-1 << ".dat";
        // read_solution_from_file(sol_prev_, ss7.str());
        std::cout << "OPEN FILE " << filename_primal << ", group "
                  << groupname_primal << ", data set "
                  << datasetname_primal[time_step - 1] << std::endl;
        sol_prev_.ReadHDF5(filename_primal, groupname_primal,
                           datasetname_primal[time_step - 1]);

        sol_.CloneFrom(sol_prev_);
      } else {
        // initial conditions for primal solution
        prepare_ic();
        prepare_vorticity(sol_);

        // visualize initial solution
        std::cout << "->  Visualize solution at time step " << time_step
                  << " with current time " << delta_t_ * time_step << std::endl;
        visualize_solution(0, adaption_counter);

        // store solution for primal problem
        // std::stringstream ss;
        // ss << "out/primal." << time_step << ".dat";
        // save_solution_to_file(sol_, ss.str());

        // Write solution and take care for HDF5 data set names
        std::stringstream sp;
        sp << prefix_primal << 0;
        sol_.WriteHDF5(filename_primal, groupname_primal, sp.str());
        datasetname_primal.push_back(sp.str());
      }

      interminable_assert(time_step > 0);
      interminable_assert(time_step <= num_time_steps_);

      //       // Prepare storm track
      //       compute_vorticity();
      //       prepare_storm_track(adaption_counter);
      //       compute_storm_track();
      //       write_storm_track(0, adaption_counter);

      // Time Stepping Loop

      while (time_step <= num_time_steps_) {
        std::cout << "ntimestep: " << time_step << std::endl;
        std::cout << "ntimesteps: " << num_time_steps_ << std::endl;

        my_local_asm_.set_prev_solution(sol_prev_);

        // solve nonlinear problem
        success = solve_nlp();

        // store solution for dual problem
        // std::stringstream ss2;
        // ss2 << "out/primal." << time_step << ".dat";
        // save_solution_to_file(sol_, ss2.str());

        //         Timer my_timer;
        //         my_timer.start();

        // Write solution and take care for HDF5 data set names
        std::stringstream sp2;
        sp2 << prefix_primal << time_step;
        sol_.WriteHDF5(filename_primal, groupname_primal, sp2.str());
        datasetname_primal.push_back(sp2.str());

        //         my_timer.stop();
        //         std::cout << "Time to write to HDF5: " <<
        //         my_timer.get_duration() << std::endl;

        if (!success) {
          std::cout << "ERROR: NL problem wasn't solved!" << std::endl;
          break;
        }

        sol_prev_.CloneFrom(sol_);

        if ((time_step % visual_step_ == 0) || (time_step == num_time_steps_) ||
            (time_step % track_step_ == 0) ||
            ((time_step % indicator_step_ == 0) &&
             (refinement_criterion_ == "vorticity"))) {
          //           my_timer.reset();
          //           my_timer.start();
          prepare_vorticity(sol_);
          //           my_timer.stop();
          //           std::cout << "Time to write compute vorticity: " <<
          //           my_timer.get_duration() << std::endl;
        }

        if ((time_step % visual_step_ == 0) || (time_step == num_time_steps_)) {
          // visualize initial solution
          std::cout << "->  Visualize solution at time step " << time_step
                    << " with current time " << delta_t_ * time_step
                    << std::endl;
          //           my_timer.reset();
          //           my_timer.start();
          visualize_solution(time_step, adaption_counter);
          //           my_timer.stop();
          //           std::cout << "Time to visualize: " <<
          //           my_timer.get_duration() << std::endl;
        }

        // 	  if (time_step%indicator_step_ == 0)
        // 	    {
        // 	      estimate_indicators();
        // 	      write_temp_indicators(adaption_counter, time_step);
        // 	      write_error(time_step, adaption_counter);
        // 	    }
        // 	  if (time_step%track_step_ == 0)
        // 	    {
        // 	      compute_storm_track();
        // 	      write_storm_track(time_step, adaption_counter);
        // 	    }

        ++time_step;
      }
      --time_step;
    } else if (primal_mode == "Load") {
      std::cout << "-> PrimalMode = " << primal_mode << std::endl;
      time_step = num_time_steps_;

      // Take care for HDF5 data set names
      for (int i = 0; i <= time_step; ++i) {
        std::stringstream sp;
        sp << prefix_primal << i;
        datasetname_primal.push_back(sp.str());
      }
    } else {
      std::cout << "-> PrimalMode = " << primal_mode << std::endl;
      std::cout << "-> ERROR: " << primal_mode
                << " not a treated parameter choice!" << std::endl;
      quit_program();
    }

    timer.stop();
    log_data["Primal Problem, Time"] = timer.get_duration();

    // ************************************************************
    // POD (Primal problem)
    // ************************************************************

    std::string primal_pod_mode =
        params_["WeightedErrorEstimator"]["PrimalPODMode"].get< std::string >();

    const double tol = 1.0e-16; // kleinster EW, der nicht Null ist

    if (primal_pod_mode == "Calculate") {
      std::cout << "-> PrimalPODMode = " << primal_pod_mode << std::endl;

      std::cout << "***************************************************"
                << std::endl;
      std::cout << "NOW THE PRIMAL POD BASIS (AC=" << adaption_counter << ")"
                << std::endl;
      std::cout << "***************************************************"
                << std::endl;

      timer.reset();
      timer.start();

      std::cout << "-> Preparation" << std::endl;

      POD primal_pod(master_rank_, comm_, &pod_inner_product_, mesh_, degrees,
                     false, filename_primal, groupname_primal,
                     datasetname_primal, out_filename_primal,
                     out_groupname_primal, out_prefix_primal, BLAS, CPU, tol);

      // calculate

      std::cout << "-> Calculation of the basis" << std::endl;
      primal_pod.run();
      // TODO: RELOAD OPTION
      basis_length_primal = primal_pod.get_basis_size();

      // postprocess

      std::cout << "-> Postprocess" << std::endl;

      std::vector< double > ev_primal;
      primal_pod.get_eigenvalues(ev_primal);
      std::string log_ev_primal_filename = "log/ev_primal.txt";
      std::ofstream log_ev_primal(log_ev_primal_filename.c_str(),
                                  std::ofstream::out);
      log_ev_primal.precision(18);
      std::vector< double >::const_iterator it;
      log_ev_primal << "#Eigenvalues of primal POD basis vectors" << std::endl;
      for (it = ev_primal.begin(); it != ev_primal.end(); ++it)
        log_ev_primal << *it << std::endl;
      log_ev_primal.close();

      std::stringstream sv;
      sv << "primal." << adaption_counter;
      primal_pod.get_meanvalue(pod_basis_);
      pp_pod_basis_.UpdateValues();
      prepare_vorticity(pod_basis_);
      visualize_pod_basis(
          sv.str(), -1); // mean value to have index "-1", basis is from 0 to N

      for (int i = 0; i < primal_pod.get_basis_size(); ++i) {
        pod_basis_.Zeros();
        primal_pod.get_basis_vector(i, pod_basis_);
        pp_pod_basis_.UpdateValues();
        visualize_pod_basis(sv.str(), i);
      }

      timer.stop();
      log_data["Primal POD Basis, Time"] = timer.get_duration();

      std::cout << "***************************************************"
                << std::endl;
      std::cout << "PRIMAL POD BASIS COMPUTED (AC=" << adaption_counter << ")"
                << std::endl;
      std::cout << "***************************************************"
                << std::endl;
    }

    // ************************************************************
    // TWICE IS ENOUGH FOR PRIMAL POD BASIS
    // ************************************************************

    const bool doTwice = params_["Twice"]["DoTwice"].get< bool >();

    if (doTwice) {
      std::cout << "***************************************************"
                << std::endl;
      std::cout << "NOW 'TWICE IS ENOUGH' FOR PRIMAL POD BASIS" << std::endl;
      std::cout << "***************************************************"
                << std::endl;

      timer.reset();
      timer.start();

      // Output parameters

      std::string out_filename_twice =
          params_["Twice"]["PrimalBasis"].get< std::string >();
      std::string out_groupname_twice =
          params_["Twice"]["PrimalBasisGroup"].get< std::string >();
      std::string out_prefix_twice =
          params_["Twice"]["PrimalBasisPrefix"].get< std::string >();

      // help vector
      LAD::VectorType help;
      help.Init(comm_, couplings_, la_sys_.Platform, APP_LINALG_IMPLEMENTATION);
      help.InitStructure();
      help.Zeros();

      // Length of Primal Basis
      const int basisLength =
          params_["Twice"]["PrimalBasisLength"].get< int >();

      // Prepare and read primal POD basis
      std::cout << "Prepare and load primal POD basis" << std::endl;

      std::vector< LAD::VectorType * > basis;
      basis.resize(basisLength, new Vector);
      for (int i = 0; i < basisLength; ++i) {
        basis[i] = new LAD::VectorType;
        basis[i]->Init(comm_, couplings_, la_sys_.Platform,
                       APP_LINALG_IMPLEMENTATION);
        basis[i]->InitStructure();
        basis[i]->Zeros();
        basis[i]->UpdateCouplings();
      }

      std::stringstream pb;
      pb << out_prefix_primal << -1;
      basis[0]->ReadHDF5(out_filename_primal, out_groupname_primal, pb.str());

      std::stringstream tb;
      tb << out_prefix_twice << -1;
      basis[0]->WriteHDF5(out_filename_twice, out_groupname_twice, tb.str());

      for (int i = 0; i < basisLength; ++i) {
        std::stringstream pb1;
        pb1 << out_prefix_primal << i;
        basis[i]->ReadHDF5(out_filename_primal, out_groupname_primal,
                           pb1.str());
      }

      // Twice is enough algorithm
      std::cout << "Starting 'Twice is enough' algorithm" << std::endl;

      // loop over all basis vectors
      for (int i = 0; i < basisLength; ++i) {
        // normalize current vector
        help.Zeros();
        pod_inner_product_.VectorMult(*basis[i], &help);
        const double norm_current = std::sqrt(help.Dot(*basis[i]));
        basis[i]->Scale(1. / norm_current);

        // Write result to HDF5
        std::stringstream tb2;
        tb2 << out_prefix_twice << i;
        basis[i]->WriteHDF5(out_filename_twice, out_groupname_twice, tb2.str());

        // apply twice orthogonalization to current vector to all later vectors
        // premultiply current vector with inner product matrix
        help.Zeros();
        pod_inner_product_.VectorMult(*basis[i], &help);
        for (int j = i + 1; j < basisLength; ++j) {
          // first orthogonalization
          const double factor1 = basis[j]->Dot(help);
          basis[j]->Axpy(*basis[i], -1. * factor1);

          // second orthogonalization
          const double factor2 = basis[j]->Dot(help);
          basis[j]->Axpy(*basis[i], -1. * factor2);
        }
      }

      std::cout << "'Twice is enough' algorithm finished" << std::endl;

      // check orthonormality
      std::cout << "Check for orthonormality" << std::endl;
      std::vector< std::string > names;
      names.push_back("i");
      names.push_back("j");
      names.push_back("Inner product");

      std::vector< std::vector< double > > data;
      data.resize(3);
      data[0].resize(0);
      data[1].resize(0);
      data[2].resize(0);

      double check = 0.;

      for (int i = 0; i < basisLength; ++i) {
        help.Zeros();
        pod_inner_product_.VectorMult(*basis[i], &help);
        for (int j = i; j < basisLength; ++j) {
          check = basis[j]->Dot(help);
          data[0].push_back(i);
          data[1].push_back(j);
          data[2].push_back(check);
          data[0].push_back(j);
          data[1].push_back(i);
          data[2].push_back(check);
        }
      }

      std::ostringstream fileerr;
      fileerr << "./out/inner_product_basis.csv";
      if (rank_ == master_rank_) {
        write_csv(names, data, fileerr.str());
      }

      // Visualize computed basis
      std::cout << "Visualize computed basis" << std::endl;

      std::stringstream sv;
      sv << "twice." << adaption_counter;

      for (int i = 0; i < basisLength; ++i) {
        pod_basis_.Zeros();
        pod_basis_.CloneFrom(*basis[i]);
        pp_pod_basis_.UpdateValues();
        visualize_pod_basis(sv.str(), i);
      }

      timer.stop();
      log_data["'Twice is enough', Time"] = timer.get_duration();
    }

    // ************************************************************
    // DUAL PROBLEM
    // ************************************************************

    std::cout << "***************************************************"
              << std::endl;
    std::cout << "NOW THE DUAL PROBLEM (AC=" << adaption_counter << ")"
              << std::endl;
    std::cout << "***************************************************"
              << std::endl;

    timer.reset();
    timer.start();

    std::string dual_mode =
        params_["WeightedErrorEstimator"]["DualMode"].get< std::string >();

    if (dual_mode == "Calculate") {
      std::cout << "-> DualMode = " << dual_mode << std::endl;

      my_local_asm_.set_mode_to_dual();

      if (params_["WeightedErrorEstimator"]["DualOffset"].get< int >()) {
        // DualOffset is the index of the last existing file
        time_step =
            params_["WeightedErrorEstimator"]["DualOffset"].get< int >();
        std::cout << "-> Resuming Simulation at step " << time_step - 1
                  << std::endl;

        // Take care for HDF5 data set names
        for (int i = num_time_steps_; i >= time_step; --i) {
          std::stringstream sp;
          sp << prefix_dual << i;
          datasetname_dual.push_back(sp.str());
        }

        // read last dual solution from file
        // std::stringstream ss7;
        // ss7 << "out/primal." << time_step-1 << ".dat";
        // read_solution_from_file(sol_prev_, ss7.str());
        std::cout << "OPEN FILE " << filename_dual << ", group "
                  << groupname_dual << ", data set "
                  << datasetname_dual[num_time_steps_ - time_step] << std::endl;
        sol_prev_dual_.ReadHDF5(
            filename_primal, groupname_primal,
            datasetname_primal[num_time_steps_ - time_step]);
      } else {
        // initial conditions for dual solution

        // read last primal solution from file
        // std::stringstream ss3;
        // ss3 << "out/primal." << num_time_steps_ << ".dat";
        // read_solution_from_file(sol_, ss3.str());
        sol_.ReadHDF5(filename_primal, groupname_primal,
                      datasetname_primal[time_step]);

        // initial conditions for dual solution
        prepare_ic_dual();

        // GOAL-FUNTIONAL-5: Use primal solution as initial condition for dual
        // problem
        //         std::cout << std::endl;
        //         std::cout << "GOAL-FUNCTIONAL-5:" << std::endl;
        //         std::cout << "-> Global Kinetic Energy" << std::endl;
        //         std::cout << "-> J(u) = 1/2 * (u(T),u(T))_\\Omega" <<
        //         std::endl; std::cout << "-> z(T) = J'(u) = u(T)" <<
        //         std::endl; std::cout << std::endl;
        //         sol_prev_dual_.CloneFrom(sol_);
        //         sol_dual_     .CloneFrom(sol_prev_dual_);

        // store solution for dual problem
        // std::stringstream ss6;
        // ss6 << "out/dual." << time_step << ".dat";
        // save_solution_to_file(sol_prev_dual_, ss6.str());

        // Write solution and take care for HDF5 data set names
        std::stringstream sd;
        sd << prefix_dual << time_step;
        sol_prev_dual_.WriteHDF5(filename_dual, groupname_dual, sd.str());
        datasetname_dual.push_back(sd.str());

        // visualize initial solution
        std::cout << "->  Visualize dual solution at time step " << time_step
                  << " with current time " << delta_t_ * time_step << std::endl;
        visualize_solution_dual(time_step, adaption_counter);

        // Time Stepping Loop

        std::cout << "time_step: " << time_step << std::endl;
        std::cout << "num_time_steps: " << num_time_steps_ << std::endl;
        interminable_assert(time_step == num_time_steps_);
      }

      while (time_step > 0) {
        std::cout << "ntimestep: " << time_step << std::endl;
        std::cout << "ntimesteps: " << num_time_steps_ << std::endl;

        // prepare primal solution

        sol_.CloneFrom(sol_prev_);
        // std::stringstream ss4;
        // ss4 << "out/primal." << time_step-1 << ".dat";
        // read_solution_from_file(sol_prev_, ss4.str());
        sol_prev_.ReadHDF5(filename_primal, groupname_primal,
                           datasetname_primal[time_step - 1]);

        my_local_asm_.set_solP(sol_);
        my_local_asm_.set_solP_prev(sol_prev_);

        // prepare dual solution, i.e. copy sol_prev_dual_ to sol_dual_

        sol_dual_.CloneFrom(sol_prev_dual_);

        my_local_asm_.set_solD(sol_dual_);
        my_local_asm_.set_solD_prev(sol_prev_dual_);

        // solve dual problem

        solve_lp_dual();

        // store solution for dual problem
        // std::stringstream ss5;
        // ss5 << "out/dual." << time_step-1 << ".dat";
        // save_solution_to_file(sol_prev_dual_, ss5.str());

        // Write solution and take care for HDF5 data set names
        std::stringstream sd2;
        sd2 << prefix_dual << time_step - 1;
        sol_prev_dual_.WriteHDF5(filename_dual, groupname_dual, sd2.str());
        datasetname_dual.push_back(sd2.str());

        success = true; // some more meaningful way?
        if (!success) {
          std::cout << "ERROR: Dual problem wasn't solved!" << std::endl;
          break;
        }

        if ((time_step - 1) % visual_step_ == 0) {
          std::cout << "->  Visualizing dual solution at time index "
                    << time_step - 1
                    << " with current time: " << delta_t_ * (time_step - 1)
                    << std::endl;
          visualize_solution_dual(time_step - 1, adaption_counter);
        }

        --time_step;
      }

      my_local_asm_.set_mode_to_primal();
    } else if (dual_mode == "Load") {
      std::cout << "-> DualMode = " << dual_mode << std::endl;

      // Take care for HDF5 data set names
      for (int i = num_time_steps_; i >= 0; --i) {
        std::stringstream sp;
        sp << prefix_dual << i;
        datasetname_dual.push_back(sp.str());
      }
    } else {
      std::cout << "-> DualMode = " << dual_mode << std::endl;
      std::cout << "-> ERROR: " << primal_mode
                << " not a treated parameter choice!" << std::endl;
      quit_program();
    }

    timer.stop();
    log_data["Dual Problem, Time"] = timer.get_duration();

    std::cout << "***************************************************"
              << std::endl;
    std::cout << "DUAL PROBLEM SOLVED ... (AC=" << adaption_counter << ")"
              << std::endl;
    std::cout << "***************************************************"
              << std::endl;

    // ************************************************************
    // POD (Dual problem)
    // ************************************************************

    std::string dual_pod_mode =
        params_["WeightedErrorEstimator"]["DualPODMode"].get< std::string >();

    if (dual_pod_mode == "Calculate") {
      std::cout << "-> DualPODMode = " << dual_pod_mode << std::endl;

      std::cout << "***************************************************"
                << std::endl;
      std::cout << "NOW THE DUAL POD BASIS (AC=" << adaption_counter << ")"
                << std::endl;
      std::cout << "***************************************************"
                << std::endl;

      timer.reset();
      timer.start();

      std::cout << "-> Preparation" << std::endl;

      POD dual_pod(master_rank_, comm_, &pod_inner_product_dual_, mesh_,
                   degrees_dual, false, filename_dual, groupname_dual,
                   datasetname_dual, out_filename_dual, out_groupname_dual,
                   out_prefix_dual, BLAS, CPU, tol);

      // calculate

      std::cout << "-> Calculation of the basis" << std::endl;
      dual_pod.run();
      // TODO: RELOAD OPTION
      basis_length_dual = dual_pod.get_basis_size();

      // postprocess

      std::cout << "-> Postprocess" << std::endl;

      std::vector< double > ev_dual;
      dual_pod.get_eigenvalues(ev_dual);
      std::string log_ev_dual_filename = "log/ev_dual.txt";
      std::ofstream log_ev_dual(log_ev_dual_filename.c_str(),
                                std::ofstream::out);
      log_ev_dual.precision(18);
      std::vector< double >::const_iterator it2;
      log_ev_dual << "#Eigenvalues of dual POD basis vectors" << std::endl;
      for (it2 = ev_dual.begin(); it2 != ev_dual.end(); ++it2)
        log_ev_dual << *it2 << std::endl;
      log_ev_dual.close();

      std::stringstream sv2;
      sv2 << "dual." << adaption_counter;
      dual_pod.get_meanvalue(pod_basis_);
      pp_pod_basis_.UpdateValues();
      visualize_pod_basis(sv2.str(), -1);

      for (int i = 0; i < dual_pod.get_basis_size(); ++i) {
        pod_basis_.Zeros();
        dual_pod.get_basis_vector(i, pod_basis_);
        pp_pod_basis_.UpdateValues();
        visualize_pod_basis(sv2.str(), i);
      }

      timer.stop();
      log_data["Dual POD Basis, Time"] = timer.get_duration();

      std::cout << "***************************************************"
                << std::endl;
      std::cout << "DUAL POD BASIS COMPUTED (AC=" << adaption_counter << ")"
                << std::endl;
      std::cout << "***************************************************"
                << std::endl;
    }

    // ************************************************************
    // DUAL PROBLEM IN POD
    // ************************************************************

    // switch whether dual problem shall be solved in POD or not
    bool solve_dual_pod = params_["DualPOD"]["DoPOD"].get< bool >();

    if (solve_dual_pod) {
      std::cout << "***************************************************"
                << std::endl;
      std::cout << "NOW THE DUAL PROBLEM IN POD (AC=" << adaption_counter << ")"
                << std::endl;
      std::cout << "***************************************************"
                << std::endl;

      timer.reset();
      timer.start();

      // Timing pre-processing
      Timer timer_pre;
      double time_pre = 0.;

      // Timing solution processing
      Timer timer_solve;
      double time_solve = 0.;

      // Timing post-processing
      Timer timer_post;
      double time_post = 0.;

      timer_pre.reset();
      timer_pre.start();

      // set time-step
      time_step = num_time_steps_;

      // filename for output
      std::string filename_dual_pod =
          params_["DualPOD"]["OutputFile"].get< std::string >();

      // Prepare primal and dual projection class

      // Get and set basis length
      int primal_length = params_["DualPOD"]["PrimalBasisLength"].get< int >();
      int dual_length = params_["DualPOD"]["DualBasisLength"].get< int >();

      if (primal_length > 0) {
        basis_length_primal = primal_length;
      }

      if (dual_length > 0) {
        basis_length_dual = dual_length;
      }

      std::string pod_filename_primal =
          params_["DualPOD"]["PrimalBasis"].get< std::string >();
      std::string pod_groupname_primal =
          params_["DualPOD"]["PrimalBasisGroup"].get< std::string >();
      std::string pod_prefix_primal =
          params_["DualPOD"]["PrimalBasisPrefix"].get< std::string >();

      std::string pod_filename_dual =
          params_["DualPOD"]["DualBasis"].get< std::string >();
      std::string pod_groupname_dual =
          params_["DualPOD"]["DualBasisGroup"].get< std::string >();
      std::string pod_prefix_dual =
          params_["DualPOD"]["DualBasisPrefix"].get< std::string >();

      // Prepare datasetnames of basis
      std::vector< std::string > basis_primal;
      for (int i = -1; i < basis_length_primal; ++i) {
        std::stringstream pb;
        pb << pod_prefix_primal << i;
        basis_primal.push_back(pb.str());
      }

      std::vector< std::string > basis_dual;
      for (int i = -1; i < basis_length_dual; ++i) {
        std::stringstream db;
        db << pod_prefix_dual << i;
        basis_dual.push_back(db.str());
      }

      proj_ = new Projection(master_rank_, comm_, &pod_inner_product_, mesh_,
                             degrees, pod_filename_primal, pod_groupname_primal,
                             basis_primal, OPENMP, CPU);
      proj_->init();

      proj_dual_ = new Projection(master_rank_, comm_, &pod_inner_product_dual_,
                                  mesh_, degrees_dual, pod_filename_dual,
                                  pod_groupname_dual, basis_dual, OPENMP, CPU);
      proj_dual_->init();

      // Prepare dual pod data strtuctures
      prepare_dual_pod();

      // read last primal solution from file
      sol_.ReadHDF5(filename_primal, groupname_primal,
                    datasetname_primal[time_step]);
      // project primal solution to primal POD space
      proj_->compute_coefficients(sol_, sol_pod_);
      for (int i = 0; i < basis_length_primal; ++i) {
        sol_prev_pod_[i] = sol_pod_[i];
      }

      // initial conditions for dual solution
      // prepare_ic_dual();
      // TODO: switch back to dual initial condition once prepare_ic_dual()
      // works project initial condition to dual POD space
      sol_dual_.ReadHDF5(filename_dual, groupname_dual, datasetname_dual[0]);
      proj_dual_->compute_coefficients(sol_dual_, sol_dual_pod_);
      proj_dual_->compute_coefficients(sol_dual_, sol_prev_dual_pod_);
      // proj_dual_->compute_coefficients(sol_, sol_dual_pod_);
      // proj_dual_->compute_coefficients(sol_, sol_prev_dual_pod_);

      // recompute projected solutions
      proj_dual_->compute_linear_combination(sol_dual_, sol_dual_pod_);
      proj_dual_->compute_linear_combination(sol_prev_dual_,
                                             sol_prev_dual_pod_);

      // Write solution
      std::stringstream sd;
      sd << prefix_dual << time_step;
      sol_prev_dual_.WriteHDF5(filename_dual_pod, groupname_dual, sd.str());

      // visualize initial solution
      std::cout << "->  Visualize dual solution at time step " << time_step
                << " with current time " << delta_t_ * time_step << std::endl;
      // TODO: better name for POD visualization
      visualize_solution_dual(time_step, -1);

      // prepare nonlinear part of last time-step
      if (rank_ == master_rank_) {
        for (int i = 0; i < basis_length_primal; ++i) {
          nonlinear_dual_sum_prev_pod_.Axpy(*(nonlinear_dual_pod_[i]),
                                            sol_pod_[i]);
        }
      }

      // Time Stepping Loop

      std::cout << "time_step: " << time_step << std::endl;
      std::cout << "num_time_steps_: " << num_time_steps_ << std::endl;

      // Prepare time-stepping parameters
      double theta1 = 0.;
      double theta2 = 0.;
      double theta3 = 0.;
      double theta4 = 0.;
      double theta5 = 0.;
      double theta6 = 0.;

      const int sub_steps = params_["DualPOD"]["SubSteps"].get< int >();

      std::string scheme =
          params_["DualPOD"]["TimeSteppingScheme"].get< std::string >();

      if (scheme == "CrankNicolson") {
        theta1 = 0.5;
        theta2 = 0.5;
        theta3 = 0.5;
        theta4 = 0.5;
        theta5 = 0.5;
        theta6 = 0.5;
        std::cout << "Time-stepping scheme: Crank-Nicolson" << std::endl;
      } else if (scheme == "ImplicitEuler") {
        theta1 = 1.;
        theta3 = 1.;
        theta5 = 1.;
        std::cout << "Time-stepping scheme: Implicit Euler" << std::endl;
      } else if (scheme == "ExplicitEuler") {
        theta2 = 1.;
        theta4 = 1.;
        theta6 = 1.;
        std::string filename_dual_pod =
            params_["DualPOD"]["OutputFile"].get< std::string >();
        std::cout << "Time-stepping scheme: Explicit Euler" << std::endl;
      } else if (scheme == "SemiImplicit") {
        theta1 = 0.5;
        theta2 = 0.5;
        theta3 = 0.5;
        theta4 = 0.5;
        theta6 = 1.;
        std::cout << "Time-stepping scheme: Semi Implicit" << std::endl;
      } else {
        std::cout << "ERROR: Invalid time-stepping scheme!" << std::endl;
        exit(-1);
      }

      const double fac1 = delta_t_ * theta1 / static_cast< double >(sub_steps);
      const double fac2 = -delta_t_ * theta2 / static_cast< double >(sub_steps);
      const double fac3 = -delta_t_ * theta3 / static_cast< double >(sub_steps);
      const double fac4 = -delta_t_ * theta4 / static_cast< double >(sub_steps);
      const double fac5 = delta_t_ * theta5 / static_cast< double >(sub_steps);
      const double fac6 = -delta_t_ * theta6 / static_cast< double >(sub_steps);

      timer_pre.stop();
      time_pre += timer_pre.get_duration();

      while (time_step > 0) {
        timer_pre.reset();
        timer_pre.start();

        // prepare primal solution
        for (int i = 0; i < basis_length_primal; ++i) {
          sol_pod_[i] = sol_prev_pod_[i];
          sol_prev_pod_[i] = 0.;
        }
        sol_prev_.ReadHDF5(filename_primal, groupname_primal,
                           datasetname_primal[time_step - 1]);
        proj_->compute_coefficients(sol_prev_, sol_prev_pod_);

        // prepare interpolated primal solution
        std::vector< double > sol_prev_ip_pod_;
        sol_prev_ip_pod_.resize(basis_length_primal, 0.);

        timer_pre.stop();
        time_pre += timer_pre.get_duration();

        // loop through sub_steps
        for (int j = 0; j < sub_steps; ++j) {
          timer_solve.reset();
          timer_solve.start();

          // prepare dual solution, i.e. copy sol_prev_dual_ to sol_dual_
          for (int i = 0; i < basis_length_dual; ++i) {
            sol_dual_pod_[i] = sol_prev_dual_pod_[i];
          }

          const double w1 =
              static_cast< double >(j + 1) / static_cast< double >(sub_steps);
          const double w2 = static_cast< double >(sub_steps - (j + 1)) /
                            static_cast< double >(sub_steps);
          for (int i = 0; i < basis_length_primal; ++i) {
            sol_prev_ip_pod_[i] = w1 * sol_prev_pod_[i] + w2 * sol_pod_[i];
          }

          //***************************************************
          // solve dual problem
          // TODO: source out into own function?

          // copy source_prev to source
          for (int i = 0; i < basis_length_dual; ++i) {
            source_dual_pod_[i] = source_prev_dual_pod_[i];
          }

          // assemble source_prev
          // TODO: adapt for more general goal functionals
          for (int i = 0; i < basis_length_dual; ++i) {
            source_prev_dual_pod_[i] = 0.;
          }

          if (rank_ == master_rank_) {
            // copy nonlinear part of last time-step
            nonlinear_dual_sum_pod_.Zeros();
            nonlinear_dual_sum_pod_.Add(nonlinear_dual_sum_prev_pod_);

            // prepare right hand side
            std::vector< double > help;
            help.resize(basis_length_dual, 0.);
            rhs_dual_pod_.resize(basis_length_dual, 0.);

            mass_dual_pod_.VectorMult(sol_dual_pod_, rhs_dual_pod_);

            stiff_dual_pod_.VectorMult(sol_dual_pod_, help);
            for (int i = 0; i < basis_length_dual; ++i) {
              rhs_dual_pod_[i] += fac2 * help[i];
            }

            help.resize(basis_length_dual, 0.);
            nonlinear_dual_sum_pod_.VectorMult(sol_dual_pod_, help);
            for (int i = 0; i < basis_length_dual; ++i) {
              rhs_dual_pod_[i] += fac6 * help[i];
            }

            help.resize(basis_length_dual, 0.);
            mass_dual_pod_.VectorMult(source_dual_pod_, help);
            for (int i = 0; i < basis_length_dual; ++i) {
              rhs_dual_pod_[i] += fac4 * help[i];
            }

            help.resize(basis_length_dual, 0.);
            mass_dual_pod_.VectorMult(source_prev_dual_pod_, help);
            for (int i = 0; i < basis_length_dual; ++i) {
              rhs_dual_pod_[i] += fac3 * help[i];
            }

            // compute nonlinear part of current time-step
            nonlinear_dual_sum_prev_pod_.Zeros();
            for (int i = 0; i < basis_length_primal; ++i) {
              nonlinear_dual_sum_prev_pod_.Axpy(*(nonlinear_dual_pod_[i]),
                                                sol_prev_ip_pod_[i]);
            }

            // compute iteration matrix
            iteration_dual_pod_.Zeros();
            iteration_dual_pod_.Add(mass_dual_pod_);
            iteration_dual_pod_.Axpy(stiff_dual_pod_, fac1);
            iteration_dual_pod_.Axpy(nonlinear_dual_sum_prev_pod_, fac5);

            // solve linear system
            sol_prev_dual_pod_.resize(basis_length_dual, 0.);
            iteration_dual_pod_.Solve(rhs_dual_pod_, sol_prev_dual_pod_);
          }

          timer_solve.stop();
          time_solve += timer_solve.get_duration();

          timer_post.reset();
          timer_post.start();

          // Communicate current solution
          MPI_Bcast(&sol_prev_dual_pod_[0], basis_length_dual, MPI_DOUBLE,
                    master_rank_, comm_);

          // recompute dual solution in full FEM system
          sol_prev_dual_.Zeros();
          proj_dual_->compute_linear_combination(sol_prev_dual_,
                                                 sol_prev_dual_pod_);

          //***************************************************

          // Write solution
          std::stringstream sd2;
          sd2 << prefix_dual << time_step - 1 << "." << sub_steps - j;
          sol_prev_dual_.WriteHDF5(filename_dual_pod, groupname_dual,
                                   sd2.str());

          timer_post.stop();
          time_post += timer_post.get_duration();
        }

        timer_post.reset();
        timer_post.start();

        if ((time_step - 1) == 0) {
          std::cout << "->  Visualizing dual solution at time index "
                    << time_step - 1
                    << " with current time: " << delta_t_ * (time_step - 1)
                    << std::endl;
          visualize_solution_dual(time_step - 1, -1);
        }

        timer_post.stop();
        time_post += timer_post.get_duration();

        --time_step;
      }

      timer.stop();
      log_data["Dual Problem in POD, Time"] = timer.get_duration();
      log_data["Dual Problem in POD, Pre-processing Time"] = time_pre;
      log_data["Dual Problem in POD, Solving Time"] = time_solve;
      log_data["Dual Problem in POD, Post-processing Time"] = time_post;
    }

    // *************************************************************
    // ERROR COMPUTATION (POD SOLUTION TO GIVEN REFERENCE SOLUTION)
    // *************************************************************

    // read parameter
    bool doError = params_["ErrorComputation"]["DoError"].get< bool >();

    if (doError) {
      std::cout << "***********************************************************"
                   "***********************"
                << std::endl;
      std::cout << "COMPUTE L2-ERROR OF POD SOLUTION TO GIVEN REFERENCE "
                   "SOLUTION (DUAL PROBLEM) (AC="
                << adaption_counter << ")" << std::endl;
      std::cout << "***********************************************************"
                   "***********************"
                << std::endl;

      timer.reset();
      timer.start();

      // Construct and initialize projection to primal POD space
      std::vector< double > sol_dual_projection;
      sol_dual_projection.resize(basis_length_primal, 0.);

      // Prepare data structures for output
      std::vector< std::string > names;
      names.push_back("Time");
      names.push_back("L2 error (projection, relative)");
      names.push_back("L2 error (absolute)");
      names.push_back("L2 error (relative)");
      names.push_back("L2 error (scaled)");
      names.push_back("Scaling factor (lambda)");

      std::vector< std::vector< double > > data;
      data.resize(6);
      data[0].resize(0);
      data[1].resize(0);
      data[2].resize(0);
      data[3].resize(0);
      data[4].resize(0);
      data[5].resize(0);

      // parameters for reference solution
      std::string filename_reference =
          params_["ErrorComputation"]["ReferenceSnapshots"]
              .get< std::string >();
      std::string groupname_reference =
          params_["ErrorComputation"]["ReferenceSnapshotsGroup"]
              .get< std::string >();
      std::string prefix_reference =
          params_["ErrorComputation"]["ReferenceSnapshotsPrefix"]
              .get< std::string >();

      const int timesteps_reference =
          params_["ErrorComputation"]["ReferenceNumberOfTimesteps"]
              .get< int >();
      const double delta_t_reference =
          params_["ErrorComputation"]["ReferenceTimeStep"].get< double >();

      const double max_time_reference = timesteps_reference * delta_t_reference;

      // read parameters of POD solution
      const int sub_steps = params_["DualPOD"]["SubSteps"].get< int >();
      std::string filename_dual_pod =
          params_["DualPOD"]["OutputFile"].get< std::string >();

      // help vectors
      LAD::VectorType sol_dual_reference, sol_prev_dual_reference,
          sol_dual_reference_current, difference, help, projection;

      sol_dual_reference.Init(comm_, couplings_dual_, la_sys_.Platform,
                              APP_LINALG_IMPLEMENTATION);
      sol_dual_reference.InitStructure();
      sol_dual_reference.Zeros();

      sol_prev_dual_reference.Init(comm_, couplings_dual_, la_sys_.Platform,
                                   APP_LINALG_IMPLEMENTATION);
      sol_prev_dual_reference.InitStructure();
      sol_prev_dual_reference.Zeros();

      sol_dual_reference_current.Init(comm_, couplings_dual_, la_sys_.Platform,
                                      APP_LINALG_IMPLEMENTATION);
      sol_dual_reference_current.InitStructure();
      sol_dual_reference_current.Zeros();

      difference.Init(comm_, couplings_dual_, la_sys_.Platform,
                      APP_LINALG_IMPLEMENTATION);
      difference.InitStructure();
      difference.Zeros();

      help.Init(comm_, couplings_dual_, la_sys_.Platform,
                APP_LINALG_IMPLEMENTATION);
      help.InitStructure();
      help.Zeros();

      projection.Init(comm_, couplings_dual_, la_sys_.Platform,
                      APP_LINALG_IMPLEMENTATION);
      projection.InitStructure();
      projection.Zeros();

      // is error computation done?
      bool done = false;

      // loop through time-steps of POD solution
      for (int i = 0; i <= num_time_steps_; ++i) {
        if (!done) {
          // loop through sub-steps per time-step
          for (int j = 0; j < sub_steps; ++j) {
            // compute current time
            double current_time = i * delta_t_ + j * (delta_t_ / sub_steps);

            // check whether current_time is within time-interval of reference
            // solution
            if (current_time < max_time_reference) {
              // load POD solution at current time
              std::stringstream sd;
              sd << prefix_dual << i << "." << (j + 1);
              sol_dual_.Zeros();
              sol_dual_.ReadHDF5(filename_dual_pod, groupname_dual, sd.str());

              // compute indices of reference snapshots embracing current time
              int lower_index = static_cast< int >(
                  std::floor(current_time / delta_t_reference));
              int upper_index = lower_index + 1;

              // load embracing snapshots
              std::stringstream lower;
              lower << prefix_reference << lower_index;
              std::stringstream upper;
              upper << prefix_reference << upper_index;

              sol_prev_dual_reference.Zeros();
              sol_prev_dual_reference.ReadHDF5(
                  filename_reference, groupname_reference, lower.str());
              sol_dual_reference.Zeros();
              sol_dual_reference.ReadHDF5(filename_reference,
                                          groupname_reference, upper.str());

              // compute interpolation weights
              const double time_upper = upper_index * delta_t_reference;
              const double w2 = (time_upper - current_time) / delta_t_reference;
              const double w1 = 1. - w2;

              // compute interpolation
              sol_dual_reference_current.Zeros();
              sol_dual_reference_current.CloneFrom(sol_prev_dual_reference);
              sol_dual_reference_current.Scale(w2);
              sol_dual_reference_current.Axpy(sol_dual_reference, w1);

              // compute projection
              sol_dual_projection.resize(full_basis_length, 0.);
              projection.Zeros();
              proj_->compute_coefficients(sol_dual_reference_current,
                                          sol_dual_projection);
              proj_->compute_linear_combination(projection,
                                                sol_dual_projection);

              difference.Zeros();
              difference.CloneFrom(projection);
              difference.Axpy(sol_dual_reference_current, -1.);
              help.Zeros();
              pod_inner_product_dual_.VectorMult(difference, &help);
              const double projection_abs = std::sqrt(help.Dot(difference));

              // compute L2-error
              difference.Zeros();
              difference.CloneFrom(sol_dual_);
              difference.Axpy(sol_dual_reference_current, -1.);

              // compute absolute L2-error
              help.Zeros();
              pod_inner_product_dual_.VectorMult(difference, &help);
              const double L2_abs = std::sqrt(help.Dot(difference));

              // compute relative L2-error
              help.Zeros();
              pod_inner_product_dual_.VectorMult(sol_dual_reference_current,
                                                 &help);
              const double L2_ref_norm =
                  std::sqrt(help.Dot(sol_dual_reference_current));

              const double L2_rel = L2_abs / L2_ref_norm;
              const double projection_rel = projection_abs / L2_ref_norm;

              // compute scaled L2-error
              help.Zeros();
              pod_inner_product_dual_.VectorMult(sol_dual_, &help);
              const double L2_current = std::sqrt(help.Dot(sol_dual_));

              const double lambda = L2_ref_norm / L2_current;

              difference.Zeros();
              difference.CloneFrom(sol_dual_reference_current);
              difference.Axpy(sol_dual_, -1. * lambda);

              help.Zeros();
              pod_inner_product_dual_.VectorMult(difference, &help);
              const double L2_scaled =
                  std::sqrt(help.Dot(difference)) / L2_ref_norm;

              data[0].push_back(current_time);
              data[1].push_back(projection_rel);
              data[2].push_back(L2_abs);
              data[3].push_back(L2_rel);
              data[4].push_back(L2_scaled);
              data[5].push_back(lambda);
            } else if (current_time == max_time_reference) {
              // load POD solution at current time
              std::stringstream sd;
              sd << prefix_dual << i << "." << (j + 1);
              sol_dual_.Zeros();
              sol_dual_.ReadHDF5(filename_dual_pod, groupname_dual, sd.str());

              // load embracing snapshots
              std::stringstream upper;
              upper << prefix_reference << timesteps_reference;
              sol_dual_reference.Zeros();
              sol_dual_reference.ReadHDF5(filename_reference,
                                          groupname_reference, upper.str());

              // compute projection
              sol_dual_projection.resize(basis_length_primal, 0.);
              projection.Zeros();
              proj_->compute_coefficients(sol_dual_reference_current,
                                          sol_dual_projection);
              proj_->compute_linear_combination(projection,
                                                sol_dual_projection);

              difference.Zeros();
              difference.CloneFrom(projection);
              difference.Axpy(sol_dual_reference_current, -1.);
              help.Zeros();
              pod_inner_product_dual_.VectorMult(difference, &help);
              const double projection_abs = std::sqrt(help.Dot(difference));

              // compute L2-error
              difference.Zeros();
              difference.CloneFrom(sol_dual_);
              difference.Axpy(sol_dual_reference, -1.);

              // compute absolute L2-error
              help.Zeros();
              pod_inner_product_dual_.VectorMult(difference, &help);
              const double L2_abs = std::sqrt(help.Dot(difference));

              // compute relative L2-error
              help.Zeros();
              pod_inner_product_dual_.VectorMult(sol_dual_reference, &help);
              const double L2_ref_norm =
                  std::sqrt(help.Dot(sol_dual_reference));

              const double L2_rel = L2_abs / L2_ref_norm;
              const double projection_rel = projection_abs / L2_ref_norm;

              // compute scaled L2-error
              help.Zeros();
              pod_inner_product_dual_.VectorMult(sol_dual_, &help);
              const double L2_current = std::sqrt(help.Dot(sol_dual_));

              const double lambda = L2_ref_norm / L2_current;

              help.Zeros();
              difference.CloneFrom(sol_dual_reference_current);
              difference.Axpy(sol_dual_, -1. * lambda);

              help.Zeros();
              pod_inner_product_dual_.VectorMult(difference, &help);
              const double L2_scaled =
                  std::sqrt(help.Dot(difference)) / L2_ref_norm;

              data[0].push_back(current_time);
              data[1].push_back(projection_rel);
              data[2].push_back(L2_abs);
              data[3].push_back(L2_rel);
              data[4].push_back(L2_scaled);
              data[5].push_back(lambda);
            } else {
              done = true;
            }
          }
        }
      }

      std::ostringstream fileerr;
      std::string filename_error =
          params_["ErrorComputation"]["CSVOutput"].get< std::string >();
      fileerr << filename_error;
      if (rank_ == master_rank_) {
        write_csv(names, data, fileerr.str());
      }

      timer.stop();
      log_data["Error Computation (POD Solution), Time"] = timer.get_duration();
    }

    // write time-statistics

    std::string log_filename = "log/time_statistics.txt";
    std::ofstream log(log_filename.c_str(), std::ofstream::out);
    if (log.is_open() == false)
      std::cout << "IOPROBLEM: cannot open file '" << log_filename << "'!"
                << std::endl;
    log.precision(18);
    std::map< std::string, double >::const_iterator it3;
    log << "#Time-Statistics" << std::endl;
    for (it3 = log_data.begin(); it3 != log_data.end(); ++it3)
      log << "'" << it3->first << "',\t" << it3->second << std::endl;
    log.close();

    quit_program();

    // COMPUTE INDICATORS

    write_indicators(adaption_counter);
    visualize_indicators(adaption_counter);
    write_log_file(adaption_counter);

    ++adaption_counter;

    if (!success) {
      int restart_point;
      if (time_step > 2 * indicator_step_)
        restart_point =
            time_step - (time_step % indicator_step_) - indicator_step_;
      else
        restart_point = time_step - (time_step % indicator_step_);

      std::cout << " LOADING STATE OF KILLED JOB at " << restart_point
                << "ind count " << adaption_counter << "\n";
      load_temp_indicator(adaption_counter - 1, restart_point);
    }

    if (adaption_counter <= num_adaptions_) {
      // Adapt mesh
      if (refinement_criterion_ == "ns_aposteriori") {
        // adapt_mesh_const_num_strategy(adaption_counter, ns_apost);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, ns_apost);
      } else if (refinement_criterion_ == "aposteriori") {
        // adapt_mesh_const_num_strategy(adaption_counter, apost);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, apost);
      } else if (refinement_criterion_ == "energy") {
        // adapt_mesh_const_num_strategy(adaption_counter, energy);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, energy);
      } else if (refinement_criterion_ == "vorticity") {
        // adapt_mesh_const_num_strategy(adaption_counter, vort);
        adapt_mesh_fixed_fraction_strategy(adaption_counter, vort);
      } else {
        std::cout << "Refinement Criterion unknown!\n";
        break;
      }
    }
  }
}

void MetFlowCCIApp::postprocessing_run() {
  // Primal

  // For HDF5
  std::string filename_primal =
      params_["WeightedErrorEstimator"]["PrimalSnapshots"].get< std::string >();
  std::string groupname_primal = "solution_primal";
  std::string prefix_primal = "snapshot_primal_";
  std::vector< std::string > datasetname_primal;

  // setup finite element ansatz (PRIMAL PROBLEM)
  const int velocity_deg =
      params_["FiniteElements"]["VelocityDegree"].get< int >();
  const int pressure_deg =
      params_["FiniteElements"]["PressureDegree"].get< int >();
  std::vector< int > degrees(DIM + 1, velocity_deg);
  degrees[DIM] = pressure_deg;

  // Prepare application
  prepare();

  // DO SOMETHING
  datasetname_primal.push_back("snapshot_primal_0");
  std::cout << "READING file '" << filename_primal << "'\t groupname '"
            << groupname_primal << "'\t datasetname '" << datasetname_primal[0]
            << "'" << std::endl;
  sol_.ReadHDF5(filename_primal, groupname_primal, datasetname_primal[0]);

  std::cout << "WRITING file with index " << -1 << std::endl;
  pp_sol_.UpdateValues();
  prepare_vorticity();
  compute_vorticity();
  MetFlowApp::visualize_solution(-1);

  // Dual

  // For HDF5
  std::string filename_dual =
      params_["WeightedErrorEstimator"]["DualSnapshots"].get< std::string >();
  std::string groupname_dual = "solution_dual";
  std::string prefix_dual = "snapshot_dual_";
  std::vector< std::string > datasetname_dual;

  // setup finite element ansatz (DUAL PROBLEM)
  const int velocity_deg_dual =
      params_["FiniteElementsDual"]["VelocityDegree"].get< int >();
  const int pressure_deg_dual =
      params_["FiniteElementsDual"]["PressureDegree"].get< int >();
  std::vector< int > degrees_dual(DIM + 1, velocity_deg_dual);
  degrees_dual[DIM] = pressure_deg_dual;

  // Prepare application
  prepare_dual();

  // DO SOMETHING
  datasetname_dual.push_back("snapshot_dual_0");
  std::cout << "READING file '" << filename_dual << "'\t groupname '"
            << groupname_dual << "'\t datasetname '" << datasetname_dual[0]
            << "'" << std::endl;
  sol_prev_dual_.ReadHDF5(filename_dual, groupname_dual, datasetname_dual[0]);
  sol_dual_.CloneFrom(sol_prev_dual_);

  std::cout << "WRITING file with index " << -1 << std::endl;
  pp_sol_dual_.UpdateValues();
  compute_vorticity(sol_prev_dual_);
  MetFlowDualApp::visualize_solution_dual(-1);
}

// restart after a given number of adaptation cycles.

void MetFlowCCIApp::restart() {
  prepare_parameters();
  int restart_at_ = params_["Mesh"]["RestartAt"].get< int >();
  for (int count = 0; count < restart_at_; count++) {
    // load indicator file in memory
    std::cout << "Loading old state of adaption count " << count << "\n";
    load_indicator(count);
    std::cout << "indicator has size " << apost.size() << " and mesh is "
              << mesh_->num_entities(DIM) << "\n\n\n";
    // Prepare application
    MetFlowApp::prepare_parameters();
    MetFlowApp::prepare_space();
    periodify_space();

    // Adapt mesh
    if (refinement_criterion_ == "ns_aposteriori") {
      // adapt_mesh_const_num_strategy(count, ns_apost);
      adapt_mesh_fixed_fraction_strategy(count, ns_apost);
    } else if (refinement_criterion_ == "aposteriori") {
      // adapt_mesh_const_num_strategy(count, apost);
      adapt_mesh_fixed_fraction_strategy(count, apost);
    } else if (refinement_criterion_ == "energy") {
      // adapt_mesh_const_num_strategy(count, energy);
      adapt_mesh_fixed_fraction_strategy(count, energy);
    } else if (refinement_criterion_ == "vorticity") {
      // adapt_mesh_const_num_strategy(count, vort);
      adapt_mesh_fixed_fraction_strategy(count, vort);
    } else {
      std::cout << "Refinement Criterion unknown!\n";
    }
    visualize_mesh(count);
  }
}

// Prepare parameters

void MetFlowCCIApp::prepare_parameters() {
  MetFlowApp::prepare_parameters();

  std::cout << "reading MetFlowCCIApp-parameters" << std::endl;
  Um_ = params_["CCI"]["BackroundDrift"].get< double >();
  distance_ = params_["CCI"]["Distance"].get< double >();
  delta_t_ = params_["TimeDiscretization"]["DeltaT"].get< double >();
  num_time_steps_ = params_["TimeDiscretization"]["Timesteps"].get< int >();
  visual_step_ = params_["CCI"]["VisualStep"].get< int >();
  track_step_ = params_["CCI"]["StormTrackStep"].get< int >();
  indicator_step_ = params_["CCI"]["IndicatorStep"].get< int >();
  refinement_criterion_ =
      (params_["Mesh"]["AdaptationCriterion"].get< std::string >()).c_str();

  // set up time parameters
  method_ =
      (params_["TimeDiscretization"]["Method"].get< std::string >()).c_str();

  // set goal-functional
  goal_contribution_ = params_["Goal"]["Contribution"].get< std::string >();
}

/// Prepare and set initial conditions of vortex-vortex-interaction scenario,
/// including projection into divergence-free space (sol_ and sol_prev_).

void MetFlowCCIApp::prepare_ic() {
  for (int var = 0; var < DIM; ++var) {
    const int tdim = space_.mesh().tdim();
    for (mesh::EntityIterator it = space_.mesh().begin(tdim),
                              end_it = space_.mesh().end(tdim);
         it != end_it; ++it) {
      std::vector< int > global_dof_ids;
      space_.GetDofIndices(var, *it, &global_dof_ids);
      int num_dofs = global_dof_ids.size();
      std::vector< double > values;
      values.resize(num_dofs);

      // TODO: handle initial conditions at periodic boundaries, e.g. take
      // arithmetic mean instead of Master's value as currently implemented.
      std::vector< doffem::Coord > coords;
      space_.dof().get_coord_on_cell(var, it->index(), coords);
      for (int j = 0; j < num_dofs; j++) {
        Coord pt = coords[j];
        double _distance = distance_;
        double x = pt[0];
        double y = pt[1];
        double xv_1 = -_distance / 2.;
        double xv_2 = _distance / 2.;
        double yv = 0.;

        double velocity = 71.521 / 1000.; // [km/sec]
        double radius = 100.0;            // [km]
        double a = 0.3398;
        double b = 5.377e-4;
        double result = 0;

        switch (var) {
        case 0: {
          result += add_sud_vortex(x, y, xv_1, yv, velocity, radius, a, b, 0);
          result += add_sud_vortex(x, y, xv_2, yv, velocity, radius, a, b, 0);
          result += Um_ / 1000.;
          break;
        }
        case 1: {
          result += add_sud_vortex(x, y, xv_1, yv, velocity, radius, a, b, 1);
          result += add_sud_vortex(x, y, xv_2, yv, velocity, radius, a, b, 1);
          break;
        }
        case 2: {
          result += 0.0;
          break;
        }
        default: {
          quit_program();
          break;
        }
        }

        values[j] = result;
      }

      sol_.SetValues(vec2ptr(global_dof_ids), num_dofs, vec2ptr(values));
    }
  }

  // project sol_ into space of divergence free solutions
  div_free_velocity();

  // copy sol_ to sol_prev_
  sol_prev_.CloneFrom(sol_);
}

/// prepare and set initial conditions of the dual problem,
/// including projection into divergence-free space (sol_dual_ and
/// sol_prev_dual_).

void MetFlowCCIApp::prepare_ic_dual() {
  if (goal_functional_->initial_type_active() == false) {
    // initialize with zero solution
    sol_prev_dual_.Zeros();
  } else {
    // L2-representation of initial condition as given by goal-functional
    // followed by projection into divergence-free space

    // assemble mass matrix
    SparsityStructure sparsity;
    global_asm_.compute_sparsity_structure(space_dual_, sparsity);
    LAD::MatrixType mass_dual;
    mass_dual.Init(comm_, couplings_dual_, la_sys_.Platform,
                   MetFlowDualApp::APP_LINALG_IMPLEMENTATION,
                   APP_MATRIX_FORMAT);
    mass_dual.InitStructure(
        vec2ptr(sparsity.diagonal_rows), vec2ptr(sparsity.diagonal_cols),
        sparsity.diagonal_rows.size(), vec2ptr(sparsity.off_diagonal_rows),
        vec2ptr(sparsity.off_diagonal_cols), sparsity.off_diagonal_rows.size());
    L2InnerProduct L2_asm;
    global_asm_.assemble_matrix(space_dual_, L2_asm, mass_dual);

    // assemble rhs
    DualRHSAssembler< DIM > DualRHS_asm(goal_functional_);
    DualRHS_asm.set_solP(sol_);
    global_asm_.assemble_vector(space_dual_, DualRHS_asm, rhs_dual_);
    rhs_dual_.UpdateCouplings();

    // solve linear system
    ScopedPtr< GMRES< LAD > >::Type linear_solver(new GMRES< LAD >());
    linear_solver->InitControl(500, 1.e-12, 1.e-12, 1.e6);
    linear_solver->InitParameter(70, "NoPreconditioning");
    linear_solver->SetupOperator(mass_dual);
    linear_solver->Solve(rhs_dual_, &sol_prev_dual_);
    interpolate_constrained_vector(space_dual_, sol_prev_dual_);
    sol_prev_dual_.UpdateCouplings();

    // project into space of divergence free solutions
    div_free_dual_velocity(sol_prev_dual_);
  }

  sol_dual_.CloneFrom(sol_prev_dual_);

  /*
  for (int var = 0; var < DIM; ++var)
  {
    const int tdim = space_.mesh().tdim();
    for (mesh::EntityIterator it = space_.mesh().begin(tdim), end_it =
  space_.mesh().end(tdim); it != end_it;
         ++it)
    {
      std::vector<int>  global_dof_ids;
      space_.GetDofIndices(var, *it, &global_dof_ids);
      int num_dofs = global_dof_ids.size();
      std::vector<double> values;
      values.resize(num_dofs);

      // TODO: handle initial conditions at periodic boundaries, e.g. take
  arithmetic mean instead of Master's value as currently implemented.
      std::vector<doffem::Coord> coords;
      space_.dof().get_coord_on_cell (var, it->index(), coords);
      for (int j= 0; j < num_dofs; j++)
      {
        Coord pt = coords[j];
        double _distance = distance_;
        double x = pt[0];
        double y = pt[1];
        double xv_1 = - _distance/2.;
        double xv_2 =  _distance/2.;
        double yv   = 0.;

        // MANIMUPATION: use only half of the distance!!! and other place!!!
        xv_1 = - _distance/4.;
        xv_2 =  _distance/4.;
        xv_1-=1000.;
        xv_2-=1000.;
        yv  -=1000.;

        double velocity = 71.521/1000.; // [km/sec]
        double radius   = 100.0;       // [km]
        double a        = 0.3398;
        double b        = 5.377e-4;
        double result   = 0;

        switch(var)
        {
          case 0:
          {
            result += add_sud_vortex(x, y, xv_1, yv, velocity, radius, a, b, 0);
            result += add_sud_vortex(x, y, xv_2, yv, velocity, radius, a, b, 0);
            result += Um_/1000.;
            break;
          }
          case 1:
          {
            result += add_sud_vortex(x, y, xv_1, yv, velocity, radius, a, b, 1);
            result += add_sud_vortex(x, y, xv_2, yv, velocity, radius, a, b, 1);
            break;
          }
          case 2:
          {
            result += 0.0;
            break;
          }
          default:
          {
            quit_program();
            break;
          }
        }

        values[j] = result;
      }

      sol_prev_dual_.SetValues(vec2ptr(global_dof_ids), num_dofs,
  vec2ptr(values));
    }
  }

  // project into space of divergence free solutions
  div_free_dual_velocity(sol_prev_dual_);
  sol_dual_.CloneFrom(sol_prev_dual_);*/
}

/// add velocity profile of an idealized sud vortex to velocity profile

double MetFlowCCIApp::add_sud_vortex(double x,        // coord x
                                     double y,        // coord y
                                     double xv,       // vortex center
                                     double yv,       // vortex center
                                     double velocity, // reference velocity
                                     double radius,   // radius of maximal wind
                                     double a,        // parameter ?
                                     double b,        // parameter ?
                                     int var)         // component
{
  double r = sqrt((x - xv) * (x - xv) +
                  (y - yv) * (y - yv)); // distance to vortex center
  double s = r / radius;                // relative distance to vortex center

  // Choke sud_vortex so that vortex profile equals zero if r > rcut
  double rcut = 1000.;
  double choke;
  if (r < rcut) {
    choke = 1. - exp(-(pow((r - rcut), 2.0)) / pow(100., 2.0));
  } else {
    choke = 0.;
  }
  choke = 1.;
  double v_tan = choke * velocity * (s * (1 + (6 * b / (2 * a)) * pow(s, 4))) /
                 pow(1 + a * s * s + b * pow(s, 6), 2);
  double alpha = atan2(y - yv, x - xv);
  double result;

  switch (var) {
  case 0:
    result = -v_tan * sin(alpha);
    break;
  case 1:
    result = v_tan * cos(alpha);
    break;
  case 2:
    result = 0.0;
    break;
  default:
    quit_program();
    break;
  }

  return result;
}

/// structure for the construction of dirichlet boundaries

struct DirichletBC2d {

  DirichletBC2d(int var, double Um) : var_(var), Um_(Um) {
    assert(var_ == 0 || var_ == 1);
  }

  std::vector< double >
  evaluate(const Entity &face,
           const std::vector< Coord > &coords_on_face) const {
    std::vector< double > values;
    values.resize(coords_on_face.size());
    for (int i = 0; i < coords_on_face.size(); ++i) {
      values[i] = 0;
    }
    return values;
  }

  const int var_;
  const double Um_;
};

/// prepare dirichlet boundary conditions

void MetFlowCCIApp::prepare_bc() {
  if (period.size() == 2)
    return dirichlet_dofs_.clear();
  dirichlet_values_.clear();

  DirichletBC2d bc[2] = {DirichletBC2d(0, Um_), DirichletBC2d(1, Um_)};

  for (int var = 0; var < DIM; ++var) {
    compute_dirichlet_dofs_and_values(bc[var], space_, var, dirichlet_dofs_,
                                      dirichlet_values_);
  }
}

void MetFlowCCIApp::prepare_indicators() {
  int num_cells = space_.mesh().num_entities(space_.mesh().tdim());
  vort.clear();
  energy.clear();
  apost.clear();
  ns_apost.clear();
  sum_ns_apost.clear();
  sum_apost.clear();
  vort.resize(num_cells, 0.);
  energy.resize(num_cells, 0.);
  apost.resize(num_cells, 0.);
  ns_apost.resize(num_cells, 0.);
  sum_apost.resize(num_cells, 0.);
  sum_ns_apost.resize(num_cells, 0.);
}

// Estimate error indicators for residual error
// in Poisson-problem and Navier-Stokes

void MetFlowCCIApp::estimate_indicators() {

  // A posteriori residual error estimator
  const int tdim = space_.mesh().tdim();
  EntityCount num_cells = space_.mesh().num_entities(tdim);

  std::vector< double > po_error_indicators, ns_error_indicators,
      vorticity_indicators, energy_indicators, barycenter_x, barycenter_y;

  po_error_indicators.resize(num_cells, 0.);
  ns_error_indicators.resize(num_cells, 0.);
  vorticity_indicators.resize(num_cells, 0.);
  energy_indicators.resize(num_cells, 0.);
  barycenter_x.resize(num_cells, 0.);
  barycenter_y.resize(num_cells, 0.);

  double alpha_ = params_["CCI"]["Alpha"].get< double >();
  AdaptationIndicator< DIM > indicator(space_, sol_);
  indicator.evaluate_vorticity_indicator(vorticity_indicators, alpha_);
  indicator.evaluate_energy_indicator(energy_indicators, alpha_);
  indicator.evaluate_poisson_error_indicator(po_error_indicators);
  indicator.evaluate_ns_error_indicator(ns_error_indicators, nu_);
  indicator.evaluate_barycenter(barycenter_x, 0);
  indicator.evaluate_barycenter(barycenter_y, 1);

  for (mesh::EntityIterator it = space_.mesh().begin(tdim),
                            end_it = space_.mesh().end(tdim);
       it != end_it; ++it) {
    sum_ns_apost.at(it->index()) += ns_error_indicators.at(it->index());
    sum_apost.at(it->index()) += po_error_indicators.at(it->index());

    // Find maximum in time!
    if (ns_apost.at(it->index()) < ns_error_indicators.at(it->index()))
      ns_apost.at(it->index()) = ns_error_indicators.at(it->index());
    if (apost.at(it->index()) < po_error_indicators.at(it->index()))
      apost.at(it->index()) = po_error_indicators.at(it->index());
    if (vort.at(it->index()) < vorticity_indicators.at(it->index()))
      vort.at(it->index()) = vorticity_indicators.at(it->index());
    if (energy.at(it->index()) < energy_indicators.at(it->index()))
      energy.at(it->index()) = energy_indicators.at(it->index());
  }
}

// Load an existing set of indicator values in memory

void MetFlowCCIApp::load_indicator(int initial_adaption_counter) {

  std::stringstream filename;
  filename << "PP/indicators.adapt_" << initial_adaption_counter << ".txt";
  std::ifstream indicator_file;
  std::string line;
  vort.clear();
  energy.clear();
  apost.clear();
  ns_apost.clear();
  bool vort_count, en_count, ap_count, ns_count;
  vort_count = false;
  en_count = false;
  ap_count = false;
  ns_count = false;
  indicator_file.open(filename.str().c_str());
  while (std::getline(indicator_file, line)) {
    if (vort_count)
      vort.push_back(std::atof(line.c_str()));
    if (en_count)
      energy.push_back(std::atof(line.c_str()));
    if (ap_count)
      apost.push_back(std::atof(line.c_str()));
    if (ns_count)
      ns_apost.push_back(std::atof(line.c_str()));

    if (line.compare(0, 11, "Vorticity^2") == 0) {
      vort_count = true;
    }
    if (line.compare(0, 6, "Energy") == 0) {
      vort_count = false;
      en_count = true;
    }
    if (line.compare(0, 12, "A posteriori") == 0) {
      en_count = false;
      ap_count = true;
    }
    if (line.compare(0, 15, "NS A posteriori") == 0) {
      ap_count = false;
      ns_count = true;
    }
  }

  vort.erase(vort.end() - 1, vort.end());
  energy.erase(energy.end() - 1, energy.end());
  apost.erase(apost.end() - 1, apost.end());

  indicator_file.close();
}

// write indicator values to file

void MetFlowCCIApp::write_indicators(int adaption_counter) {
  std::ofstream indicator_file;
  std::stringstream filename;

  filename << "PP/indicators.adapt_" << adaption_counter << ".txt";

  indicator_file.open(filename.str().c_str());
  indicator_file << "Vorticity^2 \n";
  for (int i = 0; i < vort.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << vort.at(i) << "\n";
  indicator_file << "Energy \n";
  for (int i = 0; i < energy.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << energy.at(i) << "\n";
  indicator_file << "A posteriori\n";
  for (int i = 0; i < apost.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << apost.at(i) << "\n";
  indicator_file << "NS A posteriori\n";
  for (int i = 0; i < ns_apost.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << ns_apost.at(i) << "\n";

  indicator_file.close();
}

// load temporary indicator values

void MetFlowCCIApp::load_temp_indicator(int adaption_counter, int time) {
  std::stringstream filename;
  filename << "PP/indicators.adapt_" << adaption_counter << ".timestep_" << time
           << ".txt";
  std::ifstream indicator_file;
  std::string line;
  vort.clear();
  energy.clear();
  apost.clear();
  ns_apost.clear();
  bool vort_count, en_count, ap_count, ns_count;
  vort_count = false;
  en_count = false;
  ap_count = false;
  ns_count = false;
  indicator_file.open(filename.str().c_str());
  while (std::getline(indicator_file, line)) {
    if (vort_count)
      vort.push_back(std::atof(line.c_str()));
    if (en_count)
      energy.push_back(std::atof(line.c_str()));
    if (ap_count)
      apost.push_back(std::atof(line.c_str()));
    if (ns_count)
      ns_apost.push_back(std::atof(line.c_str()));

    if (line.compare(0, 11, "Vorticity^2") == 0) {
      vort_count = true;
    }
    if (line.compare(0, 6, "Energy") == 0) {
      vort_count = false;
      en_count = true;
    }
    if (line.compare(0, 12, "A posteriori") == 0) {
      en_count = false;
      ap_count = true;
    }
    if (line.compare(0, 15, "NS A posteriori") == 0) {
      ap_count = false;
      ns_count = true;
    }
  }

  vort.erase(vort.end() - 1, vort.end());
  energy.erase(energy.end() - 1, energy.end());
  apost.erase(apost.end() - 1, apost.end());

  indicator_file.close();
}

// write temporary indicator values

void MetFlowCCIApp::write_temp_indicators(int adaption_counter, int time) {
  std::ofstream indicator_file;
  std::stringstream filename;

  filename << "PP/indicators.adapt_" << adaption_counter << ".timestep_" << time
           << ".txt";

  indicator_file.open(filename.str().c_str());
  indicator_file << "Vorticity^2 \n";
  for (int i = 0; i < vort.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << vort.at(i) << "\n";
  indicator_file << "Energy \n";
  for (int i = 0; i < energy.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << energy.at(i) << "\n";
  indicator_file << "A posteriori\n";
  for (int i = 0; i < apost.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << apost.at(i) << "\n";
  indicator_file << "NS A posteriori\n";
  for (int i = 0; i < ns_apost.size(); ++i)
    indicator_file << std::setfill(' ') << std::setprecision(20)
                   << std::setw(25) << ns_apost.at(i) << "\n";

  indicator_file.close();
}

void MetFlowCCIApp::visualize_solution(int step, int adaption_counter) {
  std::stringstream input;
  if (step < 10)
    input << "CCI.adapt_" << adaption_counter << ".000" << step;
  else if (step < 100)
    input << "CCI.adapt_" << adaption_counter << ".00" << step;
  else if (step < 1000)
    input << "CCI.adapt_" << adaption_counter << ".0" << step;
  else
    input << "CCI.adapt_" << adaption_counter << "." << step;

  if (num_partitions_ > 1)
    input << ".pvtu";
  else
    input << ".vtu";

  std::string visu_filename = "out/" + input.str();

  visu_names_.push_back("u");
  visu_names_.push_back("v");
  visu_names_.push_back("p");

  pp_sol_.UpdateValues();

  CellVisualization visu(space_, 1);
  for (int i = 0; i < space_.get_nb_var(); i++)
    visu.visualize(EvalFeFunction(space_, pp_sol_, i), visu_names_.at(i));
  visu.visualize(EvalFeFunction(space_, visu_sol_vort_, 0), "vorticity");
  visu.write(visu_filename);
}

void MetFlowCCIApp::visualize_solution_dual(int step, int adaption_counter) {
  // Filename
  std::stringstream input;
  if (step < 10)
    input << "CCI.adapt_dual_" << adaption_counter << ".000" << step;
  else if (step < 100)
    input << "CCI.adapt_dual_" << adaption_counter << ".00" << step;
  else if (step < 1000)
    input << "CCI.adapt_dual_" << adaption_counter << ".0" << step;
  else
    input << "CCI.adapt_dual_" << adaption_counter << "." << step;

  if (num_partitions_ > 1)
    input << ".pvtu";
  else
    input << ".vtu";
  std::string visu_filename = "out/" + input.str();

  // Names of the data
  std::vector< std::string > visu_names_dual;
  visu_names_dual.push_back("u");
  visu_names_dual.push_back("v");
  if (DIM == 3)
    visu_names_dual.push_back("w");
  visu_names_dual.push_back("p");

  // Update values in postprocessing vector
  pp_sol_dual_.UpdateValues();

  // Visualize
  CellVisualization visu(space_dual_, 1);
  for (int i = 0; i < space_dual_.get_nb_var(); i++)
    visu.visualize(EvalFeFunction(space_dual_, pp_sol_dual_, i),
                   visu_names_dual.at(i));
  visu.write(visu_filename);
}

void MetFlowCCIApp::visualize_mesh(int adaption_counter) {
  PVtkWriter writer(MPI::COMM_WORLD);
  std::stringstream output_file;
  output_file << "mesh/adapt_mesh." << adaption_counter << ".pvtu";
  writer.add_all_attributes(*mesh_, true);
  writer.write(output_file.str().c_str(), *mesh_);
}

void MetFlowCCIApp::visualize_indicators(int adaption_counter) {
  std::stringstream input;
  input << "PP/indicators.adapt_" << adaption_counter << ".vtu";

  CellVisualization visu(space_, 1);
  visu.visualize_cell_data(vort, "vort");
  visu.visualize_cell_data(energy, "energy");
  visu.visualize_cell_data(apost, "aposteriori");
  visu.visualize_cell_data(ns_apost, "ns_aposteriori");
  visu.write(input.str());
}

void MetFlowCCIApp::save_solution_to_file(LAD::VectorType const &sol,
                                          std::string const &filename) {
  // convert CoupledVector to std::vector

  std::vector< LAD::DataType > backup_vec(sol.size_global(), 1.e20);

  std::vector< int > dof_ids;
  std::vector< double > values;

  PpVector< LAD > pp_sol;
  pp_sol.Init(MPI_COMM_WORLD, space_.dof(), sol);
  pp_sol.InitStructure();
  pp_sol.UpdateValues();
  pp_sol.GetDofsAndValues(dof_ids, values);

  for (int i = 0; i < values.size(); ++i)
    backup_vec.at(dof_ids[i]) = values.at(i);

  // store values to file

  std::ofstream out_file;
  out_file.open(filename.c_str());

  // write number of entries in the very first line
  out_file << values.size() << std::endl;

  // write the data, one value per line
  for (int i = 0; i < backup_vec.size(); ++i)
    out_file << std::setfill(' ') << std::setprecision(20) << std::setw(25)
             << backup_vec.at(i) << "\n";

  out_file.close();

  std::cout << "Writing file '" << filename << "' done." << std::endl;
}

void MetFlowCCIApp::read_solution_from_file(LAD::VectorType &sol,
                                            std::string const &filename) {
  std::string strline;
  double *data;

  // open file

  std::cout << "Reading file '" << filename << "' done." << std::endl;
  std::ifstream file(filename.c_str());
  if (!file) {
    std::cerr << "ERROR: Can't open file : " << filename << std::endl;
    quit_program();
  }

  // resize vector (use dimension given in first line)

  getline(file, strline);
  std::istringstream istl1(strline.c_str());
  int n;
  istl1 >> n;
  assert(n > 0);
  data = new double[n];

  // read vector data

  for (int k = 0; k < n; ++k) {
    getline(file, strline);
    std::istringstream istlc(strline.c_str());

    istlc >> data[k];
  }

  // set data to CoupledVector

  assert(sol.size_global() == n);
  sol.SetValues(data);

  delete data;
}

void MetFlowCCIApp::write_log_file(int adaption_counter) {
  std::ofstream log_file;
  std::stringstream filename;
  filename << "log/log.adapt_" << adaption_counter << ".txt";
  log_file.open(filename.str().c_str());
  log_file << "#Adaption\n" << adaption_counter << "\n";
  log_file << "#Cells\n" << mesh_->num_entities(DIM) << "\n";
  log_file << "#DOFs\n" << space_.dof().ndofs_global() << "\n";
  log_file << "#ErrorKM\n" << compute_position_error() << "\n";
  log_file << "#ErrorNS\n" << compute_global_error(sum_ns_apost) << "\n";
  log_file << "#Error\n" << compute_global_error(sum_apost) << "\n";
  log_file.close();
}

double MetFlowCCIApp::compute_position_error() {
  double x = fabs(storm_coord_one.at(0)) - 1043.678;
  double y = fabs(storm_coord_one.at(1)) - 153.365;
  return pow((pow(x, 2.) + pow(y, 2.)), .5);
}

void MetFlowCCIApp::write_error(int timestep, int adaption_counter) {
  std::ofstream apost_file;
  std::stringstream filename;
  filename << "log/error.adapt_" << adaption_counter << ".txt";
  if (timestep == 0) {
    apost_file.open(filename.str().c_str());
    apost_file
        << "Time step, Global Aposteriori Error, Global NS Aposteriori Error\n";
    apost_file << std::setprecision(5) << timestep << ", " << std::setfill(' ')
               << std::setprecision(10) << compute_global_error(sum_apost)
               << ", " << std::setfill(' ') << std::setprecision(10)
               << compute_global_error(sum_ns_apost) << "\n";
  } else {
    apost_file.open(filename.str().c_str(), std::fstream::app);
    apost_file << std::setprecision(5) << timestep << ", " << std::setfill(' ')
               << std::setprecision(10) << compute_global_error(sum_apost)
               << ", " << std::setfill(' ') << std::setprecision(10)
               << compute_global_error(sum_ns_apost) << "\n";
  }
}

double MetFlowCCIApp::compute_global_error(
    std::vector< double > const &error_indicators) {
  double norm_error_indicators = 0.;
  for (int i = 0; i < error_indicators.size(); ++i)
    norm_error_indicators += pow(error_indicators.at(i), 2);
  norm_error_indicators = sqrt(norm_error_indicators);
  return norm_error_indicators;
}

/// CCI specific tools

void MetFlowCCIApp::compute_storm_track() {
  storm_coord_one = storm_position(storm_coord_one);
  storm_coord_two = storm_position(storm_coord_two);

  std::cout << "-> Storm 1 is at ( " << storm_coord_one.at(0) << ",\t "
            << storm_coord_one.at(1) << " )\n";
  std::cout << "-> Storm 2 is at ( " << storm_coord_two.at(0) << ",\t "
            << storm_coord_two.at(1) << " )\n";
  std::cout << "\n";
}

// find maximal vorticity on dofs

void MetFlowCCIApp::max_vorticity(std::vector< double > old_coord,
                                  double max_distance, double &max_value,
                                  std::vector< double > &coord) {
  const int tdim = vorticity->space()->mesh().tdim();
  max_value = -9999.;
  coord.resize(DIM, 0.);
  for (mesh::EntityIterator it = vorticity->space()->mesh().begin(tdim),
                            end_it = vorticity->space()->mesh().end(tdim);
       it != end_it; ++it) {
    AssemblyAssistant< DIM > aa_;
    Quadrature< LAD::DataType > q;
    q.set_cell_type(3);
    q.set_quadrature_by_order("GaussQuadrilateral", 6);
    aa_.initialize_for_element(Element(*(vorticity->space()), it->index()), q);

    FunctionValues< double > vort_val_;
    vort_val_.clear();
    aa_.evaluate_fe_function(vort_, 0, vort_val_);

    // loop over quadrature points
    for (int q = 0; q < aa_.num_quadrature_points(); ++q) {
      double vort = vort_val_[q];
      double distance = 0.;
      for (int j = 0; j < DIM; j++)
        distance += pow((aa_.x(q)[j] - old_coord[j]), 2);
      distance = pow(distance, 0.5);

      if ((vort_val_[q] > max_value) && (distance < max_distance)) {
        max_value = vort_val_[q];
        for (int i = 0; i < DIM; ++i)
          coord.at(i) = aa_.x(q)[i];
      }
    }
  }
}

// calculate storm position as weighted barycenter of vorticity

std::vector< double >
MetFlowCCIApp::storm_position(std::vector< Coordinate > old_coord) {
  // Find maximal value of vorticity in vicinity of old storm position
  double max_vorticity_ = -9999.;
  std::vector< double > max_coord;
  double max_distance_old_position =
      params_["CCI"]["MaxDistanceStorm"].get< double >();

  // std::cout << "Determine Vorticity Bound for Calculation of Storm Position"
  //	    << std::endl;
  max_vorticity(old_coord, max_distance_old_position, max_vorticity_,
                max_coord);
  // std::cout << "-> Max Vorticity on Quadrature Points: " << max_vorticity_ <<
  // std::endl; std::cout << "-> at: X = " << max_coord.at(0) << " and Y = " <<
  // max_coord.at(1) << "\n";
  double vorticity_bound = max_vorticity_ * 0.5;

  if (max_vorticity_ < vorticity_bound * 1.05)
    vorticity_bound = 0.75 * vorticity_bound;
  // std::cout << "-> Vorticity bound:       " << vorticity_bound << std::endl;

  // Calculate Storm Position
  double upper_vort = -999.;
  double upper_x = 0.;
  double upper_y = 0.;

  double weighted_upper_x = 0.;
  double weighted_upper_y = 0.;
  double weight_upper = 0.;
  double r = params_["CCI"]["RadiusStorm"].get< double >();

  // loop over all cells
  const int tdim = vorticity->space()->mesh().tdim();
  for (mesh::EntityIterator it = vorticity->space()->mesh().begin(tdim),
                            end_it = vorticity->space()->mesh().end(tdim);
       it != end_it; ++it) {

    AssemblyAssistant< DIM > aa_;
    Quadrature< LAD::DataType > q;
    q.set_cell_type(3);
    q.set_quadrature_by_order("GaussQuadrilateral", 6);
    aa_.initialize_for_element(Element(*(vorticity->space()), it->index()), q);

    FunctionValues< double > vort_val_;
    vort_val_.clear();
    aa_.evaluate_fe_function(vort_, 0, vort_val_);

    // loop over quadrature points
    for (int q = 0; q < aa_.num_quadrature_points(); ++q) {
      double vort = vort_val_[q];
      double x_tilde = aa_.x(q)[0];
      double y_tilde = aa_.x(q)[1];
      std::vector< double > temp_coord(DIM, 0.);

      double distance = 0.;
      for (int j = 0; j < DIM; j++)
        distance += pow((aa_.x(q)[j] - max_coord.at(j)), 2);

      distance = pow(distance, 0.5);

      if ((vort > vorticity_bound) && (distance < r)) {
        double weight = vort - vorticity_bound;
        weight = weight * (aa_.detJ(q) * aa_.detJ(q));
        weighted_upper_x += x_tilde * weight;
        weighted_upper_y += y_tilde * weight;
        weight_upper += weight;
      }
    }
  }
  weighted_upper_x = weighted_upper_x / weight_upper;
  weighted_upper_y = weighted_upper_y / weight_upper;

  weighted_upper_x = weighted_upper_x;
  weighted_upper_y = weighted_upper_y;

  std::vector< Coordinate > coord(DIM, 0);
  coord[0] = weighted_upper_x;
  coord[1] = weighted_upper_y;

  return coord;
}

void MetFlowCCIApp::prepare_storm_track(int adaption_counter) {
  // Prepare storm tracking
  storm_coord_one.clear();
  storm_coord_two.clear();
  storm_coord_one.resize(DIM, 0.);
  storm_coord_one.at(0) = -distance_ / 2.;
  storm_coord_two.resize(DIM, 0.);
  storm_coord_two.at(0) = distance_ / 2.;

  std::ofstream track_file;
  std::stringstream filename;
  filename << "PP/stormtrack.adapt_" << adaption_counter << ".txt";
  track_file.open(filename.str().c_str());
  track_file << "Timestep,\tX-Coord 1\tY-Coord 1,\tX-Coord 2,\tY-Coord 2\n";
  track_file.close();
}

void MetFlowCCIApp::write_storm_track(int time_step, int adaption_counter) {
  std::ofstream track_file;
  std::stringstream filename;
  filename << "PP/stormtrack.adapt_" << adaption_counter << ".txt";
  track_file.open(filename.str().c_str(), std::ios::app);
  track_file << time_step << ",\t" << storm_coord_one.at(0) << ",\t"
             << storm_coord_one.at(1) << ",\t" << storm_coord_two.at(0) << ",\t"
             << storm_coord_two.at(1) << "\n";
  track_file.close();
}
