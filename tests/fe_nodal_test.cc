// Copyright (C) 2011-2017 Vincent Heuveline
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

/// \author Philipp Gerstner
/// TODO: uncomment missing elements

#define BOOST_TEST_MODULE fe_transformation

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include <mpi.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "common/macros.h"
#include "fem/reference_cell.h"
#include "fem/fe_reference.h"
#include "fem/fe_transformation.h"

#include "fem/ansatz/ansatz_p_line_lagrange.h"
#include "fem/ansatz/ansatz_p_tri_lagrange.h"
#include "fem/ansatz/ansatz_aug_p_tri_mono.h"
#include "fem/ansatz/ansatz_sum.h"
#include "fem/ansatz/ansatz_p_tet_lagrange.h"
#include "fem/ansatz/ansatz_q_quad_lagrange.h"
#include "fem/ansatz/ansatz_q_hex_lagrange.h"
//#include "fem/ansatz/ansatz_pyr_lagrange.h"

#include "dof/dof_impl/dof_container_lagrange.h"
#include "dof/dof_impl/dof_container_rt_bdm.h"


#include "test.h"

using namespace hiflow::doffem;

BOOST_AUTO_TEST_CASE(fe_nodal) {

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,  &rank);
  
  if (rank >= 1)
  {
    return;
  }
  
  /// Test if shape functions satisfy nodal property, 

  /// create reference cells 
  ConstRefCellPtr<double, 1> ref_cell_line = ConstRefCellPtr<double, 1> ( new RefCellLineStd<double, 1>);
  ConstRefCellPtr<double, 2> ref_cell_tri = ConstRefCellPtr<double, 2> ( new RefCellTriStd <double, 2>);
  ConstRefCellPtr<double, 2> ref_cell_quad = ConstRefCellPtr<double, 2> ( new RefCellQuadStd<double, 2>);
  ConstRefCellPtr<double, 3> ref_cell_tet = ConstRefCellPtr<double, 3> ( new RefCellTetStd <double, 3>);
  ConstRefCellPtr<double, 3> ref_cell_hex = ConstRefCellPtr<double, 3> ( new RefCellHexStd <double, 3>);
  ConstRefCellPtr<double, 3> ref_cell_pyr = ConstRefCellPtr<double, 3> ( new RefCellPyrStd <double, 3>);

  /// create Lagrange ansatz spaces
  PLineLag<double, 1> ansatz_line_lag(ref_cell_line);
  PTriLag <double, 2> ansatz_tri_lag (ref_cell_tri);
  PTetLag <double, 3> ansatz_tet_lag (ref_cell_tet);
  QQuadLag<double, 2> ansatz_quad_lag(ref_cell_quad);
  QHexLag <double, 3> ansatz_hex_lag (ref_cell_hex);
//PyrLag  <double, 3> ansatz_pyr_lag(ref_cell_pyr);

  /// create Lagrange dof containers
  DofContainerLagrange<double, 1> dof_line_lag(ref_cell_line);
  DofContainerLagrange<double, 2> dof_tri_lag (ref_cell_tri);
  DofContainerLagrange<double, 3> dof_tet_lag (ref_cell_tet);
  DofContainerLagrange<double, 2> dof_quad_lag(ref_cell_quad);
  DofContainerLagrange<double, 3> dof_hex_lag (ref_cell_hex);
//DofContainerLagrange<double, 3> dof_pyr_lag (ref_cell_pyr);
  
  /// create ansatz spaces for RT and BDM elements
  PTriLag <double, 2> ansatz_tri_lag_2 (ref_cell_tri);
  
  PTriLag<double, 2> ansatz_tri_rt_1 (ref_cell_tri);
  AugPTriMono<double, 2> ansatz_tri_rt_2 (ref_cell_tri);
  AnsatzSpaceSum<double, 2> ansatz_tri_rt (ref_cell_tri);
      
  /// create BDM dof containers
  DofContainerRTBDM<double, 2> dof_tri_bdm (ref_cell_tri);
  DofContainerRTBDM<double, 2> dof_tri_rt (ref_cell_tri);
  
  const int min_bdm_deg = 1;
  const int max_bdm_deg = 2;  // TODO: this should be 3
  
  const int min_rt_deg = 0;
  const int max_rt_deg = 2;  // TODO: this should be 3
  
  const size_t lag_nb_comp = 2;
   
  for (int deg = 0; deg < 4; ++deg) 
  {
    CONSOLE_OUTPUT(rank, " ===================== ");

#if 1
    bool modal_basis = false;

    /// Lagrange line
    ansatz_line_lag.init(deg, lag_nb_comp);
    dof_line_lag.init(deg, lag_nb_comp);
    CONSOLE_OUTPUT(rank, "initialize Lag(" << deg << ") line element with " << dof_line_lag.nb_dof_on_cell() << " dofs and " 
              << ansatz_line_lag.dim() << " ansatz functions ");
    RefElement<double, 1> fe_line_lag;
    FETransformationStandard<double, 1> fe_trafo_line_lag(&fe_line_lag);
    fe_line_lag.init(&ansatz_line_lag,&dof_line_lag,&fe_trafo_line_lag,modal_basis, FE_TYPE_LAGRANGE);
    
    /// Lagrange triangle
    ansatz_tri_lag.init(deg, lag_nb_comp);
    dof_tri_lag.init(deg, lag_nb_comp);
    CONSOLE_OUTPUT(rank, "initialize Lag(" << deg << ") tri  element with " << dof_tri_lag.nb_dof_on_cell() << " dofs and " 
              << ansatz_tri_lag.dim() << " ansatz functions ");
    RefElement<double, 2> fe_tri_lag; 
    FETransformationStandard<double, 2> fe_trafo_tri_lag (&fe_tri_lag);
    fe_tri_lag.init (&ansatz_tri_lag, &dof_tri_lag, &fe_trafo_tri_lag, modal_basis, FE_TYPE_LAGRANGE);
    
    /// Lagrange tet
    ansatz_tet_lag.init(deg, lag_nb_comp);
    dof_tet_lag.init(deg, lag_nb_comp);
    std::cout << "initialize Lag(" << deg << ") tet  element with " << dof_tet_lag.nb_dof_on_cell() << " dofs and " 
              << ansatz_tet_lag.dim() << " ansatz functions " << std::endl;
    RefElement<double, 3> fe_tet_lag;
    FETransformationStandard<double, 3> fe_trafo_tet_lag (&fe_tet_lag);
    fe_tet_lag.init (&ansatz_tet_lag, &dof_tet_lag, &fe_trafo_tet_lag, modal_basis, FE_TYPE_LAGRANGE);

    /// Lagrange quad
    ansatz_quad_lag.init(deg, lag_nb_comp);
    dof_quad_lag.init(deg, lag_nb_comp);
    std::cout << "initialize Lag(" << deg << ") quad element with " << dof_quad_lag.nb_dof_on_cell() << " dofs and " 
              << ansatz_quad_lag.dim() << " ansatz functions " << std::endl;
    RefElement<double, 2> fe_quad_lag;
    FETransformationStandard<double, 2> fe_trafo_quad_lag(&fe_quad_lag);
    fe_quad_lag.init(&ansatz_quad_lag,&dof_quad_lag,&fe_trafo_quad_lag,modal_basis, FE_TYPE_LAGRANGE);

    ///  Lagrange hex
    ansatz_hex_lag.init(deg, lag_nb_comp);
    dof_hex_lag.init(deg, lag_nb_comp);
    CONSOLE_OUTPUT(rank, "initialize Lag(" << deg << ") hex  element with " << dof_hex_lag.nb_dof_on_cell() << " dofs and " 
              << ansatz_hex_lag.dim() << " ansatz functions ");
    RefElement<double, 3> fe_hex_lag;
    FETransformationStandard<double, 3> fe_trafo_hex_lag (&fe_hex_lag);
    fe_hex_lag.init (&ansatz_hex_lag, &dof_hex_lag, &fe_trafo_hex_lag, modal_basis, FE_TYPE_LAGRANGE);

    /// Lagrange pyr
//  ansatz_pyr_lag.init(deg, lag_nb_comp);
//  dof_pyr_lag.init(deg, lag_nb_comp);
//  std::cout << "initialize Lag(" << deg << ") pyramid element with " << dof_pyr_lag.nb_dof_on_cell() << " dofs and " 
//            << ansatz_pyr_lag.dim() << " ansatz functions " << std::endl;
//  RefElement<double, 3> fe_pyr_lag;
//  FETransformationStandard<double, 3> fe_trafo_pyr_lag (&fe_pyr_lag);
//  fe_pyr_lag.init (&ansatz_pyr_lag, &dof_pyr_lag, &fe_trafo_pyr_lag, modal_basis, FE_TYPE_LAGRANGE);
#endif

    /// BDM
    RefElement<double, 2> fe_tri_bdm;
    FETransformationContraPiola<double, 2> fe_trafo_tri_bdm (&fe_tri_bdm);
    if (deg >= min_bdm_deg && deg <= max_bdm_deg)
    {
      ansatz_tri_lag_2.init(deg, 2);
      dof_tri_bdm.init(deg, DOF_CONTAINER_BDM);
      CONSOLE_OUTPUT(rank, "initialize BDM(" << deg << ") tri  element with " << dof_tri_bdm.nb_dof_on_cell() << " dofs and " 
                << ansatz_tri_lag_2.dim() << " ansatz functions ");
      fe_tri_bdm.init (&ansatz_tri_lag_2, &dof_tri_bdm, &fe_trafo_tri_bdm, false, FE_TYPE_BDM);
    }

    /// RT
    RefElement<double, 2> fe_tri_rt;
    FETransformationContraPiola<double, 2> fe_trafo_tri_rt (&fe_tri_rt);
    if (deg >= min_rt_deg && deg <= max_rt_deg)
    {
      ansatz_tri_rt_1.init(deg,2);
      ansatz_tri_rt_2.init(deg);
      ansatz_tri_rt.init(&ansatz_tri_rt_1, &ansatz_tri_rt_2, ANSATZ_RT);
      dof_tri_rt.init(deg, DOF_CONTAINER_RT);  
      CONSOLE_OUTPUT(rank, "initialize RT(" << deg << ")  tri  element with " << dof_tri_rt.nb_dof_on_cell() << " dofs and " 
              << ansatz_tri_rt.dim() << " ansatz functions ");

      fe_tri_rt.init (&ansatz_tri_rt, &dof_tri_rt, &fe_trafo_tri_rt, false, FE_TYPE_RT);
      BOOST_TEST(fe_tri_rt.dim() == (deg+1)*(deg+3));
      BOOST_TEST(dof_tri_rt.nb_dof_on_cell() == ansatz_tri_rt.dim());
    }
    
#if 1
    /// Test Lagrange HEXAHEDRON
    std::vector< DofID > dof_ids_hex_lag  (dof_hex_lag.nb_dof_on_cell(),0);
    for (size_t i = 0; i < dof_hex_lag.nb_dof_on_cell(); ++i)
    {
      dof_ids_hex_lag[i] = i;
    }
    std::vector< std::vector<double> > dof_val_hex_lag;
    dof_hex_lag.evaluate (&fe_hex_lag, dof_ids_hex_lag, dof_val_hex_lag);
    
    for (size_t i = 0; i < dof_val_hex_lag.size(); ++i)
    {
      for (size_t j = 0; j < dof_val_hex_lag[i].size(); ++j)
      {
        double val = dof_val_hex_lag[i][j];
        if (i != j)
        {
          BOOST_TEST(std::abs(val) < 1.e-10);
        }
        else
        {
          BOOST_TEST(std::abs(val - 1.0) < 1.e-10);
        }
      }
    }
    CONSOLE_OUTPUT(rank, "hex  lag(" << deg << ") test: passed ");

    /// Test Lagrange TETRAHEDRON
    std::vector< DofID > dof_ids_tet_lag  (dof_tet_lag.nb_dof_on_cell(),0);
    for (size_t i = 0; i < dof_tet_lag.nb_dof_on_cell(); ++i)
    {
      dof_ids_tet_lag[i] = i;
    }
    std::vector< std::vector<double> > dof_val_tet_lag;
    dof_tet_lag.evaluate (&fe_tet_lag, dof_ids_tet_lag, dof_val_tet_lag);
    
    for (size_t i = 0; i < dof_val_tet_lag.size(); ++i)
    {
      for (size_t j = 0; j < dof_val_tet_lag[i].size(); ++j)
      {
        double val = dof_val_tet_lag[i][j];
        if (i != j)
        {
          BOOST_TEST(std::abs(val) < 1.e-10);
        }
        else
        {
          BOOST_TEST(std::abs(val - 1.0)< 1.e-10);
        }
      }
    }
    std::cout << "tet  lag(" << deg << ") test: passed " << std::endl;

    /// Test Lagrange TRIANGLE
    std::vector< DofID > dof_ids_tri_lag  (dof_tri_lag.nb_dof_on_cell(),0);
    for (size_t i = 0; i < dof_tri_lag.nb_dof_on_cell(); ++i)
    {
      dof_ids_tri_lag[i] = i;
    }
    std::vector< std::vector<double> > dof_val_tri_lag;
    dof_tri_lag.evaluate (&fe_tri_lag, dof_ids_tri_lag, dof_val_tri_lag);
    
    for (size_t i = 0; i < dof_val_tri_lag.size(); ++i)
    {
      for (size_t j = 0; j < dof_val_tri_lag[i].size(); ++j)
      {
        double val = dof_val_tri_lag[i][j];
        if (i != j)
        {
          BOOST_TEST(std::abs(val) < 1.e-10);
        }
        else
        {
          BOOST_TEST(std::abs(val - 1.0)< 1.e-10);
        }
      }
    }
    CONSOLE_OUTPUT(rank, "tri  lag(" << deg << ") test: passed ");

    /// Test Lagrange QUADRILATERAL

    std::vector< DofID > dof_ids_quad_lag (dof_quad_lag.nb_dof_on_cell(),0);
    for (size_t i = 0; i < dof_quad_lag.nb_dof_on_cell(); ++i)
    {
      dof_ids_quad_lag[i] = i;
    }
    std::vector< std::vector<double> > dof_val_quad_lag;
    dof_quad_lag.evaluate (&fe_quad_lag, dof_ids_quad_lag, dof_val_quad_lag);
    
    for (size_t i = 0; i < dof_val_quad_lag.size(); ++i)
    {
      for (size_t j = 0; j < dof_val_quad_lag[i].size(); ++j)
      {
        double val = dof_val_quad_lag[i][j];
        if (i != j)
        {
          BOOST_TEST(std::abs(val) < 1.e-10);
        }
        else
        {
          BOOST_TEST(std::abs(val - 1.0)< 1.e-10);
        }
      }
    }
    std::cout << "quad lag(" << deg << ") test: passed " << std::endl;
 
    /// Test Lagrange LINE
    std::vector< std::vector<double> > dof_val_line_lag;
    std::vector< DofID > dof_ids_line_lag (dof_line_lag.nb_dof_on_cell(),0);
    for (size_t i = 0; i < dof_line_lag.nb_dof_on_cell(); ++i)
    {
      dof_ids_line_lag[i] = i;
    }
    dof_line_lag.evaluate (&fe_line_lag, dof_ids_line_lag, dof_val_line_lag);
    
    for (size_t i = 0; i < dof_val_line_lag.size(); ++i)
    {
      for (size_t j = 0; j < dof_val_line_lag[i].size(); ++j)
      {
        double val = dof_val_line_lag[i][j];
        if (i != j)
        {
          BOOST_TEST(std::abs(val) < 1.e-10);
        }
        else
        {
          BOOST_TEST(std::abs(val - 1.0) < 1.e-10);
        }
      }
    }
    CONSOLE_OUTPUT(rank, "line lag(" << deg << ") test: passed ");


    /// Test Lagrange PYRAMID
/*
    std::vector< DofID > dof_ids_pyr_lag  (dof_pyr_lag.nb_dof_on_cell(),0);
    for (size_t i = 0; i < dof_pyr_lag.nb_dof_on_cell(); ++i)
    {
      dof_ids_pyr_lag[i] = i;
    }
    std::vector< std::vector<double> > dof_val_pyr_lag;
    dof_pyr_lag.evaluate (&fe_pyr_lag, dof_ids_pyr_lag, dof_val_pyr_lag);
    
    for (size_t i = 0; i < dof_val_pyr_lag.size(); ++i)
    {
      for (size_t j = 0; j < dof_val_pyr_lag[i].size(); ++j)
      {
        double val = dof_val_pyr_lag[i][j];
        if (i != j)
        {
          TEST_EQUAL_EPS(val, 0.0, 1.e-10);
        }
        else
        {
          TEST_EQUAL_EPS(val, 1.0, 1.e-10);
        }
      }
    }
    std::cout << "pyr  lag(" << deg << ") test: passed " << std::endl;
*/
#endif
    
  
    /// Test BDM TRIANGLE
    if (deg >= min_bdm_deg && deg <= max_bdm_deg)
    {
      std::vector< DofID > dof_ids_tri_bdm;
      dof_ids_tri_bdm.resize(dof_tri_bdm.nb_dof_on_cell(), 0);
      for (size_t i = 0; i < dof_tri_bdm.nb_dof_on_cell(); ++i)
      {
        dof_ids_tri_bdm[i] = i;
      }
      std::vector< std::vector<double> > dof_val_tri_bdm;
      dof_tri_bdm.evaluate (&fe_tri_bdm, dof_ids_tri_bdm, dof_val_tri_bdm);
      
      for (size_t i = 0; i < dof_val_tri_bdm.size(); ++i)
      {
        for (size_t j = 0; j < dof_val_tri_bdm[i].size(); ++j)
        {
          double val = dof_val_tri_bdm[i][j];
          if (i != j)
          {
            BOOST_TEST(std::abs(val) < 1.e-10);
          }
          else
          {
            BOOST_TEST(std::abs(val - 1.0) < 1.e-10);
          }
        }
      }
      CONSOLE_OUTPUT(rank, "tri  BDM(" << deg << ") test: passed ");
    }
    
    /// Test RT TRIANGLE
    if (deg >= min_rt_deg && deg <= max_rt_deg)
    {
      std::vector< DofID > dof_ids_tri_rt;
      dof_ids_tri_rt.resize(dof_tri_rt.nb_dof_on_cell(), 0);
      for (size_t i = 0; i < dof_tri_rt.nb_dof_on_cell(); ++i)
      {
        dof_ids_tri_rt[i] = i;
      }
      std::vector< std::vector<double> > dof_val_tri_rt;
      dof_tri_rt.evaluate (&fe_tri_rt, dof_ids_tri_rt, dof_val_tri_rt);
        
      for (size_t i = 0; i < dof_val_tri_rt.size(); ++i)
      {
        for (size_t j = 0; j < dof_val_tri_rt[i].size(); ++j)
        {
          double val = dof_val_tri_rt[i][j];
          if (i != j)
          {
            BOOST_TEST(std::abs(val) < 1.e-10);
          }
          else
          {
            BOOST_TEST(std::abs(val - 1.0) < 1.e-10);
          }
        }
      }
      CONSOLE_OUTPUT(rank, "tri  RT(" << deg << ")  test: passed ");
    }
  }
  MPI_Finalize();
  return;
}
