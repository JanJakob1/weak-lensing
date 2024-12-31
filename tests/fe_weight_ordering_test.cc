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

/// \author Michael Schick, Martin Baumann

#define BOOST_TEST_MODULE fe_weight_ordering

#include <mpi.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
namespace utf = boost::unit_test;

#include <cmath>
#include <iostream>
#include <vector>

#include "common/macros.h"

#include "fem/reference_cell.h"
#include "fem/fe_reference.h"
#include "fem/fe_transformation.h"

#include "fem/ansatz/ansatz_p_line_lagrange.h"
#include "fem/ansatz/ansatz_p_tri_lagrange.h"
#include "fem/ansatz/ansatz_p_tet_lagrange.h"
#include "fem/ansatz/ansatz_q_quad_lagrange.h"
#include "fem/ansatz/ansatz_q_hex_lagrange.h"
//#include "fem/ansatz/ansatz_pyr_lagrange.h"

#include "dof/dof_impl/dof_container_lagrange.h"
#include "dof/dof_impl/dof_functional.h"
#include "dof/dof_impl/dof_functional_point.h"

#include "test.h"

using namespace hiflow::doffem;

BOOST_AUTO_TEST_CASE(fe_weight_ordering) {

  int rank = 0;
  
  // Test to check if the numbering in the fem module is consistent
  // i.e. check if weights are equal to 1.0 on the desired dof

  // create reference cells 
  ConstRefCellPtr<double, 1> ref_cell_line = ConstRefCellPtr<double, 1> ( new RefCellLineStd<double, 1>);
  ConstRefCellPtr<double, 2> ref_cell_tri = ConstRefCellPtr<double, 2> ( new RefCellTriStd <double, 2>);
  ConstRefCellPtr<double, 2> ref_cell_quad = ConstRefCellPtr<double, 2> ( new RefCellQuadStd<double, 2>);
  ConstRefCellPtr<double, 3> ref_cell_tet = ConstRefCellPtr<double, 3> ( new RefCellTetStd <double, 3>);
  ConstRefCellPtr<double, 3> ref_cell_hex = ConstRefCellPtr<double, 3> ( new RefCellHexStd <double, 3>);
  ConstRefCellPtr<double, 3> ref_cell_pyr = ConstRefCellPtr<double, 3> ( new RefCellPyrStd <double, 3>);

  // create ansatz spaces
  PLineLag<double, 1> ansatz_line_lag(ref_cell_line);
  PTriLag <double, 2> ansatz_tri_lag (ref_cell_tri);
  PTetLag <double, 3> ansatz_tet_lag (ref_cell_tet);
  QQuadLag<double, 2> ansatz_quad_lag(ref_cell_quad);
  QHexLag <double, 3> ansatz_hex_lag (ref_cell_hex);
//PyrLag  <double, 3> ansatz_pyr_lag(ref_cell_pyr);

  // create dof containers
  DofContainerLagrange<double, 1> dof_line_lag(ref_cell_line);
  DofContainerLagrange<double, 2> dof_tri_lag (ref_cell_tri);
  DofContainerLagrange<double, 3> dof_tet_lag (ref_cell_tet);
  DofContainerLagrange<double, 2> dof_quad_lag(ref_cell_quad);
  DofContainerLagrange<double, 3> dof_hex_lag (ref_cell_hex);
//DofContainerLagrange<double, 3> dof_pyr_lag (ref_cell_pyr);


  for (int deg = 0; deg < 9; ++deg) {
    ansatz_line_lag.init(deg, 1);
    ansatz_tri_lag.init(deg, 1);
    ansatz_tet_lag.init(deg, 1);
    ansatz_quad_lag.init(deg, 1);
    ansatz_hex_lag.init(deg, 1);
//  ansatz_pyr_lag.init(deg, 1);
    
    dof_line_lag.init(deg, 1);
    dof_tri_lag.init(deg, 1);
    dof_tet_lag.init(deg, 1);
    dof_quad_lag.init(deg, 1);
    dof_hex_lag.init(deg, 1);
//  dof_pyr_lag.init(deg, 1);
    
    RefElement<double, 1> fe_line_lag;
    RefElement<double, 2> fe_tri_lag;
    RefElement<double, 3> fe_tet_lag;
    RefElement<double, 2> fe_quad_lag;
    RefElement<double, 3> fe_hex_lag;
//  RefElement<double, 3> fe_pyr_lag;
    
    FETransformationStandard<double, 1> fe_trafo_line_lag(&fe_line_lag);
    FETransformationStandard<double, 2> fe_trafo_tri_lag(&fe_tri_lag);
    FETransformationStandard<double, 3> fe_trafo_tet_lag(&fe_tet_lag);
    FETransformationStandard<double, 2> fe_trafo_quad_lag(&fe_quad_lag);
    FETransformationStandard<double, 3> fe_trafo_hex_lag(&fe_hex_lag);
//  FETransformationStandard<double, 3> fe_trafo_pyr_lag(&fe_pyr_lag);
    
    bool modal_basis = false;
    fe_line_lag.init(&ansatz_line_lag,&dof_line_lag,&fe_trafo_line_lag,modal_basis, FE_TYPE_LAGRANGE);
    fe_tri_lag.init (&ansatz_tri_lag, &dof_tri_lag, &fe_trafo_tri_lag, modal_basis, FE_TYPE_LAGRANGE);
    fe_tet_lag.init (&ansatz_tet_lag, &dof_tet_lag, &fe_trafo_tet_lag, modal_basis, FE_TYPE_LAGRANGE);
    fe_quad_lag.init(&ansatz_quad_lag,&dof_quad_lag,&fe_trafo_quad_lag,modal_basis, FE_TYPE_LAGRANGE);
    fe_hex_lag.init (&ansatz_hex_lag, &dof_hex_lag, &fe_trafo_hex_lag, modal_basis, FE_TYPE_LAGRANGE);
//  fe_pyr_lag.init (&ansatz_pyr_lag, &dof_pyr_lag, &fe_trafo_pyr_lag, modal_basis, FE_TYPE_LAGRANGE);

    for (int ind_dof = 0; ind_dof < fe_line_lag.nb_dof_on_cell(); ++ind_dof) 
    {
      std::vector< double > weight(fe_line_lag.nb_dof_on_cell());
      DofPointEvaluation<double, 1> const * lag_dof = dynamic_cast<DofPointEvaluation<double, 1> const *>(fe_line_lag.dof_container()->get_dof(ind_dof));
      hiflow::Vec<1,double> dof_coord = lag_dof->get_point();
      fe_line_lag.N(dof_coord, weight);

      for (int i = 0; i < static_cast< int >(weight.size()); ++i)
        if (std::abs(weight[i] - 1.0) < 1.0e-15) {
          if (ind_dof != i)
            CONSOLE_OUTPUT(rank, "========== TEST FAILED FOR LINE! with Degree "
                      << deg << " ==========");

          BOOST_CHECK_EQUAL(ind_dof, i);
        }
    }
    CONSOLE_OUTPUT(rank, "Line Lagrange(" << deg << ") : passed ");
    
    for (int ind_dof = 0; ind_dof < fe_tri_lag.nb_dof_on_cell(); ++ind_dof) 
    {
      std::vector< double > weight(fe_tri_lag.nb_dof_on_cell());
      DofPointEvaluation<double, 2> const * lag_dof = dynamic_cast<DofPointEvaluation<double, 2> const *>(fe_tri_lag.dof_container()->get_dof(ind_dof));
      hiflow::Vec<2,double> dof_coord = lag_dof->get_point();
      fe_tri_lag.N(dof_coord, weight);

      for (int i = 0; i < static_cast< int >(weight.size()); ++i)
        if (std::abs(weight[i] - 1.0) < 1.0e-15) {
          if (ind_dof != i)
            CONSOLE_OUTPUT(rank, "========== TEST FAILED FOR TRIANGLE! with Degree "
                      << deg << " ==========");

          BOOST_TEST(ind_dof == i);
        }
    }
    CONSOLE_OUTPUT(rank, "Tri  Lagrange(" << deg << ") : passed ");
    

    for (int ind_dof = 0; ind_dof < fe_tet_lag.nb_dof_on_cell(); ++ind_dof) 
    {
      std::vector< double > weight(fe_tet_lag.nb_dof_on_cell());
      DofPointEvaluation<double, 3> const * lag_dof = dynamic_cast<DofPointEvaluation<double, 3> const *>(fe_tet_lag.dof_container()->get_dof(ind_dof));
      hiflow::Vec<3,double> dof_coord = lag_dof->get_point();
      fe_tet_lag.N(dof_coord, weight);

      for (int i = 0; i < static_cast< int >(weight.size()); ++i)
        if (std::abs(weight[i] - 1.0) < 1.0e-15) {
          if (ind_dof != i)
            std::cout << "========== TEST FAILED FOR TETRAHEDRON! with Degree "
                      << deg << " ==========" << std::endl;

          BOOST_CHECK_EQUAL(ind_dof, i);
        }
    }
    std::cout << "Tet  Lagrange(" << deg << ") : passed " << std::endl;
    
    for (int ind_dof = 0; ind_dof < fe_quad_lag.nb_dof_on_cell(); ++ind_dof) 
    {
      std::vector< double > weight(fe_quad_lag.nb_dof_on_cell());
      DofPointEvaluation<double, 2> const * lag_dof = dynamic_cast<DofPointEvaluation<double, 2> const *>(fe_quad_lag.dof_container()->get_dof(ind_dof));
      hiflow::Vec<2,double> dof_coord = lag_dof->get_point();
      fe_quad_lag.N(dof_coord, weight);

      for (int i = 0; i < static_cast< int >(weight.size()); ++i)
        if (std::abs(weight[i] - 1.0) < 1.0e-15) {
          if (ind_dof != i)
            std::cout << "========== TEST FAILED FOR QUADRILATERAL! with Degree "
                      << deg << " ==========" << std::endl;

          BOOST_CHECK_EQUAL(ind_dof, i);
        }
    }
    std::cout << "Quad Lagrange(" << deg << ") : passed " << std::endl;

    for (int ind_dof = 0; ind_dof < fe_hex_lag.nb_dof_on_cell(); ++ind_dof) 
    {
      std::vector< double > weight(fe_hex_lag.nb_dof_on_cell());
      DofPointEvaluation<double, 3> const * lag_dof = dynamic_cast<DofPointEvaluation<double, 3> const *>(fe_hex_lag.dof_container()->get_dof(ind_dof));
      hiflow::Vec<3,double> dof_coord = lag_dof->get_point();
      fe_hex_lag.N(dof_coord, weight);

      for (int i = 0; i < static_cast< int >(weight.size()); ++i)
        if (std::abs(weight[i] - 1.0) < 1.0e-15) {
          if (ind_dof != i)
            CONSOLE_OUTPUT(rank, "========== TEST FAILED FOR HEXAHEDRON! with Degree "
                      << deg << " ==========");

          BOOST_CHECK_EQUAL(ind_dof, i);
        }
    }
    CONSOLE_OUTPUT(rank, "Hex  Lagrange(" << deg << ") : passed "); 
    
/*
    for (int ind_dof = 0; ind_dof < fe_pyr_lag.nb_dof_on_cell(); ++ind_dof) 
    {
      std::vector< double > weight(fe_pyr_lag.nb_dof_on_cell());
      DofPointEvaluation<double, 3> const * lag_dof = dynamic_cast<DofPointEvaluation<double, 3> const *>(fe_pyr_lag.dof_container()->get_dof(ind_dof));
      hiflow::Vec<3,double> dof_coord = lag_dof->get_point();
      fe_pyr_lag.N(dof_coord, weight);

      for (int i = 0; i < static_cast< int >(weight.size()); ++i)
        if (std::abs(weight[i] - 1.0) < 1.0e-15) {
          if (ind_dof != i)
            std::cout << "========== TEST FAILED FOR PYRAMID! with Degree "
                      << deg << " ==========" << std::endl;

          BOOST_CHECK_EQUAL(ind_dof, i);
        }
    }
    std::cout << "Pyr  Lagrange(" << deg << ") : passed " << std::endl; 
*/
  }
}
