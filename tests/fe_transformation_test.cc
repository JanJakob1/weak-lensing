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

#define BOOST_TEST_MODULE fe_transformation

#include <mpi.h>
#include <string>
#include <vector>
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "hiflow.h"
#include "test.h"

using namespace std;
using namespace hiflow;
using namespace hiflow::mesh;

/// FE tranformation test
///
/// \brief For a list of files it is analysed, whether the FE transformation applied 
/// to its inverse yields the identity operator
///


static const char *datadir = MESH_DATADIR;

template <class DataType, int DIM>
void check_test_points (std::vector<Vec<DIM,DataType> >& test_points, 
                        RefElement<DataType, DIM> * ref_fe,  
                        MappingPhys2Ref <DataType, DIM, MappingRef2Phys<DataType, DIM, RefElement<DataType, DIM> > > * eval_phys_on_ref)
{       
  for (int q=0; q<test_points.size(); ++q)
  {
    Vec<DIM,DataType> pt = test_points[q];
    size_t weight_size = ref_fe->dim() * ref_fe->nb_comp();
                        
    std::vector< DataType > weight_unmapped (weight_size, 0.);
    std::vector< DataType > weight_mapped (weight_size, 0.);
               
    ref_fe->N(pt, weight_unmapped);
    eval_phys_on_ref->evaluate(pt, weight_mapped);
/*            
std::cout << " ========================= " << std::endl;
std::cout << "point " << pt[0] << " , " << pt[1] << std::endl;
std::cout << std::scientific << std::setprecision(2);

    for (size_t i=0; i<weight_size; ++i)
    {
std::cout << std::setw(5) << weight_unmapped[i] << " ";
    }
std::cout << std::endl;
    for (size_t i=0; i<weight_size; ++i)
    {
std::cout << std::setw(5) << weight_mapped[i] << " ";
    }
std::cout << std::endl;
*/
    for (size_t i=0; i<weight_size; ++i)
    {
      BOOST_TEST(std::abs(weight_unmapped[i]- weight_mapped[i])<1e-10);
    }
  }
}

template <class DataType, int DIM>
void check_dofs (DofContainer<DataType, DIM> const * dof, 
                 RefElement<DataType, DIM> * ref_fe,  
                 MappingPhys2Ref <DataType, DIM, MappingRef2Phys<DataType, DIM, RefElement<DataType, DIM> > > * eval_phys_on_ref)
{
  size_t nb_dof = ref_fe->nb_dof_on_cell();
  std::vector<DofID> all_dofs(nb_dof);
  for (size_t k=0; k<nb_dof; ++k)
  {
    all_dofs[k] = k;
  }
            
  // evaluate all dofs for unmapped shape functions
  std::vector< std::vector<DataType> > dof_vals_unmapped;
  dof->evaluate (ref_fe, all_dofs, dof_vals_unmapped);
                        
  // evaluate all dofs at inversly mapped shape functions
  // -> should yield identity matrix 
  std::vector< std::vector<DataType> > dof_vals_mapped;
  dof-> evaluate (eval_phys_on_ref, all_dofs, dof_vals_mapped);
        
  BOOST_TEST(dof_vals_unmapped.size() == dof_vals_mapped.size()); 
            
//std::cout << std::scientific << std::setprecision(2);

  for (size_t l=0; l<dof_vals_unmapped.size(); ++l)
  {
//    std::cout << " dof " << l << std::endl;
    BOOST_TEST(dof_vals_unmapped[l].size() == dof_vals_mapped[l].size());
/*
    for (size_t i=0; i<dof_vals_unmapped[l].size(); ++i)
    {
std::cout << std::setw(5) << dof_vals_unmapped[l][i] << " ";
    }
    std::cout << std::endl;
    for (size_t i=0; i<dof_vals_unmapped[l].size(); ++i)
    {
std::cout << std::setw(5) << dof_vals_mapped[l][i] << " ";
    }
    std::cout << std::endl << std::endl;
*/
    for (size_t i=0; i<dof_vals_unmapped[l].size(); ++i)
    {
      BOOST_TEST(std::abs(dof_vals_unmapped[l][i]-dof_vals_mapped[l][i])<1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(fe_transformation) {

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  
  int init_level = 2;
  int num_level = 3;
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,  &rank);

  // Which files should be checked?

  std::vector< std::string > filenames;
  std::vector< TDim > tdims;
  std::vector< GDim > gdims;
  std::vector< RefCellType > cell_types;
  
  filenames.push_back(std::string(datadir) + std::string("two_triangles_2d.inp"));
  tdims.push_back(2);
  gdims.push_back(2);
  cell_types.push_back(REF_CELL_TRI_STD);
   
  filenames.push_back(std::string(datadir) + std::string("two_quads_2d.inp"));
  tdims.push_back(2);
  gdims.push_back(2);
  cell_types.push_back(REF_CELL_QUAD_STD);
  
  filenames.push_back(std::string(datadir) + std::string("two_tetras_3d.inp"));
  tdims.push_back(3);
  gdims.push_back(3);
  cell_types.push_back(REF_CELL_TET_STD);
  
  filenames.push_back(std::string(datadir) + std::string("two_hexas_3d.inp"));
  tdims.push_back(3);
  gdims.push_back(3);
  cell_types.push_back(REF_CELL_HEX_STD);
  
  std::vector< Vec<2,double> > tri_test_points (6);
  tri_test_points[0][0] = 0.;
  tri_test_points[0][1] = 0.;
  tri_test_points[1][0] = 0.5;
  tri_test_points[1][1] = 0.;
  tri_test_points[2][0] = 1.;
  tri_test_points[2][1] = 0.;
  tri_test_points[3][0] = 0.;
  tri_test_points[3][1] = 0.5;
  tri_test_points[4][0] = 0.5;
  tri_test_points[4][1] = 0.5;
  tri_test_points[5][0] = 0.;
  tri_test_points[5][1] = 1.;
        
  std::vector< std::vector< int > > lag_degrees(3);
  for (int l=0; l < lag_degrees.size(); ++l)
  {
    lag_degrees[l].resize(1,l);
  }
  std::vector< std::vector< int > > rt_degrees(3);
  for (int l=0; l < rt_degrees.size(); ++l)
  {
    rt_degrees[l].resize(1,l);
  }
  std::vector< std::vector< int > > bdm_degrees(2);
  for (int l=0; l < bdm_degrees.size(); ++l)
  {
    bdm_degrees[l].resize(1,l+1);
  }
      
      
  for (int test_number = 0; test_number < filenames.size(); ++test_number) 
  {

    std::string filename = filenames.at(test_number);
    TDim tdim = tdims.at(test_number);
    GDim gdim = gdims.at(test_number);

    /////////////////////////////////////
    // mesh

    std::vector<MeshPtr> master_mesh(num_level);
    
    if (rank == 0) 
    {
      master_mesh[0] = read_mesh_from_file(filename.c_str(), gdim, gdim, 0);
      int cur_level = 0;
      while (cur_level < init_level)
      {
        master_mesh[0] = master_mesh[0]->refine();
        cur_level++;
      }
      for (int i=1; i<num_level; ++i)
      {
        master_mesh[i] = master_mesh[i-1]->refine();
      }
    }
    
    std::vector< MeshPtr > mesh(num_level);
    for (int i=0; i<num_level; ++i)
    {
      int num_ref_seq_steps;
      MeshPtr l_mesh = partition_and_distribute(master_mesh[i], 0, MPI_COMM_WORLD, num_ref_seq_steps);
      assert(l_mesh != 0);
  
      SharedVertexTable shared_verts;
      mesh[i] = compute_ghost_cells(*l_mesh, MPI_COMM_WORLD, shared_verts);
    }

    for (int i = 0; i < num_level; ++i) 
    {
      // tests

      CONSOLE_OUTPUT(rank, "Testing " << filename << " on mesh level: " << init_level + i);
#if 1
      // LAG element
      if (gdim == 2 && cell_types[test_number] == REF_CELL_TRI_STD)
      { 
        for (int l=0; l<rt_degrees.size(); ++l)
        {
          int deg = lag_degrees[l][0];
          CONSOLE_OUTPUT(rank, "test LAG space of degree " << deg);

          ConstRefCellPtr<double, 2> ref_cell_tri = ConstRefCellPtr<double, 2>( new RefCellTriStd <double, 2> );
          PTriLag <double, 2> ansatz_tri_lag_2 (ref_cell_tri);
          DofContainerLagrange<double, 2> dof_tri_lag (ref_cell_tri);

          ansatz_tri_lag_2.init(deg, 2);
          dof_tri_lag.init(deg, 2);
          RefElement<double, 2> ref_fe;
          FETransformationStandard<double, 2> fe_trafo (&ref_fe);
          ref_fe.init (&ansatz_tri_lag_2, &dof_tri_lag, &fe_trafo, false, FE_TYPE_LAGRANGE);

          // loop over cells
          for (mesh::EntityIterator it = mesh[i]->begin(2), e_it = mesh[i]->end(2); it != e_it; ++it) 
          {
//          std::cout << "cell index " << it->index() << std::endl;
//          std::cout << "===============" << std::endl;
  
            // Cell
            Entity cell = mesh[i]->get_entity(2, it->index());
            
            // cell trafo
            //const CellTransformation< double, 2 > *c_trafo = rt_space2.fe_manager().get_cell_transformation(it->index());
            CellTrafoPtr< double, 2 > c_trafo = CellTrafoPtr< double, 2 >(new LinearTriangleTransformation<double, 2>(ref_cell_tri));
            std::vector<double> coord_vtx;
            it->get_coordinates(coord_vtx);
            c_trafo->reinit(coord_vtx);
                   
            // DOF container
            DofContainer<double, 2> const * dof = ref_fe.dof_container();
            
            // apply FE transformation (Piola in this case) to map reference nodal basis to physical cell 
            MappingRef2Phys<double, 2, RefElement<double, 2> > * eval_phys 
              = new MappingRef2Phys<double, 2, RefElement<double, 2> > ( &ref_fe, &fe_trafo, c_trafo );
    
            // inverse FE transformation
            MappingPhys2Ref <double, 2, MappingRef2Phys<double, 2, RefElement<double, 2> > > * eval_phys_on_ref
              = new MappingPhys2Ref <double, 2, MappingRef2Phys<double, 2, RefElement<double, 2> > > (eval_phys, &cell, &fe_trafo, c_trafo);
    
            check_test_points<double, 2> (tri_test_points, &ref_fe, eval_phys_on_ref);
            check_dofs<double, 2> (dof, &ref_fe, eval_phys_on_ref); 
          }
          CONSOLE_OUTPUT(rank, "     passed");
        }
      }
#endif
#if 1
      /// BDM element
      if (gdim == 2 && cell_types[test_number] == REF_CELL_TRI_STD)
      { 
        for (int l=0; l<bdm_degrees.size(); ++l)
        {
          int deg = bdm_degrees[l][0];
          CONSOLE_OUTPUT(rank, "test BDM space of degree " << deg);

          ConstRefCellPtr<double, 2> ref_cell_tri = ConstRefCellPtr<double, 2>( new RefCellTriStd <double, 2> );
          PTriLag <double, 2> ansatz_tri_lag_2 (ref_cell_tri);
          ansatz_tri_lag_2.init(deg, 2);
  
          DofContainerRTBDM<double, 2> dof_tri_bdm (ref_cell_tri);
          dof_tri_bdm.init(deg, DOF_CONTAINER_BDM); 
          
          RefElement<double, 2> ref_fe;
          FETransformationContraPiola<double, 2> fe_trafo (&ref_fe);
          ref_fe.init (&ansatz_tri_lag_2, &dof_tri_bdm, &fe_trafo, false, FE_TYPE_BDM);
  
          // loop over cells
          for (mesh::EntityIterator it = mesh[i]->begin(2), e_it = mesh[i]->end(2); it != e_it; ++it) 
          {
//          std::cout << "cell index " << it->index() << std::endl;
//          std::cout << "===============" << std::endl;
  
            // Cell
            Entity cell = mesh[i]->get_entity(2, it->index());
            
            // cell trafo
            //const CellTransformation< double, 2 > *c_trafo = rt_space2.fe_manager().get_cell_transformation(it->index());
            CellTrafoPtr< double, 2 > c_trafo = CellTrafoPtr< double, 2 >(new LinearTriangleTransformation<double, 2>(ref_cell_tri));
            std::vector<double> coord_vtx;
            it->get_coordinates(coord_vtx);
            c_trafo->reinit(coord_vtx);
            
            // DOF container
            DofContainer<double, 2> const * dof = ref_fe.dof_container();
            
            // apply FE transformation (Piola in this case) to map reference nodal basis to physical cell 
            MappingRef2Phys<double, 2, RefElement<double, 2> > * eval_phys 
              = new MappingRef2Phys<double, 2, RefElement<double, 2> > ( &ref_fe, &fe_trafo, c_trafo );
    
            // inverse FE transformation
            MappingPhys2Ref <double, 2, MappingRef2Phys<double, 2, RefElement<double, 2> > > * eval_phys_on_ref
              = new MappingPhys2Ref <double, 2, MappingRef2Phys<double, 2, RefElement<double, 2> > > (eval_phys, &cell, &fe_trafo, c_trafo);
    
            check_test_points<double,2> (tri_test_points, &ref_fe, eval_phys_on_ref);
            check_dofs<double,2> (dof, &ref_fe, eval_phys_on_ref);     
          }
          CONSOLE_OUTPUT(rank, "     passed");
        }
      }
#endif
      
      // RT element
      if (gdim == 2 && cell_types[test_number] == REF_CELL_TRI_STD)
      { 
        for (int l=0; l<rt_degrees.size(); ++l)
        {
          int deg = rt_degrees[l][0];
          CONSOLE_OUTPUT(rank, "test RT space of degree " << deg);

          ConstRefCellPtr<double, 2> ref_cell_tri = ConstRefCellPtr<double, 2>( new RefCellTriStd <double, 2> );
          PTriLag <double, 2> ansatz_tri_lag_2 (ref_cell_tri);
          PTriLag<double, 2> ansatz_tri_rt_1 (ref_cell_tri);
          AugPTriMono<double, 2> ansatz_tri_rt_2 (ref_cell_tri);
          AnsatzSpaceSum<double, 2> ansatz_tri_rt (ref_cell_tri);
          DofContainerRTBDM<double, 2> dof_tri_rt (ref_cell_tri);
        
          ansatz_tri_rt_1.init(deg,2);
          ansatz_tri_rt_2.init(deg);
          ansatz_tri_rt.init(&ansatz_tri_rt_1, &ansatz_tri_rt_2, ANSATZ_RT);
          dof_tri_rt.init(deg, DOF_CONTAINER_RT);  
          RefElement<double, 2> ref_fe;
          FETransformationContraPiola<double, 2> fe_trafo (&ref_fe);
          ref_fe.init (&ansatz_tri_rt, &dof_tri_rt, &fe_trafo, false, FE_TYPE_RT);
  
          // loop over cells
          for (mesh::EntityIterator it = mesh[i]->begin(2), e_it = mesh[i]->end(2); it != e_it; ++it) 
          {
//          std::cout << "cell index " << it->index() << std::endl;
//          std::cout << "===============" << std::endl;
  
            // Cell
            Entity cell = mesh[i]->get_entity(2, it->index());
            
            // cell trafo
            //const CellTransformation< double, 2 > *c_trafo = rt_space2.fe_manager().get_cell_transformation(it->index());
            CellTrafoPtr< double, 2 > c_trafo = CellTrafoPtr< double, 2 > (new LinearTriangleTransformation<double, 2>(ref_cell_tri));
            std::vector<double> coord_vtx;
            it->get_coordinates(coord_vtx);
            c_trafo->reinit(coord_vtx);
            
            // DOF container
            DofContainer<double, 2> const * dof = ref_fe.dof_container();
            
            // apply FE transformation (Piola in this case) to map reference nodal basis to physical cell 
            MappingRef2Phys<double, 2, RefElement<double, 2> > * eval_phys 
              = new MappingRef2Phys<double, 2, RefElement<double, 2> > ( &ref_fe, &fe_trafo, c_trafo );
    
            // inverse FE transformation
            MappingPhys2Ref <double, 2, MappingRef2Phys<double, 2, RefElement<double, 2> > > * eval_phys_on_ref
              = new MappingPhys2Ref <double, 2, MappingRef2Phys<double, 2, RefElement<double, 2> > > (eval_phys, &cell, &fe_trafo, c_trafo);
    
            check_test_points<double,2> (tri_test_points, &ref_fe, eval_phys_on_ref);
            check_dofs<double,2> (dof, &ref_fe, eval_phys_on_ref); 
          }
          CONSOLE_OUTPUT(rank, "     passed");
        }
      }
      
      CONSOLE_OUTPUT(rank, "===============================");
    }

  }

  MPI_Finalize();

  return;
}
