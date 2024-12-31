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

/// \author Martin Baumann, Thomas Gengenbach, Michael Schick
#define BOOST_TEST_MODULE fe_hdiv

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
using namespace hiflow::la;

/// H(div) conformity test
///
/// \brief For a list of files it is analysed, whether the FE ansatz space
/// has continuous normal component
///

const int DEBUG_OUT = 0;
const double pi = 3.14159265;

template <class DataType, int DIM>
struct TestFunction
{
  size_t nb_func() const 
  {
    return 1;
  }
  
  size_t nb_comp() const 
  {
    return DIM;
  }
  
  size_t weight_size() const 
  {
    return DIM;
  }
  
  inline size_t iv2ind (size_t i, size_t var ) const 
  {
    assert (i==0);
    assert (var < DIM);
    return var;
  }
  
  void evaluate(const Entity& cell, const Vec<DIM, DataType> & x, std::vector<DataType>& vals) const
  {
    assert (DIM >= 2);
    assert (space_ != nullptr);
    vals.clear();
    vals.resize(this->nb_comp(), 0.);
    
    //Vec<DIM, DataType> x;
    //space_->get_cell_transformation(cell.index()).transform(ref_pt, x);
    
    vals[0] = std::sin(pi * x[0]) * std::sin(pi * x[0]) 
            * std::sin(2. * pi * x[1]);

    vals[1] = - std::sin(2. * pi * x[0]) 
            * std::sin(pi * x[1]) * std::sin(pi * x[1]);
    
    if (DIM == 3)
    {
      vals[2] = 1.;
    }
  }
  
  VectorSpace<DataType, DIM>* space_;
};

template <class DataType, int DIM>
void check_normal_comp (VectorSpace<DataType,DIM> &space, 
                        MeshPtr mesh_ptr,
                        RefCellType c_type)
{
  // define fe coeff vector
  CoupledVector<DataType> sol;
  sol.Init(MPI_COMM_WORLD, space.la_couplings());
  sol.InitStructure();
  sol.Zeros();
  
      
  TestFunction<DataType, DIM> test_func;
  test_func.space_ = &space;
      
  FeInterNodal<DataType, DIM, TestFunction<DataType, DIM> > fe_inter (space, &test_func, 0);
  fe_inter.interpolate (sol);
      
  //space.interpolate_function (0, test_func, sol);
  
  /*
  int nlocal = sol.size_local();
  for (int i=0; i<nlocal; ++i)
  {
    sol.SetValue(i, static_cast<DataType> (i));
  }
  */
        
  // loop over mesh interfaces
  for (mesh::EntityIterator it = mesh_ptr->begin(DIM-1), e_it = mesh_ptr->end(DIM-1); it != e_it; ++it) 
  {
    // get incident cells
    std::vector<int> cell_ind;
    for (mesh::IncidentEntityIterator c_it = it->begin_incident(2),
         ce_it = it->end_incident(2); c_it != ce_it; ++c_it)
    {
      cell_ind.push_back(c_it->index());
    }
    if (cell_ind.size() != 2)
    {
      continue;
    }
                        
    Entity l_cell = mesh_ptr->get_entity(DIM, cell_ind[0]);
    Entity r_cell = mesh_ptr->get_entity(DIM, cell_ind[1]);
            
    std::vector<DataType> l_coords;
    std::vector<DataType> r_coords;
            
    l_cell.get_coordinates(l_coords);
    r_cell.get_coordinates(r_coords);
          
    CellTransformation<DataType, DIM> * l_trafo;
    CellTransformation<DataType, DIM> * r_trafo;
    ConstRefCellPtr<DataType, 2>  ref_cell_tri = ConstRefCellPtr<DataType, 2> (new RefCellTriStd <DataType, 2>);
    
    if (c_type == REF_CELL_TRI_STD)
    {
      l_trafo = new LinearTriangleTransformation<DataType, 2> (ref_cell_tri);
      r_trafo = new LinearTriangleTransformation<DataType, 2> (ref_cell_tri);
    }
    else
    {
      assert(false);
    }
    
    l_trafo->reinit(l_coords);
    r_trafo->reinit(r_coords);
            
    // get mid point and normal of edge
    std::vector< DataType > coords;
    it->get_coordinates(coords);
            
    Vec<DIM,DataType> mid_point;
    Vec<DIM,DataType> tangent;
    Vec<DIM,DataType> tangent2;
    Vec<DIM,DataType> n;
    
    if (c_type == REF_CELL_TRI_STD || c_type == REF_CELL_QUAD_STD)
    {
      mid_point[0] = 0.5 *(coords[0] + coords[2]);
      mid_point[1] = 0.5 *(coords[1] + coords[3]);
      
      tangent[0] = coords[2] - coords[0];
      tangent[1] = coords[3] - coords[1];
      
      n = normal(tangent);
    }
    else
    {
      assert(false);
    }
    
    if (DEBUG_OUT >= 2)
    {
      std::cout << " ============================== "<< std::endl;
      std::cout << "Mid Point ";
      for (size_t d=0; d<DIM; ++d)
      {
        std::cout << mid_point[d] << " ";
      }
      std::cout << std::endl;
      std::cout << "Tangent ";
      for (size_t d=0; d<DIM; ++d)
      {
        std::cout << tangent[d] << " ";
      }
      std::cout << std::endl;
      std::cout << "Normal ";
      for (size_t d=0; d<DIM; ++d)
      {
        std::cout << n[d] << " ";
      }
      std::cout << std::endl;
    }
    
    DataType eps = 1e-8;
    Vec<DIM,DataType> a_point = mid_point + eps*n;
    Vec<DIM,DataType> b_point = mid_point - eps*n;
    Vec<DIM,DataType> l_point_ref;
    Vec<DIM,DataType> r_point_ref;
    Vec<DIM,DataType> tmp;
    Vec<DIM,DataType> l_point;
    Vec<DIM,DataType> r_point;
            
    if (l_trafo->contains_physical_point(a_point, tmp))
    {
      l_point = a_point;
      r_point = b_point;
    }
    else
    {
      l_point = b_point;
      r_point = a_point;
    }
    
    l_trafo->inverse(l_point, l_point_ref);
    r_trafo->inverse(r_point, r_point_ref);
              
    if (DEBUG_OUT >= 2)
    {
      std::cout << "Left Point ";
      for (size_t d=0; d<DIM; ++d)
      {
        std::cout << l_point[d] << " ";
      }
      std::cout << std::endl;
      std::cout << "Right Point ";
      for (size_t d=0; d<DIM; ++d)
      {
        std::cout << r_point[d] << " ";
      }
      std::cout << std::endl;
    }

    assert (l_trafo->contains_physical_point(l_point, tmp));
    assert (r_trafo->contains_physical_point(r_point, tmp));
            
    // evaluate FE function on both sides of interface
    FeEvalCell<DataType, DIM> fe_eval(space, sol, 0); 
    std::vector<DataType> l_val; 
    std::vector<DataType> r_val; 
    fe_eval.evaluate(r_cell, r_point, r_val);
    fe_eval.evaluate(l_cell, l_point, l_val);

       
    std::vector<DataType> func_r_val;
    std::vector<DataType> func_l_val;

    test_func.evaluate (r_cell, r_point, func_r_val);
    test_func.evaluate (l_cell, l_point, func_l_val);

    assert (l_val.size() == DIM);
    assert (r_val.size() == DIM);
    assert (func_l_val.size() == DIM);
    assert (func_r_val.size() == DIM);
                
    DataType n_times_lval = 0.;
    DataType n_times_rval = 0.;
    DataType n_times_func_lval = 0.;
    DataType n_times_func_rval = 0.;
          
    for (size_t d=0; d<DIM; ++d)
    {
      n_times_lval += n[d] * l_val[d];
      n_times_rval += n[d] * r_val[d];
      
      n_times_func_lval += n[d] * func_l_val[d];
      n_times_func_rval += n[d] * func_r_val[d];
    }
            
    if (DEBUG_OUT >= 2)
    {
      for (size_t d=0; d<DIM; ++d)
      {
        std::cout << d << ": " << l_val[d] << " <> " << r_val[d] << std::endl;
      }
    }
    
    if (DEBUG_OUT >= 1)
    {
      std::cout << "  left normal comp " << std::fixed << std::setw( 7 ) << std::setprecision( 4 ) 
                << n_times_lval << " ,   right normal comp " << n_times_rval << std::endl
                << "f left normal comp " << std::fixed << std::setw( 7 ) << std::setprecision( 4 ) 
                << n_times_func_lval << " , f right normal comp " << n_times_func_rval 
                << std::endl;
    }
    BOOST_TEST(std::abs(n_times_lval - n_times_rval)< 1e-3);
  }
}
        
static const char *datadir = MESH_DATADIR;

BOOST_AUTO_TEST_CASE(fe_hdiv) {

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  
  int init_level = 3;
  int num_level = 1;
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

/*   
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
 */ 
 
  VectorSpace< double, 2 > rt_space2;
  VectorSpace< double, 3 > rt_space3; 
  VectorSpace< double, 2 > bdm_space2;
  VectorSpace< double, 3 > bdm_space3; 

  std::vector< FEType > rt_fe_ansatz (1, FE_TYPE_RT);
  std::vector< FEType > bdm_fe_ansatz (1, FE_TYPE_BDM);
  std::vector< bool > is_cg(1, true);
      
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
        
      // BDM element
      if (gdim == 2 && cell_types[test_number] == REF_CELL_TRI_STD)
      {
        for (int l=0; l<bdm_degrees.size(); ++l)
        {
          CONSOLE_OUTPUT(rank, "test BDM space of degree " << bdm_degrees[l][0]);
          bdm_space2.Init(*mesh[i], bdm_fe_ansatz, is_cg, bdm_degrees[l]);
          check_normal_comp<double, 2> (bdm_space2, mesh[i], cell_types[test_number]);
        }
      }
      // RT element
      if (gdim == 2 && cell_types[test_number] == REF_CELL_TRI_STD)
      {
        for (int l=0; l<rt_degrees.size(); ++l)
        {
          CONSOLE_OUTPUT(rank, "test RT space of degree " << rt_degrees[l][0]);
          rt_space2.Init(*mesh[i], rt_fe_ansatz, is_cg, rt_degrees[l]);
          check_normal_comp<double, 2> (rt_space2, mesh[i], cell_types[test_number]);
        }
      }
      CONSOLE_OUTPUT(rank, "===============================" );
    }
  }

  MPI_Finalize();

  return;
}
