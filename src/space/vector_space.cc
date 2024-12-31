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

#include "space/vector_space.h"
#include "dof/dof_partition.h"
#include "dof/numbering_strategy.h"
#include "dof/numbering_lagrange.h"
#include "fem/fe_reference.h"
#include "fem/fe_mapping.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "linear_algebra/vector.h"
#include "mesh/entity.h"
#include "mesh/iterator.h"
#include "mesh/mesh.h"
#include "mesh/types.h"

namespace hiflow {

template < class DataType, int DIM > 
VectorSpace< DataType, DIM >::VectorSpace() 
{
  dof_ = nullptr;
  fe_manager_ = nullptr;
  mesh_ = nullptr;
  comm_ = MPI_COMM_WORLD;
  MPI_Comm_rank(comm_, &rank_);
}

template < class DataType, int DIM >
VectorSpace< DataType, DIM >::VectorSpace(const MPI_Comm &comm) 
{
  dof_ = nullptr;
  fe_manager_ = nullptr;
  mesh_ = nullptr;
  comm_ = comm;
  MPI_Comm_rank(comm, &rank_);
}

template < class DataType, int DIM >
void VectorSpace< DataType, DIM >::get_dof_indices_on_subentity(int cell_index, 
                                                                int tdim, int sindex,
                                                                std::vector< int > *inddof) const 
{
  // TODO high level function in degree_of_freedom
  // now doing this manually

  inddof->resize(this->nb_dof_on_subentity(cell_index, tdim, sindex));
  std::vector< int > inddof_var;

  // loop on every var
  std::vector< int >::iterator iter = inddof->begin();

  for (size_t fe_ind = 0; fe_ind < this->nb_fe(); ++fe_ind) 
  {
    this->get_dof_indices_on_subentity(fe_ind, cell_index, tdim, sindex, &inddof_var);
    // copy values to inddof
    copy(inddof_var.begin(), inddof_var.end(), iter);
    iter += inddof_var.size();
  }
  assert(iter == inddof->end());
}

template < class DataType, int DIM >
void VectorSpace< DataType, DIM >::get_dof_indices(int cell_index,
                                                   std::vector< int > *inddof) const 
{
  // TODO high level function in degree_of_freedom
  // now doing this manually
  inddof->clear();
  inddof->resize(this->nb_dof_on_cell(cell_index));
  std::vector< int > inddof_var;

  // loop on every var
  std::vector< int >::iterator iter = inddof->begin();

  for (size_t fe_ind = 0; fe_ind < this->nb_fe(); ++fe_ind) 
  {
    this->get_dof_indices(fe_ind, cell_index, &inddof_var);
    // copy values to inddof
    copy(inddof_var.begin(), inddof_var.end(), iter);
    iter += inddof_var.size();
  }

  assert(iter == inddof->end());
}

/// Initialize vector space
template < class DataType, int DIM >
void VectorSpace< DataType, DIM >::Init(Mesh &mesh, 
                                        const std::vector< FEType > & fe_types,
                                        const std::vector<bool> & is_cg, 
                                        const std::vector< int > &degrees,
                                        hiflow::doffem::DOF_ORDERING order) 
{
  assert (fe_types.size() == degrees.size());
  for (int v=0; v<degrees.size(); ++v) 
  {
    LOG_INFO("Poly deg for Ansatz " << v << " ", degrees[v]);
  }
  
  std::vector< std::vector< std::vector< int > > > new_deg(degrees.size());
  
  for (size_t l=0; l<fe_types.size(); ++l)
  {
    new_deg[l].resize(1);
    new_deg[l][0].resize(1, degrees[l]);
  }
  this->Init (mesh, fe_types, is_cg, new_deg, order);
}

template < class DataType, int DIM >
void VectorSpace< DataType, DIM >::Init(Mesh &mesh, 
                                        const std::vector< FEType > & fe_types,
                                        const std::vector<bool> & is_cg, 
                                        const std::vector<std::vector< std::vector< int > > > &params,
                                        hiflow::doffem::DOF_ORDERING order) 
{
  // clear
  this->clear();
  
  interminable_assert(params.size() > 0);
  assert (params.size() == fe_types.size());
  assert (is_cg.size() == fe_types.size());

  this->is_dg_.clear();
  this->is_dg_.resize(is_cg.size());
  
  for (size_t l=0; l<is_cg.size(); ++l)
  {
    this->is_dg_[l] = !(is_cg[l]);
  }
  
  LOG_INFO("Finite Element Ansatz", string_from_range(fe_types.begin(), fe_types.end()));
  LOG_INFO("Finite Element DG", string_from_range(is_dg_.begin(), is_dg_.end()));
  
  // set parameters
  const size_t nb_fe = params.size();
  
  const int tdim = mesh.tdim();
  mesh_ = &mesh;

  // create finite element manager
  fe_manager_ = new FEManager();
  fe_manager_->set_mesh(mesh);

  // initialize FiniteElement objects for each cell in mesh
  fe_manager_->init (fe_types, this->is_dg_, params);

  // create dof partition
  dof_ = new DofPartition; // (sequential and parallel)

  // create numbering method
  doffem::NumberingStrategy<DataType, DIM>* number_strategy = new doffem::NumberingLagrange< DataType, DIM >();

  dof_->set_mpi_comm(comm_);
  dof_->set_mesh(mesh_);
  dof_->set_fe_manager(fe_manager_);
  dof_->set_numbering_strategy(number_strategy);
  
  dof_->number (order);
  
  this->setup_la_couplings();
  
  assert (this->mesh_->num_entities(DIM) == this->fe_manager_->fe_tank_size());
      
  delete number_strategy;
}

template < class DataType, int DIM >
void VectorSpace<DataType, DIM>::extract_dof_values (size_t fe_ind, 
                                                     int cell_index,
                                                     const Vector &sol,
                                                     std::vector<DataType> &dof_values) const
{
  assert (cell_index >= 0);
  assert (cell_index < this->mesh_->num_entities(DIM));
  
  const size_t nb_dof_fe = this->nb_dof_on_cell(fe_ind, cell_index);
  size_t start_dof = 0;
  for (size_t i=0; i<fe_ind; ++i)
  { 
    start_dof += this->fe_manager().get_fe(cell_index, i)->dim();
  }
  
  //dof_values.clear();
  dof_values.resize(nb_dof_fe, 0.);
  
  // get dof values from solution vector
  std::vector< doffem::DofID > global_dof_ids;
  this->get_dof_indices(fe_ind, cell_index, &global_dof_ids);
  
  assert (global_dof_ids.size() == nb_dof_fe);

  std::vector< DataType > dof_factors;
  this->dof().get_dof_factors_on_cell(cell_index, dof_factors);
    
  for (size_t i_loc = 0; i_loc < nb_dof_fe; ++i_loc)
  {
    dof_values[i_loc] = sol.GetValue(global_dof_ids[i_loc]) * dof_factors[start_dof + i_loc];

    /*
    if (this->print_)
    {
      std::cout << "Space:extract:  dof id " << global_dof_ids[i_loc] 
                << " , dof factor " << dof_factors[start_dof + i_loc] << std::endl;
    }
    */
  }
}

template < class DataType, int DIM >
void VectorSpace<DataType, DIM>::insert_dof_values (size_t fe_ind, 
                                                    int cell_index,
                                                    Vector &sol,
                                                    const std::vector<DataType> &dof_values) const
{
  assert (cell_index >= 0);
  assert (cell_index < this->mesh_->num_entities(DIM));
  
  const size_t nb_dof_fe = this->nb_dof_on_cell(fe_ind, cell_index);
  assert (dof_values.size() == nb_dof_fe);

  size_t start_dof = 0;
  for (size_t i=0; i<fe_ind; ++i)
  { 
    start_dof += this->fe_manager().get_fe(cell_index, i)->dim();
  }
   
  // get dof values from solution vector
  std::vector< doffem::DofID > global_dof_ids;
  this->get_dof_indices(fe_ind, cell_index, &global_dof_ids);
  
  assert (global_dof_ids.size() == nb_dof_fe);

  std::vector< DataType > dof_factors;
  this->dof().get_dof_factors_on_cell(cell_index, dof_factors);
    
  for (size_t i_loc = 0; i_loc < nb_dof_fe; ++i_loc)
  {
    if (this->rank_ == this->dof().owner_of_dof(global_dof_ids[i_loc]))
    {
      sol.SetValue(global_dof_ids[i_loc], dof_values[i_loc] / dof_factors[start_dof + i_loc]);
    }
  }
}

                                                          
/// Clears allocated dof and femanager

template < class DataType, int DIM > 
void VectorSpace< DataType, DIM >::clear() 
{
  if (this->dof_ != nullptr)
  {
    delete dof_;
    this->dof_ = nullptr;
  }
  
  if (this->fe_manager_!= nullptr)
  {
    delete fe_manager_;
    fe_manager_ = nullptr;
  }
  
  // don't delete mesh object as someone else owns it
  mesh_ = nullptr;
  
  this->is_dg_.clear();
}

template < class DataType, int DIM >
std::vector< int > VectorSpace< DataType, DIM >::get_dof_func(std::vector< size_t > fe_inds) const 
{
  // if standard configuration is desired, the mapping is computed for all
  // variables
  if (fe_inds.empty()) 
  {
    fe_inds.resize(this->nb_fe());
    for (size_t i = 0; i != fe_inds.size(); ++i) 
    {
      fe_inds[i] = i;
    }
  }

  // consecutive numbering of vars
  std::sort(fe_inds.begin(), fe_inds.end());
  std::map< size_t, size_t > vars2consec;
  for (size_t i = 0; i < fe_inds.size(); ++i) 
  {
    vars2consec[fe_inds[i]] = i;
  }

  std::vector< int > dof_func;
  std::map< int, int > dof_map;

  // dof_func_.resize ( space_.dof ( ).ndofs_on_sd ( rank_ ) );
  // const int my_dof_offset = space_.dof ( ).my_dof_offset ( );
  for (hiflow::mesh::EntityIterator
           it = this->mesh_->begin(this->mesh_->tdim()),
           e_it = this->mesh_->end(this->mesh_->tdim());
       it != e_it; ++it) 
  {
    // get global dof indices on current cell
    for (size_t v_i = 0, v_i_e = fe_inds.size(); v_i != v_i_e; ++v_i) 
    {
      const size_t fe_ind = fe_inds[v_i];
      std::vector< int > var_indices;
      this->dof().get_dofs_on_cell(fe_ind, it->index(), var_indices);
      
      // reduce dof indices to only subdomain indices
      for (size_t i = 0, i_e = var_indices.size(); i != i_e; ++i) 
      {
        if (this->dof().is_dof_on_subdom(var_indices[i])) {
          dof_map[var_indices[i]] = vars2consec[fe_ind];
        }
      }
    }
  }

  // create return vector
  dof_func.resize(dof_map.size());
  size_t i = 0;
  for (std::map< int, int >::const_iterator it = dof_map.begin(),
                                            it_e = dof_map.end();
       it != it_e; ++it) 
  {
    dof_func[i] = it->second;
    ++i;
  }

  dof_map.clear();
  return dof_func;
}

template < class DataType, int DIM >
std::vector< int > VectorSpace< DataType, DIM >::get_dofs(const std::vector< size_t > &fe_inds) const 
{
  // check vars vector for correct input
  assert(fe_inds.size() >= 0);

  SortedArray< int > dof_ids;
  dof_ids.reserve(this->dof().nb_dofs_on_subdom(this->dof().my_subdom()));
  for (hiflow::mesh::EntityIterator
           it = this->mesh_->begin(this->mesh_->tdim()),
           e_it = this->mesh_->end(this->mesh_->tdim());
       it != e_it; ++it) 
  {
    // get global dof indices on current cell
    for (size_t v_i = 0, v_i_e = fe_inds.size(); v_i != v_i_e; ++v_i) 
    {
      const size_t fe = fe_inds[v_i];
      std::vector< int > var_indices;
      this->dof().get_dofs_on_cell(fe, it->index(), var_indices);
      for (size_t i = 0, i_e = var_indices.size(); i != i_e; ++i) 
      {
        dof_ids.find_insert(var_indices[i]);
      }
    }
  }

  return dof_ids.data();
}

template < class DataType, int DIM >
std::vector< int > VectorSpace< DataType, DIM >::get_dofs_for_var(const std::vector< size_t > &vars) const 
{
  SortedArray< size_t > fe_inds;
    
  for (size_t l=0; l<vars.size(); ++l)
  {
    fe_inds.find_insert(this->var_2_fe(vars[l]));
  }

  return this->get_dofs(fe_inds.data());
}

template < class DataType, int DIM >
void VectorSpace< DataType, DIM >::setup_la_couplings( )
{
  int my_rank = -1;
  int num_procs = -1;

  MPI_Comm_size(this->dof().get_mpi_comm(), &num_procs);
  assert(num_procs > 0);

  MPI_Comm_rank(this->dof().get_mpi_comm(), &my_rank);
  assert(my_rank >= 0);

  std::vector< int > global_offsets(num_procs + 1, 0);
  std::vector< int > ghost_dofs;
  std::vector< int > ghost_offsets(num_procs + 1, 0);
                                                   
  // Create global_offsets
  for (int i = 0; i < num_procs; ++i) 
  {
    global_offsets[i + 1] = global_offsets[i] + this->dof().nb_dofs_on_subdom(i);
  }

  // Iterate over whole mesh and determine DoFs not owned by current process.
  // Furthermore, for those DoFs, identify owning process and fill data
  // structure
  std::vector< SortedArray< int > > ghost_dofs_per_proc(num_procs);
  typename VectorSpace< DataType, DIM >::MeshEntityIterator mesh_it = this->mesh().begin(this->tdim());
  typename VectorSpace< DataType, DIM >::MeshEntityIterator e_mesh_it = this->mesh().end(this->tdim());
  for (typename VectorSpace< DataType, DIM >::MeshEntityIterator it = mesh_it;
       it != e_mesh_it; ++it) 
  {
    for (size_t fe_ind = 0; fe_ind < this->nb_fe(); ++fe_ind) 
    {
      std::vector< int > dofs_cell_var;
      this->get_dof_indices(fe_ind, it->index(), &dofs_cell_var);
      for (size_t j = 0; j < dofs_cell_var.size(); ++j) 
      {
        if (!this->dof().is_dof_on_subdom(dofs_cell_var[j])) 
        {
          ghost_dofs_per_proc[this->dof().owner_of_dof(dofs_cell_var[j])].find_insert(dofs_cell_var[j]);
        }
      }
    }
  }

  ghost_dofs.clear();
  for (int i = 0; i < num_procs; ++i) 
  {
    ghost_dofs.insert(ghost_dofs.end(), ghost_dofs_per_proc[i].begin(),
                      ghost_dofs_per_proc[i].end());
  }

  for (int i = 0; i < num_procs; ++i) 
  {
    ghost_offsets[i + 1] = ghost_offsets[i] + ghost_dofs_per_proc[i].size();
  }
  
  this->couplings_.Clear();
  this->couplings_.Init(this->get_mpi_comm());
  this->couplings_.InitializeCouplings(global_offsets, ghost_dofs, ghost_offsets);
}

template < class DataType, int DIM >
void VectorSpace< DataType, DIM >::global_2_local_and_ghost (const std::vector<DofId>& global_ids,
                                                             std::vector<DofId>& local_ids,
                                                             std::vector<DofId>& ghost_ids)  const 
{
  const size_t nb_dof = global_ids.size();
  const int ownership_begin = this->couplings_.dof_offset(this->rank_);
  const int ownership_end = ownership_begin + this->couplings_.nb_dofs(this->rank_);
  
  local_ids.clear();
  ghost_ids.clear();
  
  if (nb_dof == 0)
  {
    return;
  }
  
  local_ids.resize(nb_dof, -1);
  ghost_ids.resize(nb_dof, -1);
  
  for (size_t i=0; i!=nb_dof; ++i)
  {
    const DofId gl_i = global_ids[i];
    if ((gl_i >= ownership_begin) && (gl_i < ownership_end))
    {
      // local dof
      local_ids[i] = gl_i - ownership_begin;
    }
    else if (this->couplings_.global2offdiag().find(gl_i) !=
             this->couplings_.global2offdiag().end())
    {
      // ghost dof
      ghost_ids[i] = this->couplings_.Global2Offdiag(gl_i);
    }
  }
}

template class VectorSpace<float, 3 >;
template class VectorSpace<float, 2 >;
template class VectorSpace<float, 1 >;

template class VectorSpace<double, 3 >;
template class VectorSpace<double, 2 >;
template class VectorSpace<double, 1 >;

template < class DataType, int DIM >
void interpolate_constrained_std_vector(const VectorSpace< DataType, DIM > &space,
                                        la::Vector< DataType > &vector) 
{
  // TODO: necessary to take sub-domain into account here?
  const doffem::DofInterpolation<DataType> &interp = space.dof().dof_interpolation();

  const size_t num_constrained_dofs = interp.size();
  
  //LOG_INFO("# constrained dofs", num_constrained_dofs);
  
  // return early if there are no constrained dofs
  if (num_constrained_dofs == 0) {
    return;
  }

  std::vector< int > constrained_dofs;
  constrained_dofs.reserve(num_constrained_dofs);
  std::vector< DataType > constrained_values;
  constrained_values.reserve(num_constrained_dofs);

  std::vector< int > dependencies;
  std::vector< DataType > coefficients;
  std::vector< DataType > dependency_values;

  for (typename doffem::DofInterpolation<DataType>::const_iterator it = interp.begin(),
                                                end = interp.end(); it != end; ++it) {
    if (space.dof().is_dof_on_subdom(it->first)) {

      const size_t num_dependencies = it->second.size();
      LOG_DEBUG(1, " num_dep " << num_dependencies);
      // probably should not happen, but we check to avoid later problems
      if (num_dependencies > 0) {

        dependencies.resize(num_dependencies, -1);
        coefficients.resize(it->second.size(), 0.);
        dependency_values.resize(num_dependencies, 0.);

        for (size_t i = 0; i != num_dependencies; ++i) {
          dependencies[i] = it->second[i].first;  // id of dependencies
          coefficients[i] = it->second[i].second; // coefficient of dependencies
        }

        LOG_DEBUG(1,
                  string_from_range(dependencies.begin(), dependencies.end()));
        LOG_DEBUG(1,
                  string_from_range(coefficients.begin(), coefficients.end()));

        // get values of dependency dofs from vector
        vector.GetValues(&dependencies.front(), num_dependencies,
                         &dependency_values.front());

        // compute dot product of dependency_values and coefficients
        DataType val = 0.;
        for (size_t i = 0; i != num_dependencies; ++i) {
          val += coefficients[i] * dependency_values[i];
        }

        // store information
        constrained_dofs.push_back(it->first);
        constrained_values.push_back(val);
      } else {
        LOG_INFO("Constrained DoF without dependencies found", true);
        quit_program();
      }
    }
  }

  // write constrained dofs to vector
  vector.SetValues(&constrained_dofs.front(), constrained_dofs.size(),
                   &constrained_values.front());
}

template void interpolate_constrained_std_vector<double, 3> (const VectorSpace< double, 3 > &, la::Vector< double > &);
template void interpolate_constrained_std_vector<double, 2> (const VectorSpace< double, 2 > &, la::Vector< double > &);
template void interpolate_constrained_std_vector<double, 1> (const VectorSpace< double, 1 > &, la::Vector< double > &);
template void interpolate_constrained_std_vector<float, 3> (const VectorSpace< float, 3 > &, la::Vector< float > &);
template void interpolate_constrained_std_vector<float, 2> (const VectorSpace< float, 2 > &, la::Vector< float > &);
template void interpolate_constrained_std_vector<float, 1> (const VectorSpace< float, 1 > &, la::Vector< float > &);

template < class DataType, int DIM >
void interpolate_constrained_vector( const VectorSpace< DataType, DIM > &space, la::Vector< DataType > &vector)
{
  interpolate_constrained_std_vector<DataType, DIM>(space, vector);
}

template void interpolate_constrained_vector<double, 3> (const VectorSpace< double, 3 > &, la::Vector< double > &);
template void interpolate_constrained_vector<double, 2> (const VectorSpace< double, 2 > &, la::Vector< double > &);
template void interpolate_constrained_vector<double, 1> (const VectorSpace< double, 1 > &, la::Vector< double > &);
template void interpolate_constrained_vector<float, 3> (const VectorSpace< float, 3 > &, la::Vector< float > &);
template void interpolate_constrained_vector<float, 2> (const VectorSpace< float, 2 > &, la::Vector< float > &);
template void interpolate_constrained_vector<float, 1> (const VectorSpace< float, 1 > &, la::Vector< float > &);


} // namespace hiflow
