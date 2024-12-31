/*
template<class DataType, int DIM>
int DofContainerBDM<DataType, DIM>::get_dof_id_on_other_cell ( size_t my_dof_id,  
                                                               mesh::Entity const & my_cell, 
                                                               CellTransformation<DataType, DIM> const * my_trans,
                                                               int my_facet_nr,
                                                               mesh::Entity const & other_cell, 
                                                               CellTransformation<DataType, DIM> const * other_trans,
                                                               int other_facet_nr,
                                                               DofContainer<DataType, DIM> const * other_dofs,
                                                               std::vector< mesh::MasterSlave > const &period ) const
{
  assert (my_dof_id < this->nb_dof_on_cell());
  assert (*this == *other_dofs);

  // same cell -> identical dof_id
  if (my_cell.index() == other_cell.index())
  {
    return my_dof_id;
  }

  // different cells
  const size_t facet_dim = this->facet_test_space_->dim();
  const size_t nb_facet_dofs = facet_dim * this->nb_facets_;

  if (my_dof_id >= nb_facet_dofs)
  {
    // cell dof on different cells -> no match
    return -1;
  }

  const size_t dof_offset = my_dof_id % facet_dim;
  return other_facet_nr * facet_dim + dof_offset; 
}
*/


/*
template<class DataType, int DIM>
DofContainerLagrange<DataType, DIM> * DofContainerLagrange<DataType, DIM>::map_to_other_cell ( int my_cell_index, 
                                                                                               CellTransformation<DataType, DIM> const * my_trans,
                                                                                               int other_cell_index, 
                                                                                               CellTransformation<DataType, DIM> const * other_trans,
                                                                                               std::vector< mesh::MasterSlave > const &period ) const
{
  DofContainerLagrange<DataType, DIM>* mapped_dofs = new DofContainerLagrange<DataType, DIM>(this->ref_cell_);
  if (my_cell_index == other_cell_index)
  {
    // identical cell -> no mapping necessary
    mapped_dofs->init(this->degree_);
  }
  else
  {
    for (size_t i=0; i<this->dim_; ++i)
    {
      DofPointEvaluation<DataType, DIM> const * cur_dof = dynamic_cast< DofPointEvaluation<DataType, DIM> const * > ( this->dofs_[i] );
      DofPointEvaluation<DataType, DIM>* mapped_dof = cur_dof->map_to_other_cell ( my_cell_index, my_trans, other_cell_index, other_trans, period);
      if (mapped_dof != nullptr)
      {
        mapped_dofs->push_back(mapped_dof);
      }
      else
      {
        delete mapped_dofs;
        return nullptr;
      }
    }
  }
  return mapped_dofs;
}

template<class DataType, int DIM>
int DofContainerLagrange<DataType, DIM>::get_dof_id_on_other_cell ( size_t my_dof_id,  
                                                                    mesh::Entity const & my_cell, 
                                                                    CellTransformation<DataType, DIM> const * my_trans,
                                                                    int my_facet_nr,
                                                                    mesh::Entity const & other_cell, 
                                                                    CellTransformation<DataType, DIM> const * other_trans,
                                                                    int other_facet_nr,
                                                                    DofContainer<DataType, DIM> const * other_dofs,
                                                                    std::vector< mesh::MasterSlave > const &period ) const
{
  assert (my_dof_id < this->nb_dof_on_cell());
  assert (*this == *other_dofs);

  if (my_cell.index() == other_cell.index())
  {
    return my_dof_id;
  }

  DofPointEvaluation<DataType, DIM> const * dof_i = dynamic_cast< DofPointEvaluation<DataType, DIM> const * > ( this->dofs_[my_dof_id] );
  assert (dof_i != nullptr);
  assert (my_trans != nullptr);
  assert (other_trans != nullptr);

  // map own dof point to other cell
  Coord other_dof_coord;
  bool found = mesh::map_ref_coord_to_other_cell<DataType, DIM> ( dof_i->get_point(),
                                                                  other_dof_coord, 
                                                                  my_trans,
                                                                  other_trans,
                                                                  period );

  int other_dof_id = -1;
  if (found)
  {
    // loop over dof points in other_dofs
    for (size_t j=0; j<other_dofs->nb_dof_on_cell(); ++j)
    {
      DofPointEvaluation<DataType, DIM> const * dof_j = dynamic_cast< DofPointEvaluation<DataType, DIM> const * > ( other_dofs->get_dof(j) );
      assert (dof_j != nullptr);
      
      // compare dof coordinates
      DataType dist = norm(other_dof_coord - dof_j->get_point());
      if (dist < this->ref_cell_->eps())
      {
        return j;
      }
    }
  }
  return other_dof_id;
}
*/


  /*
  DofPointEvaluation<DataType, DIM> * map_to_other_cell ( int my_cell_index, 
                                                          CellTransformation<DataType, DIM> const * my_trans,
                                                          int other_cell_index, 
                                                          CellTransformation<DataType, DIM> const * other_trans,
                                                          std::vector< mesh::MasterSlave > const &period ) const
  {
    if (my_cell_index == other_cell_index)
    {
      // identical cell -> no mapping necessary
      DofPointEvaluation<DataType, DIM>* ptr = new DofPointEvaluation<DataType, DIM>();
      ptr->init(this->dof_coord_, this->ref_cell_);
      return ptr;
    }
    
    assert (my_trans != nullptr);
    assert (other_trans != nullptr);

    // map own dof point to other cell
    Coord other_dof_coord;
    bool found = mesh::map_ref_coord_to_other_cell<DataType, DIM> ( this->dof_coord_,
                                                                    other_dof_coord,
                                                                    my_trans,
                                                                    other_trans,
                                                                    period );
    if (found)
    {
      // create new Point evaluation object
      DofPointEvaluation<DataType, DIM>* ptr = new DofPointEvaluation<DataType, DIM>();
      ptr->init(other_dof_coord, this->ref_cell_);
      return ptr;
    }
    else
    {
      // point is not shared by my_cell and other_cell -> no mapping possible
      return nullptr;
    }
  }
  */
  
  
  /*
  DofContainerLagrange<DataType, DIM> * map_to_other_cell ( int my_cell_index, 
                                                            CellTransformation<DataType, DIM> const * my_trans,
                                                            int other_cell_index, 
                                                            CellTransformation<DataType, DIM> const * other_trans,
                                                            std::vector< mesh::MasterSlave > const &period ) const;

  int get_dof_id_on_other_cell (size_t my_dof_id,  
                                mesh::Entity const & my_cell, 
                                CellTransformation<DataType, DIM> const * my_trans,
                                int my_facet_nr,
                                mesh::Entity const & other_cell, 
                                CellTransformation<DataType, DIM> const * other_trans,
                                int other_facet_nr,
                                DofContainer<DataType, DIM> const * other_dofs,
                                std::vector< mesh::MasterSlave > const &period ) const;
*/


/*
  DofFunctional * map_to_other_cell ( mesh::Entity const & my_cell, 
                                      mesh::Entity const & other_cell, 
                                      CellTransformation<DataType, DIM> const * my_trans,
                                      CellTransformation<DataType, DIM> const * other_trans,
                                      std::vector< mesh::MasterSlave > const &period ) const
  {
    if (my_cell.index() == other_cell.index())
    {
      // identical cell -> no mapping necessary
      DofFunctional* ptr = new DofFacetMoment<DataType, DIM>();
      static_cast<DofFacetMoment<DataType, DIM> *> (ptr)->init( this->test_space_, 
                                                               this->test_vals_,
                                                               this->trial_vals_,
                                                               this->q_weights_,
                                                               this->test_ind_);
      return ptr;
    }
    else
    {
      // cell moments are interior degrees of freedom and thus not shared by multiple cells
      return nullptr;
    }
  }
*/


/*
  DofFunctional * map_to_other_cell ( mesh::Entity const & my_cell, 
                                      mesh::Entity const & other_cell, 
                                      CellTransformation<DataType, DIM> const * my_trans,
                                      CellTransformation<DataType, DIM> const * other_trans,
                                      std::vector< mesh::MasterSlave > const &period ) const
  {
    if (my_cell.index() == other_cell.index())
    {
      // identical cell -> no mapping necessary
      DofFunctional* ptr = new DofCellMoment<DataType, DIM>();
      static_cast<DofCellMoment<DataType, DIM> *> (ptr)->init( this->test_space_, 
                                                               this->test_vals_,
                                                               this->trial_vals_,
                                                               this->q_weights_,
                                                               this->test_ind_,
                                                               this->ref_cell_);
      return ptr;
    }
    else
    {
      // cell moments are interior degrees of freedom and thus not shared by multiple cells
      return nullptr;
    }
  }
*/



  /*
  int get_dof_id_on_other_cell ( size_t my_dof_id,  
                                 mesh::Entity const & my_cell, 
                                 CellTransformation<DataType, DIM> const * my_trans,
                                 int my_facet_nr,
                                 mesh::Entity const & other_cell, 
                                 CellTransformation<DataType, DIM> const * other_trans,
                                 int other_facet_nr,
                                 DofContainer<DataType, DIM> const * other_dofs,
                                 std::vector< mesh::MasterSlave > const &period ) const;
  */
  
  
    // TODO: remove this function if not needed by DofContainerLagrange for performance reasons
  /*
  /// returns < 0 if dof functional my_dof_id is not critical for conformity
  virtual int get_dof_id_on_other_cell ( size_t my_dof_id,
                                         mesh::Entity const & my_cell, 
                                         CellTransformation<DataType, DIM> const * my_trans,
                                         int my_facet_nr,
                                         mesh::Entity const & other_cell, 
                                         CellTransformation<DataType, DIM> const * other_trans,
                                         int other_facet_nr,
                                         DofContainer<DataType, DIM> const * other_dofs,
                                         std::vector< mesh::MasterSlave > const &period ) const = 0;
  */
  
  
  
  
  
  /// for all DoFs of (cellB, ansatzB) the weights of DoFs of (cellA, ansatzA)
/// are calculated
// TODO: debug this function before using it!!
template < class DataType, int DIM >
void NumberingLagrange< DataType, DIM >::compute_weights_lagrange(
    mesh::Entity const &cellA, doffem::RefElement< DataType, DIM > const &ansatzA, int facet_nr_A,
    mesh::Entity const &cellB, doffem::RefElement< DataType, DIM > const &ansatzB, int facet_nr_B,
    std::vector< mesh::MasterSlave > const &period,
    InterpolationWeights &weights) const 
{
#if 0
  // Absolute tolerance for setting a weight to 0 or 1.
  // TODO: check reasonable tolerances
  const DataType COEFFICIENT_EPS = 1.e3 * std::numeric_limits< DataType >::epsilon();

  size_t fe_indA = this->fe_manager_->get_fe_index(&ansatzA, cellA.index());
  size_t fe_indB = this->fe_manager_->get_fe_index(&ansatzB, cellB.index()); 

  const CellTransformation< DataType, DIM > *transA = this->fe_manager_->get_cell_transformation(cellA.index());
  const CellTransformation< DataType, DIM > *transB = this->fe_manager_->get_cell_transformation(cellB.index());

  const int gdim = this->mesh_->gdim();
  const int row_length = ansatzA.nb_dof_on_cell();

  const int col_length = ansatzB.nb_dof_on_cell();
  weights.clear();
  weights.resize(col_length);

  // map current dof from cell B to cell A
  DofContainerLagrange<DataType, DIM> const * dof_B 
    = dynamic_cast< DofContainerLagrange<DataType, DIM> const * > (ansatzB.dof_container());
  
  assert (dof_B != nullptr);

  DofContainerLagrange<DataType, DIM> * dofs_A = dof_B->map_to_other_cell ( cellB.index(), transB, cellA.index(), transA, period );

  // loop over dofs of element B 
  for (size_t i = 0, e_i = col_length; i != e_i; ++i) 
  { 
    // Calculate weights by evaluating A:s shape functions at the
    // mapped dof. The rows in the matrix correspond to
    // the dofs, and the columns to the shape functions.
    weights[i].resize(row_length, 0.);
    if (dofs_A != nullptr) 
    {
      // weights[i] = dof_B_i (phi_A_1), ... , dof_A_i (phi_A_n)
      dofs_A->evaluate(&ansatzA, i, weights[i]);

      LOG_DEBUG(1, "weights " << string_from_range(weights[i].begin(), weights[i].end()));
    }
  }
  
  if (dofs_A != nullptr)
  {
    delete dofs_A;
  }

  // Filter weights -- make sure we get zeros and ones exactly correct.
std::cout << " +++++++ lagrange ++++++++++ " << std::endl; 
  for (size_t i = 0, e_i = weights.size(); i != e_i; ++i) 
  {
    for (size_t j = 0, e_j = row_length; j != e_j; ++j) 
    {
      DataType &w = weights[i][j];
std::cout << " " << w ;
      if (std::abs(w) < COEFFICIENT_EPS) 
      {
        w = 0.;
      } 
      else if (std::abs(w - 1.) < COEFFICIENT_EPS) 
      {
        w = 1.;
      }
    }
std::cout << std::endl;

  }
#endif
}



#if 0
template < class DataType, int DIM >
void NumberingLagrange< DataType, DIM >::compute_weights_general(
    mesh::Entity const &cellA, RefElement< DataType, DIM > const &ansatzA, int facet_nr_A,
    mesh::Entity const &cellB, RefElement< DataType, DIM > const &ansatzB, int facet_nr_B,
    std::vector< mesh::MasterSlave > const &period,
    InterpolationWeights &weights) const 
{
  size_t fe_indA = this->fe_manager_->get_fe_index(&ansatzA, cellA.index());
  size_t fe_indB = this->fe_manager_->get_fe_index(&ansatzB, cellB.index());

  const CellTransformation< DataType, DIM > *transA = this->fe_manager_->get_cell_transformation(cellA.index());
  const CellTransformation< DataType, DIM > *transB = this->fe_manager_->get_cell_transformation(cellB.index());

  const int row_length = ansatzA.nb_dof_on_cell();

  const int col_length = ansatzB.nb_dof_on_cell();
  weights.clear();
  weights.resize(col_length);

  // loop over dofs of element B 
  for (size_t i = 0, e_i = col_length; i != e_i; ++i) 
  { 
    // get corresponding dof on element A
    int dof_id_A = ansatzB.dof_container()->get_dof_id_on_other_cell ( i, 
                                                                           cellB, transB, facet_nr_B, 
                                                                           cellA, transA, facet_nr_A, ansatzA.dof_container(),
                                                                           period );

    weights[i].resize(row_length, 0.);
    if (dof_id_A >= 0)
    {
      assert (dof_id_A < weights[i].size());
      weights[i][dof_id_A] = 1.;
    }
  }
}
#endif 


/// start testing //////////////////////
#if 0
//std::cout << " cell trafo " << std::endl;
std::vector <Vec<DIM,DataType> > coords_A = transA->get_coordinates();
std::vector <Vec<DIM,DataType> > coords_B = transB->get_coordinates();

/*
std::cout << " cell A " << std::endl;
for (size_t l=0; l<coords_A.size(); ++l)
{
  for (size_t d=0; d<DIM; ++d)
  {
    std::cout << coords_A[l][d] << " ";
  }
  std::cout << std::endl;
}

std::cout << " cell B " << std::endl;
for (size_t l=0; l<coords_B.size(); ++l)
{
  for (size_t d=0; d<DIM; ++d)
  {
    std::cout << coords_B[l][d] << " ";
  }
  std::cout << std::endl;
}
*/
/*
  Vec<DIM, DataType> ref_pt_A;
  ref_pt_A[0] = 0.5;
  ref_pt_A[1] = 0.5;
  Vec<DIM, DataType> phys_pt;
  transA->transform(ref_pt_A, phys_pt);
  
  Vec<DIM, DataType> ref_pt_B;
  bool foundB = transB->inverse(phys_pt, ref_pt_B);
  */
  /*
  std::cout << " cell trafo test A " << std::endl;
  Vec<DIM, DataType> phys_pt;
  transA->transform(ref_pt, phys_pt);
  for (size_t d=0; d<DIM; ++d)
  std::cout << phys_pt[d] << " ";
  
  std::cout << std::endl;
  
  std::cout << " cell trafo test B " << std::endl;
  transB->transform(ref_pt, phys_pt);
  for (size_t d=0; d<DIM; ++d)
  std::cout << phys_pt[d] << " ";
  
  std::cout << std::endl;
  */
/*
  std::cout << " -------------------" << std::endl << std::endl;
  std::cout << " -------------------" << std::endl << std::endl;

  std::cout << std::endl << std::endl << std::endl << "FE Trafo Test " << std::endl;
  std::cout << " phys point " << phys_pt[0] << ", " << phys_pt[1] << std::endl;
  std::cout << " ref point A " << ref_pt_A[0] << ", " << ref_pt_A[1] << std::endl;
  std::cout << " ref point B " << ref_pt_B[0] << ", " << ref_pt_B[1] << " : found? " << foundB << std::endl;
  std::cout << " cell index A " <<  cellA.index() << std::endl;
  std::cout << " cell index B " <<  cellB.index() << std::endl;
  std::cout << std::endl;
  */
  /*
  std::cout << " shape_vals_A" << std::endl;
  std::vector< DataType> shape_vals_A(ansatzA.nb_comp() * ansatzA.dim(), 0.);
  ansatzA.N(ref_pt_A, shape_vals_A);
  
  for (int l=0; l<shape_vals_A.size(); ++l)
  {
    std::cout << shape_vals_A[l] << " ";
  }
  std::cout << std::endl;

  std::cout << " mapped_vals_A: " << std::endl;
  std::vector<DataType> mapped_vals_A(shape_vals_A.size(),0.);
  ansatzA.fe_trafo()->map_shape_function_values (*transA, ref_pt_A, 0, ansatzA.dim(), ansatzA.nb_comp(),
                                                shape_vals_A, mapped_vals_A); 
  for (size_t l=0; l < mapped_vals_A.size(); ++l)
  {
    std::cout << mapped_vals_A[l] << " ";
  }
  std::cout << std::endl;
  
  std::cout << "invers_mapped_vals_A" << std::endl;
  std::vector<DataType> inverse_mapped_vals_A(shape_vals_A.size(),0.);
  ansatzA.fe_trafo()->inverse_map_shape_function_values (*transA, ref_pt_A, 0, ansatzA.dim(), ansatzA.nb_comp(),
                                                mapped_vals_A, inverse_mapped_vals_A );
                                            
  for (size_t l=0; l < inverse_mapped_vals_A.size(); ++l)
  {
    std::cout << inverse_mapped_vals_A[l] << " ";
  }
  std::cout << std::endl;
  
  std::cout << " evalA " << std::endl;
  std::vector<DataType> evalA_vals;
  std::vector< Vec<DIM, DataType> > tmp_pts(1,phys_pt);
  evalA_vals = evalA->evaluate (cellA, tmp_pts);
  
  for (int l=0; l<evalA_vals.size(); ++l)
  {
    std::cout << evalA_vals[l] << " ";
  }
  std::cout << std::endl;

  std::cout << " evalA_on_B at " << ref_pt_B[0] << ", " << ref_pt_B[1] << std::endl;
  std::vector<DataType>  evalA_on_B_vals;
  evalA_on_B->evaluate (ref_pt_B,  evalA_on_B_vals);
  
  for (int l=0; l< evalA_on_B_vals.size(); ++l)
  {
    std::cout <<  evalA_on_B_vals[l] << " ";
  }
  std::cout << std::endl;
  
  std::cout << " shape_vals_B at " << ref_pt_B[0] << ", " << ref_pt_B[1] << std::endl;
  std::vector< DataType> shape_vals_B(ansatzB.nb_comp() * ansatzB.dim(), 0.);
  ansatzB.N(ref_pt_B, shape_vals_B);
  
  for (int l=0; l<shape_vals_B.size(); ++l)
  {
    std::cout << shape_vals_B[l] << " ";
  }
  std::cout << std::endl << std::endl;
  
  
  ref_pt_B[0] = 0.25;
  ref_pt_B[1] = 0.;
  std::cout << " evalA_on_B at " << ref_pt_B[0] << ", " << ref_pt_B[1] << std::endl;
  evalA_on_B_vals.clear();
  evalA_on_B->evaluate (ref_pt_B,  evalA_on_B_vals);
  
  for (int l=0; l< evalA_on_B_vals.size(); ++l)
  {
    std::cout <<  evalA_on_B_vals[l] << " ";
  }
  std::cout << std::endl;

  std::cout << " shape_vals_B at " << ref_pt_B[0] << ", " << ref_pt_B[1] << std::endl;
  shape_vals_B.clear();
  shape_vals_B.resize(ansatzB.nb_comp() * ansatzB.dim(), 0.);
  ansatzB.N(ref_pt_B, shape_vals_B);
  
  for (int l=0; l<shape_vals_B.size(); ++l)
  {
    std::cout << shape_vals_B[l] << " ";
  }
  std::cout << std::endl << std::endl;
  
  ref_pt_B[0] = 0.;
  ref_pt_B[1] = 0.5;
  std::cout << " evalA_on_B at " << ref_pt_B[0] << ", " << ref_pt_B[1] << std::endl;
  evalA_on_B_vals.clear();
  evalA_on_B->evaluate (ref_pt_B,  evalA_on_B_vals);
  
  for (int l=0; l< evalA_on_B_vals.size(); ++l)
  {
    std::cout <<  evalA_on_B_vals[l] << " ";
  }
  std::cout << std::endl;
  
  std::cout << " shape_vals_B at " << ref_pt_B[0] << ", " << ref_pt_B[1] << std::endl;
  shape_vals_B.clear();
  shape_vals_B.resize(ansatzB.nb_comp() * ansatzB.dim(), 0.);
  ansatzB.N(ref_pt_B, shape_vals_B);
  
  for (int l=0; l<shape_vals_B.size(); ++l)
  {
    std::cout << shape_vals_B[l] << " ";
  }
  std::cout << std::endl;
  */

/*
  std::cout << "dof eval test" << std::endl;
  std::cout << "evaluate AnsatzA on dofs of A " << std::endl << std::endl;  
  std::cout << " dof container A " << ansatzA.dof_container() << std::endl;
  std::cout << " ansatz A " << &ansatzA << std::endl;
  
  std::vector< std::vector<DataType> >dof_A_of_phi_A;
  ansatzA.dof_container()-> evaluate (&ansatzA, all_dofs,  dof_A_of_phi_A);
  
  for (size_t l = 0; l<col_length; ++l)
  {
    for (size_t k=0; k<row_length; ++k)
    {
      std::cout << dof_A_of_phi_A[l][k] << " ";
    }
    std::cout << std::endl;
  }
  
  std::cout << "evaluate AnsatzB on dofs of B " << std::endl << std::endl;  
  std::cout << " dof container B " << ansatzB.dof_container() << std::endl;
  std::cout << " ansatz B " << &ansatzB << std::endl;
  
  std::vector< std::vector<DataType> >dof_B_of_phi_B;
  ansatzB.dof_container()-> evaluate (&ansatzB, all_dofs,  dof_B_of_phi_B);
  
  for (size_t l = 0; l<col_length; ++l)
  {
    for (size_t k=0; k<row_length; ++k)
    {
      std::cout << dof_B_of_phi_B[l][k] << " ";
    }
    std::cout << std::endl;
  }
*/
  

  std::cout << std::endl << "evaluate AnsatzA on dofs of B " << std::endl << std::endl;  
/// end testing /////////////////////////////
#endif

  void compute_weights_lagrange(mesh::Entity const &cellA, doffem::RefElement< DataType, DIM > const &ansatzA, int facet_nr_A,
                                mesh::Entity const &cellB, doffem::RefElement< DataType, DIM > const &ansatzB, int facet_nr_B,
                                std::vector< mesh::MasterSlave > const &period,
                                InterpolationWeights &weights) const;
                                
                                
                                    // Note: compute_weights_lagrange is buggy -> only use the more general function compute_weights_general!
//  if (is_lagrange_fe)
    if (false)
    {
      compute_weights_lagrange(master_cell, *virtual_ansatz, virtual_facet_nr,
                               slave_cell, slave_ansatz, slave_facet_nr,
                               this->mesh_->get_period(), weights_s2v);
    }
    else
    {
      compute_weights_general(master_cell, *virtual_ansatz, virtual_facet_nr,
                              slave_cell, slave_ansatz, slave_facet_nr,
                              this->mesh_->get_period(), weights_s2v);
    }

