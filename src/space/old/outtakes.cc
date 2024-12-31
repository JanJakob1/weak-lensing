#if 0

//////////// Evaluate Scalar FE Function //////////

template < class DataType, int DIM > 
class FeEvalOnCellScalar 
{
  typedef hiflow::doffem::DofID DofID;
  typedef hiflow::Vec<DIM, DataType> Coord;
  
public:
  FeEvalOnCellScalar(const VectorSpace< DataType, DIM > &space, 
                     const hiflow::la::Vector<DataType> &fun, 
                     int var = 0);

  virtual ~FeEvalOnCellScalar() {}

  void operator()(const mesh::Entity &cell, 
                  const Coord &ref_coord, 
                  DataType &value) const
  {
    std::vector< Coord > ref_coords (1, ref_coord);
    std::vector< DataType > values(1);
    this->operator()(cell, ref_coords, values);
    value = values[0];
  }
                  
  void operator()(const mesh::Entity &cell, 
                  const std::vector< Coord > &ref_coords, 
                  std::vector< DataType > &values) const;

protected:
  const VectorSpace< DataType, DIM > &space_;
  const la::Vector<DataType> &fun_;
  int gdim_;

  // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
  std::vector< DofID > id_;
  std::vector< DataType > val_;
  
  int var_;
};

//////////// Evaluate Gradient of Scalar FE Function //////////

template < class DataType, int DIM >
class FeEvalOnCellScalarGrad 
{
  typedef hiflow::doffem::DofID DofID;
  typedef hiflow::Vec<DIM, DataType> Coord;
  
public:

  FeEvalOnCellScalarGrad(const VectorSpace< DataType, DIM > &space, 
                         const la::Vector<DataType> &fun, 
                         int var);

  virtual ~FeEvalOnCellScalarGrad() {}

  void operator()(const mesh::Entity &cell, 
                  const Coord &ref_coord,
                  Vec<DIM,DataType> & value) const
  {
    std::vector< Coord > ref_coords (1, ref_coord);
    std::vector< Vec<DIM, DataType> > values(1);
    this->operator()(cell, ref_coords, values);
    value = values[0];
  }
  
  void operator()(const mesh::Entity &cell, 
                  const std::vector< Coord > &ref_coords, 
                  std::vector< Vec<DIM,DataType> > &values) const;
  
private:
  const VectorSpace< DataType, DIM > &space_;
  const la::Vector<DataType> &fun_;

  int gdim_;

  // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
  std::vector< DofID > id_;
  std::vector< DataType > val_;
  
  int var_;
};

template < class DataType, int DIM > 
class FeEvalOnCellFlow 
{
  typedef hiflow::doffem::DofID DofID;
  typedef hiflow::Vec<DIM, DataType> Coord;

public:
  FeEvalOnCellFlow(const VectorSpace< DataType, DIM > &space, 
                  const la::Vector<DataType> &fun, 
                  std::vector<size_t> flow_vars,
                  FlowTrafo<DataType, DIM> const * map = nullptr);

  virtual ~FeEvalOnCellFlow() {}

  void operator()(const mesh::Entity &cell, 
                  const Coord &ref_coord,
                  Vec<DIM,DataType> &value) const
  {
    std::vector< Coord > ref_coords (1, ref_coord);
    std::vector< Vec<DIM, DataType> > values (1);
    this->operator()(cell, ref_coords, values);
    value = values[0];
  }

  void operator()(const mesh::Entity &cell, 
                  const std::vector< Coord > &ref_coords,
                  std::vector< Vec<DIM,DataType> > &values) const;
      
private:
  const VectorSpace< DataType, DIM > &space_;
  const la::Vector<DataType> &fun_;
  int gdim_;

  // sorted data sets, i.e. 'val' is sorted corresponding to 'id'
  std::vector< DofID > id_;
  std::vector< DataType > val_;
  
  FlowTrafo<DataType, DIM> const * flow_trafo_;
  std::vector<size_t> flow_vars_;
  std::vector<size_t> fe_inds_;
};

#endif


#if 0
///////////////////////////////////////////////////////////////////
/////////////// FeEvalOnCellScalar ////////////////////////////////
///////////////////////////////////////////////////////////////////

template < class DataType, int DIM >
FeEvalOnCellScalar<DataType, DIM>::FeEvalOnCellScalar(const VectorSpace< DataType, DIM > &space, 
                                                      const la::Vector<DataType> &fun, 
                                                      int var)
: space_(space), 
  fun_(fun), 
  var_(var) 
{
  if (space_.nb_subdom() > 1)
  {
    sort_dofs<DataType>(fun_, id_, val_);
  }
  this->gdim_ = this->space_.mesh().gdim();
}

template < class DataType, int DIM >
void FeEvalOnCellScalar<DataType, DIM>::operator()( const mesh::Entity &cell, 
                                                    const std::vector< Coord > &ref_coords, 
                                                    std::vector< DataType > &values) const
{
  const size_t num_points = ref_coords.size();
  const size_t fe_ind = this->space_.fe_manager().var_2_fe(var_);
  const size_t comp = this->space_.fe_manager().var_2_comp(var_);
    
  Element<DataType,DIM> phys_fe (this->space_, cell.index());
  
  std::vector< DataType > dof_values;
  extract_dof_values (this->space_, cell, phys_fe, this->fun_, this->id_, this->val_, fe_ind, dof_values);
  
  values.clear();
  values.resize(num_points, 0.);
  for (int i = 0; i < num_points; ++i) 
  {
    values[i] = phys_fe.evaluate_var (this->var_, ref_coords[i], dof_values);
  }
}

template class FeEvalOnCellScalar <float, 1>;
template class FeEvalOnCellScalar <float, 2>;
template class FeEvalOnCellScalar <float, 3>;
template class FeEvalOnCellScalar <double, 1>;
template class FeEvalOnCellScalar <double, 2>;
template class FeEvalOnCellScalar <double, 3>;

///////////////////////////////////////////////////////////////////
/////////////// FeEvalOnCellScalarGrad ////////////////////////////
///////////////////////////////////////////////////////////////////

template < class DataType, int DIM >
FeEvalOnCellScalarGrad<DataType, DIM>::FeEvalOnCellScalarGrad(const VectorSpace< DataType, DIM > &space, 
                                                              const la::Vector<DataType> &fun, 
                                                              int var)
: space_(space), fun_(fun), var_(var) 
{
  if (space_.nb_subdom() > 1)
  {
    sort_dofs<DataType>(fun_, id_, val_);
  }
  this->gdim_ = this->space_.mesh().gdim();
}

template < class DataType, int DIM >
void FeEvalOnCellScalarGrad<DataType, DIM>::operator()( const mesh::Entity &cell, 
                                                        const std::vector< Coord > &ref_coords, 
                                                        std::vector< Vec<DIM,DataType> > &values) const
{
  const size_t num_points = ref_coords.size();
  const size_t fe_ind = this->space_.fe_manager().var_2_fe(var_);
  const size_t comp = this->space_.fe_manager().var_2_comp(var_);
  Element<DataType,DIM> phys_fe (this->space_, cell.index());
 
  std::vector< DataType > dof_values;
  extract_dof_values (this->space_, cell, phys_fe, this->fun_, this->id_, this->val_, fe_ind, dof_values);
  
  const size_t num_dofs = dof_values.size();

 
  values.clear();
  values.resize(num_points);

  for (int i = 0; i < num_points; ++i) 
  {
    std::vector< Vec<DIM,DataType> > gradients (num_dofs);
    phys_fe.grad_N_var(ref_coords[i], this->var_, gradients);
    
    for (int j = 0; j < num_dofs; ++j) 
    {
      values[i].Axpy(gradients[j], dof_values[j]);
    }
  }
}
template class FeEvalOnCellScalarGrad <float, 1>;
template class FeEvalOnCellScalarGrad <float, 2>;
template class FeEvalOnCellScalarGrad <float, 3>;
template class FeEvalOnCellScalarGrad <double, 1>;
template class FeEvalOnCellScalarGrad <double, 2>;
template class FeEvalOnCellScalarGrad <double, 3>;

///////////////////////////////////////////////////////////////////
/////////////// FeEvalOnCellFlow //////////////////////////////////
///////////////////////////////////////////////////////////////////

template < class DataType, int DIM >
FeEvalOnCellFlow<DataType, DIM>::FeEvalOnCellFlow(const VectorSpace< DataType, DIM > &space, 
                                                  const la::Vector<DataType> &fun, 
                                                  std::vector<size_t> flow_vars,
                                                  FlowTrafo<DataType, DIM> const * map )
: space_(space), 
  fun_(fun), 
  flow_vars_(flow_vars), 
  flow_trafo_ (map) 
{
  // collect fe indices of all variables
  std::set<size_t> tmp;
  for (size_t l=0; l<flow_vars_.size(); ++l)
  {
    tmp.insert(this->space_.var_2_fe(flow_vars_[l]));
  }
  for (std::set<size_t>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
  {
    this->fe_inds_.push_back(*it);
  }

  if (space_.nb_subdom() > 1)
  {
    sort_dofs<DataType>(fun_, id_, val_);
  }
  
  this->gdim_ = this->space_.mesh().gdim();
}

template < class DataType, int DIM >
void FeEvalOnCellFlow<DataType, DIM>::operator()(const mesh::Entity &cell, 
                                                 const std::vector< Coord > &ref_coords,
                                                 std::vector< Vec<DIM,DataType> > &values) const 
{
  // either vector valued FE should be used (size==1) or tensor product of DIM scalar valued FEs
  assert (fe_inds_.size() == 1 || fe_inds_.size() == DIM);
  
  Element<DataType,DIM> phys_fe (this->space_, cell.index());
    
  // exatract dof values  
  std::vector< std::vector< DataType > > dof_values(fe_inds_.size());

  for (size_t d=0; d<fe_inds_.size(); ++d)
  {
    extract_dof_values (this->space_, cell, phys_fe, this->fun_, this->id_, this->val_, fe_inds_[d], dof_values[d]);
  }
    
  const size_t num_points = ref_coords.size();

  values.clear();
  values.resize(num_points);
  for (int i = 0; i < num_points; ++i) 
  {
    // compute physical point corresponding to reference point
    Vec<DIM, DataType> phys_pt;
    doffem::CellTransformation<DataType, DIM> const * cell_trafo = &(this->space_.get_cell_transformation(cell.index()/*, fe_inds[0]*/));
    cell_trafo->transform(ref_coords[i], phys_pt);
      
    // evaluate flow field at reference point 
    Vec<DIM, DataType> flow;
    if (fe_inds_.size() == DIM)
    {
      // flow field is tensor prodduct of DIM scalar fe function
      for (size_t d =0; d<DIM; ++d)
      {
        flow[d] = phys_fe.evaluate_var (this->flow_vars_[d], ref_coords[i], dof_values[d]);
      }
    }
    else if (fe_inds_.size() == 1)
    {
      // flow field is vector valued fe function
      assert (phys_fe.get_fe(fe_inds_[0])->nb_comp() == DIM);
      std::vector<DataType> tmp = phys_fe.evaluate_fe (this->fe_inds_[0], ref_coords[i], dof_values[0]);
      for (size_t d=0; d<DIM; ++d)
      {
        flow[d] = tmp[d];
      }
    }
    else
    {
      assert(0);
    }
    
    // transform flow field
    if (this->flow_trafo_ != nullptr)
    {
      this->flow_trafo_->operator()(phys_pt, flow, values[i]);
    }
    else
    {
      values[i] = flow;
    }
  }
}

template class FeEvalOnCellFlow <float, 1>;
template class FeEvalOnCellFlow <float, 2>;
template class FeEvalOnCellFlow <float, 3>;
template class FeEvalOnCellFlow <double, 1>;
template class FeEvalOnCellFlow <double, 2>;
template class FeEvalOnCellFlow <double, 3>;
#endif
