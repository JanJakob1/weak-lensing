template < class LAD, IMPL MESH_IMPL, int DIM >
void DynamicMeshProblem< LAD, MESH_IMPL, DIM >::setup_auxiliary_interpolation(
    VectorSpace< DataType, DIM > *old_space, VectorSpace< DataType, DIM > *new_space,
    const VectorType &old_vec) {
  assert(old_space != nullptr);
  assert(new_space != nullptr);
  assert(old_space->meshPtr() != nullptr);
  assert(new_space->meshPtr() != nullptr);

  /*
  if ( this->aux_fe_interpolator_ != nullptr )
  {
      this->aux_fe_interpolator_->clear ( );
      delete this->aux_fe_interpolator_;
  }
  */
  if (this->aux_fe_interpolator_ == nullptr) 
  {
    if (this->fe_inter_type_ == FE_INTER_LAGRANGE) 
    {
      this->aux_fe_interpolator_ = new FEInterpolationLagrange< LAD, DIM >();
    } 
#if 0
    else if (this->fe_inter_type_ == FE_INTER_INNERPROD) {
      this->aux_fe_interpolator_ =
          new FEInterpolationInnerProduct< LAD, DIM >();
    }
#endif 
    else {
      std::cout << " Unkown FE interpolation type " << std::endl;
      exit(-1);
    }
  } else {
    this->aux_fe_interpolator_->clear();
  }

#if 0
  if (this->fe_inter_type_ == FE_INTER_INNERPROD) {
    FEInterpolationInnerProduct< LAD, DIM > *tmp_inter =
        dynamic_cast< FEInterpolationInnerProduct< LAD, DIM > * >(
            this->aux_fe_interpolator_);
    assert(tmp_inter != nullptr);

    InnerProductAssembler< LAD, DIM > *local_asm = this->get_fe_inter_asm();
    assert(local_asm != nullptr);

    tmp_inter->set_inner_product_assembler(local_asm,
                                           this->get_fe_inter_coupling_vars());

    tmp_inter->set_dirichlet_bc(this->get_fe_inter_dirichlet_dofs(),
                                this->get_fe_inter_dirichlet_values());

    tmp_inter->set_solver(this->get_fe_inter_solver());
  }
#endif

  this->aux_fe_interpolator_->init(old_space, new_space, old_vec);
}

template < class LAD, IMPL MESH_IMPL, int DIM >
void DynamicMeshProblem< LAD, MESH_IMPL, DIM >::interpolate_vector(
    const typename LAD::VectorType &in_vec, typename LAD::VectorType &out_vec) {
  assert(this->aux_fe_interpolator_->is_initialized());
  this->aux_fe_interpolator_->interpolate(in_vec, out_vec);
}
