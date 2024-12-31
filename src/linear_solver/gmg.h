#ifndef HIFLOW_LINEARSOLVER_GMG_H_
#define HIFLOW_LINEARSOLVER_GMG_H_

#include "common/timer.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/coupled_vector.h"
#include "linear_solver/linear_solver.h"
#include "linear_solver/linear_solver_creator.h"
#include "linear_solver/linear_solver_factory.h"
#include "linear_solver/richardson.h"
#include "linear_solver/preconditioner_bjacobi_standard.h"
#include "linear_solver/preconditioner_bjacobi_ext.h"
#include "linear_solver/preconditioner_vanka.h"
#include "space/vector_space.h"
#include "space/fe_interpolation_map.h"
#include "assembly/assembly_assistant.h"
#include "space/fe_evaluation.h"
#include <cmath>
#include <string>
#include <vector>


namespace hiflow {
namespace la {

/// @brief GMG
///
/// Implementation of the geometric multigrid method 
/// 

template < class Application, class Res, class Pro, class LAD, int DIM >
class GMGBase : public LinearSolver< LAD, LAD > 
{
public:
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;
    
  GMGBase()
  : ap_(nullptr), 
  cycle_type_("V"), 
  pre_it_(5), 
  post_it_(5), 
  op_manually_set_(false), 
  rhs_manually_set_(false),
  iterative_initial_(true),
  nb_lvl_(-1),
  omega_(1.),
  initialized_LA_(false),
  initialized_RP_(false),
  update_bc_(true),
  interpolate_rhs_(true)
  {
    this->name_ = "GMG";
  }
  
  virtual ~GMGBase() 
  {
    this->Clear();
  }
  
  virtual void Clear();
    
  inline void set_application (Application * ap) 
  {
    assert(ap!= nullptr);
    this->ap_= ap;
    this->SetModifiedOperator(true);
  }
  
  inline void set_spaces(std::vector < VectorSpace <DataType, DIM >*  >& spaces) 
  {
    this->nb_lvl_ = spaces.size();
    const size_t nb_fe = spaces[0]->nb_fe();
    
    this->spaces_.resize(this->nb_lvl_);
    for (int i = 0; i< this->nb_lvl_; ++i) 
    {
      assert(spaces[i]!=nullptr);
      assert(nb_fe == spaces[i]->nb_fe());
      
      this->spaces_[i] = spaces[i];
    }
    this->initialized_LA_ = false;
    this->initialized_RP_ = false;
    this->SetModifiedOperator(true);
  }
  
  virtual void set_restriction(std::vector <Res *>& R) 
  {
    assert (R.size() == (this->nb_lvl_ - 1));
    
    this->R_.resize(R.size());
    for (int i = 0; i< R.size(); ++i) 
    {
      //assert(R[i]!=nullptr);
      this->R_[i] = R[i];
    }
    this->initialized_RP_ = false;
    this->SetState(false);
  }
  
  virtual void set_prolongation(std::vector < Pro *>& P) 
  {
    assert (P.size() == (this->nb_lvl_ - 1));
  
    this->P_.resize(P.size());
    for (int i = 0; i < P.size(); ++i) 
    {
      //assert(P[i]!=nullptr);
      this->P_[i] = P[i];
    }
    this->initialized_RP_ = false;
    this->SetState(false);
  }
  
  inline void set_smoothers(std::vector <Preconditioner<LAD>* >& S) 
  {
    assert (S.size() == this->nb_lvl_ - 1);
    
    this->preS_.resize(S.size()+1, nullptr);
    for (int i = 0; i< S.size(); ++i) 
    {
      assert(S[i]!=nullptr);
      this->preS_[i+1] = S[i];
    }
    this->SetState(false);
  }
  
  inline void set_coarse_solver(LinearSolver<LAD>* S) 
  {
    assert(S != nullptr);
    this->cS_ = S;
    this->SetState(false);
  }
  
  inline void set_operators(std::vector<OperatorType*>& A) 
  {    
    assert (A.size() == this->nb_lvl_);
    
    //flag
    op_manually_set_ = true;
    
    this->A_.resize(A.size());
    for (int i = 0; i< A.size(); ++i) 
    {
      assert(A[i]!=nullptr);
      this->A_[i] = A[i];
    }
    
    this->SetModifiedOperator(true);
  }
  
  inline void set_rhs (std::vector<VectorType*>& b) 
  {
    assert (b.size() == this->nb_lvl_);
          
    //flag
    rhs_manually_set_ = true;
    this->b_.resize(b.size());
    for (int i = 0; i< b.size(); ++i) 
    {
      assert(b[i]!=nullptr);
      this->b_[i] = b[i];
    }
  }
  
  virtual void InitParameter(std::string cycle_type,
                             bool nested_iteration,
                             int pre_smooth_it,
                             int post_smooth_it,
                             DataType relax_omega,
                             bool use_approx_defect_correction,
                             bool use_transpose_prolongation,
                             bool update_bc,
                             bool interpolate_rhs)
  {
    this->iterative_initial_ = nested_iteration;
   
    assert (pre_smooth_it >= 0);
    assert (post_smooth_it >= 0);
    
    this->pre_it_= pre_smooth_it;
    this->post_it_ = post_smooth_it;
  
    assert (relax_omega > 0.);
    assert (relax_omega <= 1.);
    this->omega_ = relax_omega;
  
    assert(cycle_type == "V" || cycle_type == "W" || cycle_type == "F");
    this->cycle_type_ = cycle_type;
    
    this->use_approx_defect_correction_ = use_approx_defect_correction;
    this->use_transpose_prolongation_ = use_transpose_prolongation;
    this->update_bc_ = update_bc;
    this->interpolate_rhs_ = interpolate_rhs;
  }
  
  void Restrict  (int level, const VectorType* in, VectorType* out);
  
  void Prolongate(int level, const VectorType* in, VectorType* out);
  
  
  inline DataType get_operator_duration() 
  {
    return this->timer_operator_.get_duration();
  }
  
  inline DataType get_rhs_duration() 
  {
    return this->timer_rhs_.get_duration();
  }
  
  inline DataType get_respro_duration() 
  {
    return this->timer_respro_.get_duration();
  }
  
  inline void reinit_LA()
  {
    this->initialized_LA_ = false;
  }

  inline void reinit_RP()
  {
    this->initialized_RP_ = false;
  }
  
  virtual void setup_gmg_vectors ( std::vector<VectorType*>& vec);
    
  /*
  virtual void Init( Application * ap,
                     const std::vector< VectorSpacePtr<DataType, DIM> > all_spaces,
                     const PropertyTree &gmg_param,
                     const PropertyTree &locsolver_param,
                     LinearSolver<LAD>* & coarse_solver,
                     std::vector< Preconditioner<LAD>*> & smoothers,
                     std::vector< std::vector< VectorType* > >& aux_vectors
                   ); */
protected:
  
  LinearSolverState SolveImpl(const VectorType &b, 
                              VectorType *x); 
                                                                           
  virtual void BuildImpl(VectorType const *b, 
                         VectorType *x);
  
  virtual void UpdateOperators(); 
  
  virtual void UpdateDirichletBC();
  
  virtual void InitLA(VectorType const *b, 
                      VectorType *x);
  
  virtual void InitRestriction();
    
  virtual void ApplyDirichletBC(int level, bool zeros, VectorType *x) const;
  
  virtual void BuildRhsLvl(const VectorType &b);
  

  
  bool op_manually_set_, rhs_manually_set_;
  bool iterative_initial_;
  bool initialized_LA_;
  bool initialized_RP_;
  bool use_approx_defect_correction_;
  bool use_transpose_prolongation_;
  bool update_bc_;
  bool interpolate_rhs_;
  
  Application * ap_;
  std::vector < VectorSpace <DataType , DIM >*  > spaces_;
  
  // smoothers and coarse solver
  std::vector < Richardson<LAD>* > S_;
  std::vector < Preconditioner<LAD> * > preS_;
  LinearSolver<LAD>* cS_;

  //Restriction Operator between two successive levels
  std::vector < Res* > R_;
  //Prolongation Operator between two successive levels
  std::vector < Pro* > P_;

  //coefficient matrices of linear system on each level
  std::vector < OperatorType *> A_;
  //right-hand side on each level
  std::vector < VectorType *> b_;
  
  //defect on each level
  std::vector < VectorType *> d_;
  //error on each level
  std::vector < VectorType *> e_;
  //solution on each level
  std::vector < VectorType *> x_;
  
  // residual on each level
  std::vector < VectorType *> r_;

  
  std::vector < VectorType*> tmp_;
  std::vector < VectorType*> zeros_;
  
  // dirichlet BC
  std::vector< std::vector< DataType> > dirichlet_vals_;
  std::vector< std::vector< DataType> > dirichlet_zeros_;
  std::vector< std::vector< int > > dirichlet_dofs_;
  
  //iterations of the smoother at the beginning and at the end of the scheme
  int pre_it_;
  int post_it_;
  std::string cycle_type_;
  DataType omega_;
  
  Timer timer_operator_;
  Timer timer_rhs_;
  Timer timer_respro_;
 
  int nb_lvl_;

private: 
  //Smoothers on each level
  LinearSolverState SolveImplLevel(int level, 
                                   const VectorType &b, 
                                   const VectorType &x, 
                                   VectorType *out);
 
  LinearSolverState NestedIteration(int level, const VectorType & b);
};

template <class Application, class LAD, int DIM >     
class GMGStandard : public GMGBase<Application, FeInterMapFullNodal< LAD, DIM >,
                                    FeInterMapFullNodal< LAD, DIM >, LAD, DIM > {
public:
  
  typedef typename LAD::MatrixType OperatorType;
  typedef typename LAD::VectorType VectorType;
  typedef typename LAD::DataType DataType;
  typedef GMGBase< Application, FeInterMapFullNodal< LAD, DIM >,
                   FeInterMapFullNodal< LAD, DIM >, LAD, DIM > GMGBASE;

  GMGStandard() 
  : GMGBASE()
  {
  }
  
  ~GMGStandard()
  {
    this->Clear();
  }
 
  void set_restriction(std::vector <FeInterMapFullNodal< LAD, DIM > *>& R) 
  {
    LOG_ERROR("Invalid call for this GMG implementation");
    quit_program();
  }
  
  void set_prolongation(std::vector < FeInterMapFullNodal< LAD, DIM > *>& P) 
  {
    LOG_ERROR("Invalid call for this GMG implementation");
    quit_program();
  }
  
  void Clear() 
  {
    for (int i = 0; i < this->Res_.size(); ++ i) 
    {
      if (this->Pro_[i] != nullptr) 
      {  
        delete this->Pro_[i];
        this->Pro_[i] = nullptr;
      }
      if (this->Res_[i] != nullptr) 
      {  
        delete this->Res_[i];
        this->Res_[i] = nullptr;
      }
    }  
    GMGBASE::Clear();
  }
  
protected :

  std::vector < FeInterMapFullNodal< LAD, DIM >* >  Res_;
  std::vector < FeInterMapFullNodal< LAD, DIM >* >  Pro_;


  void InitRestriction() 
  { 
    if (this->initialized_RP_)
    {
      return;
    }
    
    assert (this->spaces_.size() > 1);
    
    /// init Restriction and Prolongation Operators
    for (int i=0; i<this->Res_.size(); ++i)
    {
      if (this->Res_[i] != nullptr)
      {
        delete this->Res_[i];
      }
      if (this->Pro_[i] != nullptr)
      {
        delete this->Pro_[i];
      }
    }
    this->Res_.clear();
    this->Pro_.clear();
    
    this->Res_.resize(this->spaces_.size()-1, nullptr);
    this->Pro_.resize(this->spaces_.size()-1, nullptr);   

    for (int i= 0; i < Res_.size(); i++) 
    {
      this->Pro_[i] = new FeInterMapFullNodal<LAD, DIM > ();

      if (!this->use_transpose_prolongation_)
      {
        this->Res_[i] = new FeInterMapFullNodal<LAD, DIM > ();
      }
    }
  
    GMGBASE::set_prolongation(Pro_);
    GMGBASE::set_restriction(Res_);
    GMGBASE::InitRestriction();
  }
};

/*
template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::Init(
              Application * ap,
              const std::vector< VectorSpacePtr<DataType, DIM> > all_spaces,
              const PropertyTree &gmg_param,
              const PropertyTree &locsolver_param,
              LinearSolver<LAD>* & coarse_solver,
              std::vector< Preconditioner<LAD>* >& smoothers,
              std::vector< std::vector< VectorType* > >& aux_vectors) 
{
  PropertyTree coarse_param = gmg_param["CoarseSolver"];
  PropertyTree smoother_param = gmg_param["KrylovSmoother"];
  
  std::string smoother_type = gmg_param["SmootherType"].get<std::string>();
  
  // spaces
  int nb_lvl = all_spaces.size();
  int gmg_nb_lvl = gmg_param["NumLevel"].get<int>();
  std::vector< VectorSpace<DataType, DIM>* > gmg_spaces(gmg_nb_lvl, nullptr);

  for (int l=0; l<gmg_nb_lvl; ++l)
  {
    assert (all_spaces[nb_lvl - gmg_nb_lvl + l] != 0);
    gmg_spaces[l] = all_spaces[nb_lvl - gmg_nb_lvl + l].get();
  }
  
  // coarse solver
  if (coarse_solver != nullptr)
  {
    if (coarse_solver->GetPreconditioner() != nullptr)
    {
      delete coarse_solver->GetPreconditioner();
    }
    delete coarse_solver;
  }
  
  prepare_krylov_solver<LAD, DIM>(coarse_solver, 
                        coarse_param,
                        locsolver_param,
                        gmg_spaces[0],
                        nullptr); 
   
  // smoothers
  for (int l=0; l<smoothers.size(); ++l)
  {
    if (smoothers[l] != nullptr)
    {
      LinearSolver<LAD>* cast_sm = dynamic_cast<LinearSolver<LAD>* >(smoothers[l]);
      if (cast_sm != 0)
      {
        if (cast_sm->GetPreconditioner() != nullptr)
        {
          delete cast_sm->GetPreconditioner();
        }
      }
      delete smoothers[l];
    }
  }
  
  smoothers.clear();
  smoothers.resize(gmg_nb_lvl-1, nullptr);
    
  for (int l=0; l<gmg_nb_lvl-1; ++l)
  {
    if (smoother_type == "FSAI" || 
        smoother_type == "HiflowILU" || 
        smoother_type == "SOR" || 
        smoother_type == "SSOR" || 
        smoother_type == "Jacobi")
    {
      PreconditionerBlockJacobiStand< LAD >* sm = new PreconditionerBlockJacobiStand< LAD >();
      sm->Init(smoother_type, locsolver_param);
      smoothers[l] = sm;
    }
    else if (smoother_type == "ILUPP" || 
             smoother_type == "MklILU" || 
             smoother_type == "MklLU" || 
             smoother_type == "UmfpackLU")
    {
      PreconditionerBlockJacobiExt< LAD >* sm = new PreconditionerBlockJacobiExt< LAD >();
      sm->Init(smoother_type, locsolver_param);
      smoothers[l] = sm;
    }
    else if (smoother_type == "Vanka")
    {
      PreconditionerVanka< LAD, DIM >* sm = new PreconditionerVanka< LAD, DIM >();
      sm->InitParameter(
        *gmg_spaces[l+1],
         locsolver_param["Vanka"]["Damping"].get< DataType >(),
         locsolver_param["Vanka"]["Iterations"].get< int >(),
         locsolver_param["Vanka"]["BlockCells"].get< bool >());
      smoothers[l] = sm;
    }
    else if (smoother_type == "Krylov")
    {
      LinearSolver<LAD>* sm;
      prepare_krylov_solver<LAD, DIM>(sm, smoother_param, locsolver_param, gmg_spaces[l+1], nullptr);
      smoothers[l] = sm;
    }
  }
  
  this->set_application(ap);
  this->set_spaces(gmg_spaces);
  this->InitControl(gmg_param["Iterations"].get<int>(), 1e-20, 1e-10, 1e6);
  this->InitParameter(gmg_param["CycleType"].get<std::string>(), 
                           gmg_param["NestedIteration"].get<bool>(), 
                           gmg_param["PreSmoothingSteps"].get<int>(), 
                           gmg_param["PostSmoothingSteps"].get<int>(), 
                           gmg_param["SmootherRelax"].get<DataType>(), 
                           false, 
                           gmg_param["TransposedP"].get<bool>(),
                           false,
                           true);
  this->SetPrintLevel(gmg_param["PrintLevel"].get<int>());
  this->set_coarse_solver(coarse_solver);
  this->set_smoothers(smoothers);
  
  // vectors of previous solutions
  for (int k=0; k<aux_vectors.size(); ++k)
  {
    this->setup_gmg_vectors(aux_vectors[k]);
  }
}
*/
template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::setup_gmg_vectors ( std::vector<VectorType*>& vec)
{
  const int num_vec = vec.size();
  for (int i=0; i<num_vec-1; ++i)
  {
    if (vec[i] != nullptr)
    {
      delete vec[i];
    }
  }
  vec.clear();
  vec.resize(this->nb_lvl_, nullptr);
  vec.resize(this->nb_lvl_, nullptr);
  
  for (int i=0; i<this->nb_lvl_; ++i)
  {
    assert (this->spaces_[i] != nullptr);
    
    vec[i] = new VectorType;
    this->ap_->prepare_rhs(i, this->spaces_[i], vec[i]); 
  }
}


template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::Clear() 
{
  LinearSolver<LAD, LAD>::Clear();
  
  for (int i = 0; i < this->nb_lvl_; ++i) 
  {
    if (this->S_[i] != nullptr) 
    {  
      delete this->S_[i];
      this->S_[i] = nullptr;
    }
    if (this->x_[i] != nullptr) 
    {  
      delete this->x_[i];
      this->x_[i] = nullptr;
    }
    if (this->e_[i] != nullptr) 
    {  
      delete this->e_[i];
      this->e_[i] = nullptr;
    }
    if (this->d_[i] != nullptr) 
    {  
      delete this->d_[i];
      this->d_[i] = nullptr;
    }
    if (this->r_[i] != nullptr) 
    {  
      delete this->r_[i];
      this->r_[i] = nullptr;
    }
    if (this->tmp_[i]!= nullptr) 
    {
      delete this->tmp_[i];
      this->tmp_[i] = nullptr;
    }
    if (this->zeros_[i] != nullptr) 
    {  
      delete this->zeros_[i];
      this->zeros_[i] = nullptr;
    }
    if(!op_manually_set_) 
    {
      if (this->A_[i] != nullptr) 
      {
        if (i < this->nb_lvl_ -1)
        {
          delete this->A_[i];
          this->A_[i] = nullptr;
        }
      }
    }
    if (!rhs_manually_set_) 
    {
      if (this->b_[i] != nullptr) 
      {
        delete this->b_[i];
        this->b_[i] = nullptr;
      }
    }
  }
  this->initialized_LA_ = false;
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::BuildImpl(VectorType const *b, VectorType *x) 
{
  assert (b != nullptr);
  assert (x != nullptr);
  
  this->InitLA(b,x);
    
  this->InitRestriction();
  
  this->UpdateDirichletBC();
  
  this->UpdateOperators();
  
  // pass operator to coarse solver
  assert(this->A_[0]!=nullptr);
  this->cS_->SetupOperator(*(this->A_[0]));
  this->cS_->Build(this->b_[0],this->x_[0]);
  
  // setup smoothers  
  for (int i = 1; i < this->nb_lvl_; ++i) 
  {
    assert (this->A_[i]!=nullptr);
    assert (this->preS_[i] != nullptr);
    assert (this->S_[i] != nullptr);
    assert (this->b_[i] != nullptr);
    assert (this->x_[i] != nullptr);
    
    this->S_[i]->SetupOperator(*(this->A_[i]));
    this->S_[i]->Build(this->b_[i], this->x_[i]);
  }
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::UpdateOperators()
{
  if(this->op_manually_set_)
  {
    return;
  }
  
  //LOG_INFO("GMG", "update OP");
  for (int i = this->nb_lvl_-2; i >= 0; --i) 
  {   
    this->ap_->assemble_operator(i,
                                 this->spaces_[i], 
                                 this->dirichlet_dofs_[i],
                                 this->dirichlet_vals_[i],
                                 this->A_[i]); 
  }
  this->A_[this->nb_lvl_-1] = this->op_;
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::UpdateDirichletBC()
{
  if (this->op_manually_set_)
  {
    return;
  }
   
  //LOG_INFO("GMG", "update OP");
  for (int i = 0; i < this->nb_lvl_; ++i) 
  {
    if (update_bc_ || (this->dirichlet_dofs_[i].size() == 0))
    {
      // call application to compute Dirichlet BC
      this->ap_->prepare_dirichlet_bc_solver(i,
                                             this->spaces_[i], 
                                             this->dirichlet_dofs_[i],
                                             this->dirichlet_vals_[i]);
                                
      this->dirichlet_zeros_[i].clear();
      this->dirichlet_zeros_[i].resize(this->dirichlet_dofs_[i].size(), 0.);
    }
  }
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::InitRestriction()
{
  if (this->initialized_RP_)
  {
    return;
  }

  assert (this->R_.size() == this->spaces_.size() - 1);
  assert (this->P_.size() == this->spaces_.size() - 1);
  assert (this->spaces_.size() > 1);
  std::vector <size_t> fe_inds;
    
  for (int i = 0; i < this->spaces_[0]->nb_fe(); ++i) 
  { 
    fe_inds.push_back(i);	   
  }
    
  // note: we assume that restriction and prolongation operators are already set
  this->timer_respro_.reset();
  this->timer_respro_.start();
  for (int i= 0; i < R_.size(); i++) 
  {
    assert (this->P_[i] != nullptr);
    assert (this->spaces_[i] != nullptr);
    assert (this->spaces_[i+1] != nullptr);
    
    this->P_[i]->init(this->spaces_[i], this->spaces_[i+1], *this->b_[i], *this->b_[i+1], fe_inds, fe_inds);

    if (!this->use_transpose_prolongation_)
    {
      assert (this->R_[i] != nullptr);
      this->R_[i]->init(this->spaces_[i+1], this->spaces_[i], *this->b_[i+1], *this->b_[i], fe_inds, fe_inds );
    }
  }
  this->timer_respro_.stop();

  this->initialized_RP_ = true;
}


template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::InitLA(VectorType const *b, VectorType *x) 
{
  if (this->initialized_LA_)
  {
    return;
  }
  
  //LOG_INFO("GMG", "init LA");
  
  // BC
  this->dirichlet_dofs_.clear();
  this->dirichlet_dofs_.resize(this->nb_lvl_);
  this->dirichlet_vals_.clear();
  this->dirichlet_vals_.resize(this->nb_lvl_);
  this->dirichlet_zeros_.clear();
  this->dirichlet_zeros_.resize(this->nb_lvl_);
    
  /// init matrices and vectors
  this->b_.resize(this->nb_lvl_);
      
  if(!op_manually_set_)
  {
    this->timer_operator_.reset();
    this->timer_operator_.start();
    this->A_.resize(this->nb_lvl_);
      
    assert(this->ap_ != nullptr);
    for (int i = 0; i < this->nb_lvl_-1; ++i) 
    {
      if (this->A_[i] != nullptr)
        delete this->A_[i];
        
      this->A_[i] = new OperatorType();
      
      assert(this->spaces_[i] != 0);
      
      // call application to setup matrix object
      this->ap_->prepare_operator(i, this->spaces_[i], this->A_[i]);    
      
      // this-> b_[i] = new VectorType();
      // this->ap_->assemble_rhs(i,
      //                         spaces_[i], 
      //                         this->dirichlet_dofs_[i],
      //                         this->dirichlet_vals_[i],
      //                         b_[i]);
                                     
      this->timer_operator_.stop();
    }
  }
    
  this->timer_rhs_.reset();
  this->timer_rhs_.start();
  
  this->setup_gmg_vectors(this->b_);
  this->setup_gmg_vectors(this->r_);
  this->setup_gmg_vectors(this->x_);
  this->setup_gmg_vectors(this->d_);
  this->setup_gmg_vectors(this->tmp_);
  this->setup_gmg_vectors(this->e_);
  this->setup_gmg_vectors(this->zeros_);  
  
  this->b_[this->nb_lvl_-1]->CopyFrom(*b);
  this->b_[this->nb_lvl_-1]->Update();
  this->timer_rhs_.stop();
    
  // setup smoothers
  this->S_.resize(this->nb_lvl_, nullptr);
  for (int i=0; i<this->S_.size(); ++i)
  {
    if (this->S_[i] != nullptr)
    {
      delete this->S_[i];
    }
    if (i == 0)
      continue;
      
    assert (this->preS_[i] != nullptr);
   
    this->S_[i] = new Richardson<LAD>();
    std::string method = "Preconditioning";
    this->S_[i]->InitParameter(method, this->omega_, this->use_approx_defect_correction_, true);
    this->S_[i]->SetupPreconditioner(*(this->preS_[i]));
  }
  
  this->initialized_LA_ = true;
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::ApplyDirichletBC(int level, bool zeros, VectorType *x) const
{
  if (!(this->dirichlet_dofs_[level].empty())) 
  {
    if (zeros)
    {
      x->SetValues(vec2ptr(this->dirichlet_dofs_[level]), 
                   this->dirichlet_dofs_[level].size(), 
                   vec2ptr(this->dirichlet_vals_[level]));
    }
    else
    {
      x->SetValues(vec2ptr(this->dirichlet_dofs_[level]), 
                   this->dirichlet_dofs_[level].size(), 
                   vec2ptr(this->dirichlet_zeros_[level]));
    }
  }
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::Restrict(int level, const VectorType* in, VectorType* out)
{
  assert (in != nullptr);
  assert (out != nullptr);
  assert (level < this->R_.size());
  
  if (this->R_[level] != nullptr)
  {
    this->R_[level]->interpolate(*in, *out);
  }
  else
  {
    assert (this->P_[level] != nullptr);
    this->P_[level]->interpolate_transpose(*in, *out);
  }
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::Prolongate(int level, const VectorType* in, VectorType* out)
{
  assert (in != nullptr);
  assert (out != nullptr);
  assert (level < this->P_.size());
  
  
  if (this->P_[level] != nullptr)
  {
    this->P_[level]->interpolate(*in, *out);
  }
  else
  {
    assert (this->R_[level] != nullptr);
    this->R_[level]->interpolate_transpose(*in, *out);
  }
}

template < class Application, class Res, class Pro, class LAD, int DIM >
void GMGBase<Application, Res, Pro, LAD, DIM>::BuildRhsLvl(const VectorType &b) 
{
  int num_levels = this->nb_lvl_;
  this->b_[num_levels-1]->CopyFrom(b);
    
  for (int i = num_levels-2; i >=0; --i) 
  {
    if (this->interpolate_rhs_)
    {
      this->Restrict(i, this->b_[i+1], this->b_[i]);
    }
    else
    {
      this->ap_->assemble_rhs(i,
                              spaces_[i], 
                              this->dirichlet_dofs_[i],
                              this->dirichlet_vals_[i],
                              b_[i]);
    }
  }
}

template < class Application, class Res, class Pro, class LAD, int DIM >
LinearSolverState GMGBase<Application, Res, Pro, LAD, DIM>::SolveImpl(const VectorType &b, VectorType *x) 
{
  IterateControl::State conv = IterateControl::kIterate;
    
  int num_levels = this->nb_lvl_;
  this->iter_ = 0;
           
  this->BuildRhsLvl(b);
    
  if (this->iterative_initial_)
  {
    // fill x_ vector by computing x_0 = A_0^-1 b_0
    // and prolongation to higher levels
    LinearSolverState state = NestedIteration(num_levels-1, *(b_[num_levels-1]));
    assert (state == kSolverSuccess);
  }
  else
  {
    this->x_[num_levels-1]->CopyFrom(*x);
  }
    
  // r = b_L - A_L * x_L
  this->r_[num_levels-1]->CloneFrom(*b_[num_levels-1]);
  this->A_[num_levels-1]->VectorMultAdd(-1, *(x_[num_levels-1]), 1, this->r_[num_levels-1]);
    
  this->res_init_ = b.Norm2();
  this->res_ = this->r_[num_levels-1]->Norm2();
  this->res_rel_ = this->res_ / this->res_init_;    
  conv = this->control().Check(this->iter_, this->res_);
    
  if (this->print_level_ > 1) 
  {
    LOG_INFO(this->name_, "initial res norm   =  " << this->res_);
  }
      
  if (conv != IterateControl::kIterate) 
  {
    return kSolverSuccess;
  }
      
  while (conv == IterateControl::kIterate ) 
  {
    this->iter_ += 1;
 
    LinearSolverState state = SolveImplLevel(num_levels -1 , 
                                             *(b_[num_levels -1]), 
                                             *(x_[num_levels -1]), 
                                               x_[num_levels-1]);
      
    // r_L = b_L - A_L * x_L 
    this->r_[num_levels-1]->CopyFrom(*b_[num_levels-1]);
    this->A_[num_levels-1]->VectorMultAdd(-1, *(x_[num_levels-1]), 1, this->r_[num_levels-1]);    
      
    this->res_ = this->r_[num_levels-1]->Norm2();
    this->res_rel_ = this->res_ / this->res_init_;
    
    if (this->print_level_ > 2) 
    {
      LOG_INFO(this->name_, "residual (iteration "<< this->iter_ << "): " << this->res_);
    }
      
    conv = this->control().Check(this->iter_, this->res_);
    if (conv != IterateControl::kIterate) 
    {
      break;
    }
  }
  
  if (this->print_level_ > 1) 
  {
    LOG_INFO(this->name_, "final iterations   = " << this->iter_);
    LOG_INFO(this->name_, "final abs res norm = " << this->res_);
    LOG_INFO(this->name_, "final rel res norm = " << this->res_rel_)
  } 
  
  x->CloneFrom(*x_[num_levels - 1]);
  
  if (conv == IterateControl::kFailureDivergenceTol ||
      conv == IterateControl::kFailureMaxitsExceeded) 
  {
    return kSolverExceeded;
  }
    
  return kSolverSuccess;	  
}

template < class Application, class Res, class Pro, class LAD, int DIM >
LinearSolverState GMGBase<Application, Res, Pro, LAD, DIM>::SolveImplLevel(int level, 
                                                                           const VectorType &b, 
                                                                           const VectorType &x, 
                                                                           VectorType *out)  
{   
  assert (out != nullptr);
  assert (A_[level]!=nullptr);
  assert (d_[level]!=nullptr);
  //assert (R_[level-1]!=nullptr);
  assert (d_[level-1]!=nullptr);
  assert (e_[level-1]!=nullptr);
  assert (tmp_[level-1]!=nullptr);
  assert (S_[level]!=nullptr);
  assert (cS_!=nullptr);
  //assert (P_[level-1]!=nullptr);
    
  this->tmp_[level]->CloneFrom(x);
  this->ApplyDirichletBC(level, false, this->tmp_[level]);  // set dirichlet BC 
  
  LinearSolverState state = kSolverSuccess;
    
  // pre_it steps of pre-smoothing via preconditioned Richardson iteration
  // tmp = S^preit * x
  if (this->pre_it_ > 0)
  {
    this->S_[level]->InitControl(this->pre_it_);
    state = this->S_[level]->Solve(b, this->tmp_[level]);   
  }

  // r_l = b_l - A_l * tmp 
  this->r_[level]->CopyFrom(b);
  this->A_[level]->VectorMultAdd(-1, *this->tmp_[level], 1., this->r_[level]);
  this->ApplyDirichletBC(level, true, this->r_[level]);  // set dirichlet BC zero
  
  //DataType pre_smooth_res = this->r_[level]->Norm2();
  //LOG_INFO("MGV", "residual after pre-smoothing by testing " << pre_smooth_res);
  //LOG_INFO("MGV", "residual after pre-smoothing by solver  " << dynamic_cast<LinearSolver<LAD>*>(S_[level])->res());
  //LOG_INFO("MGV", "iterations of  pre-smoothing by solver  " << dynamic_cast<LinearSolver<LAD>*>(S_[level])->iter());
    
  this->Restrict(level-1, this->r_[level], this->d_[level-1]);
  this->ApplyDirichletBC(level-1, true, this->d_[level-1]);  // set dirichlet BC 
    
  if (level > 1) 
  {
    if( this->cycle_type_ == "V") 
    {
      state = SolveImplLevel(level - 1, *(this->d_[level-1]), *(this->zeros_[level-1]),this->e_[level-1]);
    }
    else 
    {
      VectorType eel1;
      eel1.CloneFromWithoutContent(*this->e_[level-1]);
    
      state = SolveImplLevel(level - 1, *(this->d_[level-1]), *(this->zeros_[level-1]), &eel1);
        
      state = SolveImplLevel(level - 1, *(this->d_[level-1]), eel1, this->e_[level-1]);
    }
  }
  else 
  {
    //solve directly on coarsest level
    // e = A_0^{-1} dl1
    state = this->cS_->Solve(*(this->d_[level-1]), this->e_[level-1]);
    
    if (this->print_level_ > 1)
    {
      LOG_INFO("MGV(" << level << ")", "coarse solved with res " << this->cS_->res());
      LOG_INFO("MGV(" << level << ")", "coarse solved with iter " << this->cS_->iter());
    }
    this->e_[level-1]->Update();
  }
  this->ApplyDirichletBC(level-1, false, this->e_[level-1]);  // set dirichlet BC 
    
  // out = tmp + P_l-1 * e_l-1
  this->Prolongate(level-1, this->e_[level-1], out);
  
  out->Axpy(*this->tmp_[level], 1.);
  this->ApplyDirichletBC(level, false, out);  // set dirichlet BC 
    
  // out = S^postit * out
  if (this->post_it_ > 0)
  {     
    this->S_[level]->InitControl(this->post_it_);  
    state = this->S_[level]->Solve(b, out);
  }
  this->ApplyDirichletBC(level, false, out);  // set dirichlet BC 
  
  return state;
}

template < class Application, class Res, class Pro, class LAD, int DIM >
LinearSolverState GMGBase<Application, Res, Pro, LAD, DIM>::NestedIteration(int level, const VectorType & b) 
{
  LinearSolverState state = kSolverSuccess;
    
  if (level == 0) 
  {
    //solve directly on coarsest mesh
    state = this->cS_->Solve(b, this->x_[0]);
      
    if (state != kSolverSuccess)
      return state;
        
    this->x_[0]->Update();    // necessary ??
  }
  else 
  {
    NestedIteration(level - 1, *(this->b_[level-1]));
    
    this->Prolongate(level-1, this->x_[level - 1], this->x_[level]);
    
    this->x_[level]->Update();// necessary??
      
    state = SolveImplLevel(level, b, *(this->x_[level]), this->x_[level]); 
    if (state != kSolverSuccess)
      return state;     
  }
  return state;
}


} //end namespace la
} //end namespace hiflow

#endif // HIFLOW_LINEARSOLVER_GMG_H_
