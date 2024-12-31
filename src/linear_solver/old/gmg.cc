
#include "linear_solver/gmg.h"

#include <cassert>
#include <cstdio>
#include <iomanip>
#include <vector>

#include "common/log.h"
#include "linear_algebra/block_matrix.h"
#include "linear_algebra/la_descriptor.h"
#include "linear_algebra/seq_dense_matrix.h"

namespace hiflow {
namespace la {

/// standard constructor




template < class Res, class Pro, class LAD, int DIM >
LinearSolverState GMG<Res, Pro, LAD, DIM >::SolveImplLevel(int level, const VectorType &b, const VectorType &x, VectorType *out) {
	
  LinearSolverState state = kSolverSuccess;
  VectorType bar_x;
  bar_x.CloneFromWithoutContent(x);
  state = S_.Solve(x, bar_x);
  
  VectorType &tmp;
  tmp.CloneFromWithoutContent(d_[level]);
  A_[level].VectorMult(x, tmp);  
  d_[level].CopyFrom(b[level])
  d_level.axpy(-1, tmp);
  
  R_[level-1].interpolate(d_[level], d_[level - 1]);
    
  if (level > 1) {
    VectorType zeros;
    zeros.CloneFromWithoutContent(d_[level]);
    zeros.zeros();
    SolveImplLevel(level - 1, d_[level-1], zeros, e_[level-1]);
	  
  }
  
  else {
    //Smoother corsest, use_operator(A), in BuildImpl für alle lvl
    A_[level-1].Solve(d_[level-1], e_[level-1]); 		//TODO: how to solve directly
  }
  
  P_[level-1].interpolate(e_[level-1], tmp);
  
  VectorType& x_2 = bar_x + tmp;
  
  state = S_.Solve(x_2, out);
	
  return state;
}


template < class Res, class Pro, class LAD, int DIM >
LinearSolverState GMG<Res, Pro, LAD, DIM >::Init(int level, const VectorType & b) {

  if (level == 0) {
    A_[0].Solve(b, x_[0]);
  }
  else {
    //get b on coarser level:
    R_.interpolate(b, b_[level-1]);
        
    Init(level - 1, b_[level-1]);

    P_[level-1].interpolate(x_[level - 1], x_[level]);
    SolveImplLevel(level, b, x_[level], x_[level]);
	
	  
	  
  }
	
}


//alles in .h

} //end namespace la

} // end namespace hiflow
