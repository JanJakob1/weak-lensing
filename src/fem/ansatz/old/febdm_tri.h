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

#ifndef __FEM_FEBDM_TRI_H_
#define __FEM_FEBDM_TRI_H_

#include "febdm.h"
#include "quadrature/qgaussline.h"
#include <cassert>
#include <cmath>

namespace hiflow {
namespace doffem {

///
/// \class FEBDMTri febdm_tri.h
/// \brief Brezzi-Douglas-Marini Finite Element on a Triangle
/// \author Philipp Gerstner
///

template < class DataType, int DIM >
class FEBDMTri : public FEBDM< DataType, DIM > {
public:
  typedef Vec<DIM, DataType> Coord;

  /// Default Constructor
  FEBDMTri();

  /// Default Destructor
  ~FEBDMTri();

  std::string get_name() const { return "BDMTriangle"; }

  void set_my_deg(int degree) 
  {
    assert (my_deg_ == -1);
    assert (degree == 1 || degree == 2);
    my_deg_ = degree;
  }


  void N(const Coord &pt, std::vector< DataType > &weight) const;

  void N_x(const Coord &pt, std::vector< DataType > &weight) const;
  void N_y(const Coord &pt, std::vector< DataType > &weight) const;
  void N_z(const Coord &pt, std::vector< DataType > &weight) const;

  void N_xx(const Coord &pt, std::vector< DataType > &weight) const;
  void N_xy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_xz(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yx(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_yz(const Coord &pt, std::vector< DataType > &weight) const;
  void N_zx(const Coord &pt, std::vector< DataType > &weight) const;
  void N_zy(const Coord &pt, std::vector< DataType > &weight) const;
  void N_zz(const Coord &pt, std::vector< DataType > &weight) const;


protected:
  void init_coord();

  //        |\
  //        | \
  //        |  \ edge 1 with normal n_1
  // edge 2 |   \
  // n_2    |    \
  //        |_____\
  //
  //          edge 3 with normal n_3

  // Edge functions satsifying 
  //
  // e_i * n_j = 0 for i != j, i=1,..6, j=1,2,3
  //
  // e_1 * n_1 = 1 for pt = (1-g1, g1) -> nodal point for e1
  // e_1 * n_1 = 0 for pt = (1-g2, g2)
  //
  // e_2 * n_2 = 1 for pt = (0, g1)   -> nodal point for e2
  // e_2 * n_2 = 0 for pt = (0, g2)
  //
  // e_3 * n_3 = 1 for pt = (g1, 0)   -> nodal point for e3
  // e_3 * n_3 = 0 for pt = (g2, 0)
  
  
 
  void e1 (const DataType &g1, const DataType &g2, const Coord &pt, DataType &val_x, DataType &val_y) const 
  {
    const DataType fac = std::sqrt(2) / (g2 - g1);
    val_x = fac * g2 * pt[0];
    val_y = fac * (g2 - 1.) * pt[1];
  }

  void e2 (const DataType &g1, const DataType &g2, const Coord &pt, DataType &val_x, DataType &val_y) const
  {
    const DataType fac = 1. / (g2 - g1);
    val_x = fac * ( g2 * pt[0] + pt[1] - g2);
    val_y = fac * ( (g2 - 1.) * pt[1] ); 
  }
  
  void e3 (const DataType &g1, const DataType &g2, const Coord &pt, DataType &val_x, DataType &val_y) const
  {
    const DataType fac = 1. / (g2 - g1);
    val_x = fac * (g2 - 1.) * pt[0];
    val_y = fac * ( pt[0] + g2 * pt[1] - g2);
  }
  
  void e4 (const DataType &g1, const DataType &g2, const Coord &pt, DataType &val_x, DataType &val_y) const
  {
    const DataType fac = (1. - pt[0] - pt[1]) * std::sqrt(2.) / (g2 - g1);
    val_x = fac * g2 * pt[0];
    val_y = fac * (g2 - 1.) * pt[1];
  }
  
  void e5 (const DataType &g1, const DataType &g2, const Coord &pt, DataType &val_x, DataType &val_y) const
  {
    const DataType fac = pt[0] / (g2 - g1);
    val_x = fac * ( g2 * pt[0] + pt[1] - g2);
    val_y = fac * ( (g2 - 1.) * pt[1] ); 
  }

  void e6 (const DataType &g1, const DataType &g2, const Coord &pt, DataType &val_x, DataType &val_y) const
  {
    const DataType fac = pt[1] / (g2 - g1);
    val_x = fac * (g2 - 1.) * pt[0];
    val_y = fac * ( pt[0] + g2 * pt[1] - g2);
  }
  
  // bubble functions
  void b1 (const Coord &pt, DataType &val_x, DataType &val_y) const 
  {
    val_x = (1. - pt[0] - pt[1]) * pt[0] * pt[1];
    val_y = 0.;
  } 

  void b2 (const Coord &pt, DataType &val_x, DataType &val_y) const 
  {
    val_x = 0.;
    val_y = (1. - pt[0] - pt[1]) * pt[0] * pt[1];
  } 

  // Linear Lagrange polynomial on [g1, g2] with g(s) = 1 if s = g1   
  DataType l (const DataType &g1, const DataType &g2, const DataType &s) const 
  {
    return (s - g2) / (g1 - g2);
  } 


private:
  /// index ordering form plane (x=i, y=j) to vector
  int ij2ind(int i, int j) const;
  
  /// Gauss points on [0,1]
  std::vector<DataType> g_;


};

template < class DataType, int DIM > 
FEBDMTri< DataType, DIM >::FEBDMTri() 
{
  this->my_type_ = FEType< DataType, DIM >::BDM_TRI;

  // initialize reference cell

  assert(this->ref_cell_ == NULL);
  this->ref_cell_ = &(mesh::CellType::get_instance(mesh::CellType::TRIANGLE));

}

template < class DataType, int DIM > 
FEBDMTri< DataType, DIM >::~FEBDMTri() 
{}

template < class DataType, int DIM > 
void FEBDMTri< DataType, DIM >::init_coord() 
{
  assert (DIM == 2);

  // so far, only elements of degree 1 and 2 are implemented
  assert (this->my_deg_ >= 1);
  assert (this->my_deg_ <= 2); 

  // set topological degree
  this->tdim_ = 2;

  // init Gauss points on [0,1]
  if (this->my_deg_ == 1)
  {
    this->g_.resize(2);
    this->g_[0] = 0.5 - std::sqrt(3.) / 6.;
    this->g_[0] = 0.5 + std::sqrt(3.) / 6.;
  }
  else if (this->my_deg_ == 2)
  {
    this->g_.resize(3);
    this->g_[0] = 0.5 - std::sqrt(15.) / 10.;
    this->g_[1] = 0.5;
    this->g_[2] = 0.5 + std::sqrt(15.) / 10.;
  }
  else 
  {
    NOT_YET_IMPLEMENTED;
  }
  

  // Lexicographical ordering
  if (this->my_deg_ == 1)
  {



  }
  else if (this->my_deg_ == 2)
  {


  }
  else
  {
    NOT_YET_IMPLEMENTED;
  }
  














    DataType offset = (1.0 / this->my_deg_);

    this->coord_.clear();
    int nb_dof_on_cell = (((this->my_deg_ + 1) * (this->my_deg_ + 2)) / 2);
    this->coord_.resize(nb_dof_on_cell);

    // Filling coord vector for full triangle by lexicographical strategy
    // and filtering the line coordinates by given mesh ordering strategy
    // which was computed in first step

    int nb_dof_line = this->my_deg_ + 1;

    for (int j = 0; j < nb_dof_line; ++j) { // y axis
      const DataType j_offset = j * offset;
      for (int i = 0; i < nb_dof_line - j; ++i) // x axis
      {
        Coord coord;

        coord[0] = i * offset;
        coord[1] = j_offset;

        this->coord_[ij2ind(i, j)] = coord;
      }
    }
  }
}

/// \details The counting of the dofs on the reference cell is done by the
///          lexicographical numbering strategy. Therefore, beginning on the
///          x coordinate, then continuing on the y coordinate this is achieved
///          by computing the corresponding offsets to consider the restriction
///          given by the triangle which reads y in [0,1], x < 1 - y

template < class DataType, int DIM >
int FEBDMTri< DataType, DIM >::ij2ind(int i, int j) const 
{
  assert (DIM == 2);
  int offset = 0;
  const int nb_dof_line = this->my_deg_ + 1;

  for (int n = 0; n < j; ++n)
    offset += nb_dof_line - n;

  return (i + offset);
}

/// \details 

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N(const Coord &pt,
                                  std::vector< DataType > &weight) const 
{
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell() * 2);

  const DataType deg = static_cast< DataType >(this->my_deg_);

  switch (deg)
  {
    case 1:
      // edge shape functions -> Normal Basis of dimension 6
      // on edge 1
      this->e1(this->g_[0], this->g_[1], pt, weight[0], weight[1]);   // -> dof point s1 = (1-g0, g0)
      this->e1(this->g_[1], this->g_[0], pt, weight[2], weight[3]);   // -> dof point s2 = (1-g1, g1)
      // on edge 2
      this->e2(this->g_[1], this->g_[0], pt, weight[4], weight[5]);   // -> dof point s3 = (0, g1)
      this->e2(this->g_[0], this->g_[1], pt, weight[6], weight[7]);   // -> dof point s4 = (0, g0)
      // on edge 3
      this->e3(this->g_[0], this->g_[1], pt, weight[8], weight[9]);   // -> dof point s5 = (g0, 0)
      this->e3(this->g_[1], this->g_[0], pt, weight[10], weight[11]); // -> dof point s6 = (g1, 0)

      break;
    case 2:
      // edge shape functions -> Normal Basis of dimension 9
      // on edge 1
      this->e1(this->g_[0], this->g_[1], pt, weight[0], weight[1]);   // -> dof point s1 = (1-g0, g0)
      weight[0] *= this->l(this->g_[0], this->g_[2], pt[1]);
      weight[1] *= this->l(this->g_[0], this->g_[2], pt[1]);

      this->e1(this->g_[1], this->g_[2], pt, weight[2], weight[3]);   // -> dof point s2 = (1-g1, g1)
      weight[2] *= this->l(this->g_[1], this->g_[0], pt[1]);
      weight[3] *= this->l(this->g_[1], this->g_[0], pt[1]);
      
      this->e1(this->g_[2], this->g_[0], pt, weight[4], weight[5]);   // -> dof point s3 = (1-g2, g2)
      weight[4] *= this->l(this->g_[2], this->g_[1], pt[1]);
      weight[5] *= this->l(this->g_[2], this->g_[1], pt[1]);

      // on edge 2
      this->e2(this->g_[2], this->g_[0], pt, weight[6], weight[7]);   // -> dof point s4 = (0, g2)
      weight[6] *= this->l(this->g_[2], this->g_[1], pt[1]);
      weight[7] *= this->l(this->g_[2], this->g_[1], pt[1]);

      this->e2(this->g_[1], this->g_[2], pt, weight[8], weight[9]);   // -> dof point s5 = (0, g1)
      weight[8] *= this->l(this->g_[1], this->g_[0], pt[1]);
      weight[9] *= this->l(this->g_[1], this->g_[0], pt[1]);

      this->e2(this->g_[0], this->g_[1], pt, weight[10], weight[11]); // -> dof point s6 = (0, g0)
      weight[10] *= this->l(this->g_[0], this->g_[2], pt[1]);
      weight[11] *= this->l(this->g_[0], this->g_[2], pt[1]);

      // on edge 3
      this->e3(this->g_[0], this->g_[1], pt, weight[12], weight[13]); // -> dof point s7 = (g0, 0)
      weight[12] *= this->l(this->g_[0], this->g_[2], pt[0]);
      weight[13] *= this->l(this->g_[0], this->g_[2], pt[0]);

      this->e3(this->g_[1], this->g_[2], pt, weight[14], weight[15]); // -> dof point s8 = (g1, 0)
      weight[14] *= this->l(this->g_[1], this->g_[0], pt[0]);
      weight[15] *= this->l(this->g_[1], this->g_[0], pt[0]);

      this->e3(this->g_[2], this->g_[0], pt, weight[16], weight[17]); // -> dof point s9 = (g2, 0)
      weight[16] *= this->l(this->g_[2], this->g_[1], pt[0]);
      weight[17] *= this->l(this->g_[2], this->g_[1], pt[0]);

      // interior shape functions -> Non-Normal Basis of dimension 3
      this->e4(this->g_[0], this->g_[1], pt, weight[18], weight[19]);
      this->e5(this->g_[0], this->g_[1], pt, weight[20], weight[21]);
      this->e6(this->g_[0], this->g_[1], pt, weight[22], weight[23]);

      break;
    default:
      NOT_YET_IMPLEMENTED;
      break;
  }
}

/// \details 

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_x(const Coord &pt,
                                    std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  const DataType deg = static_cast< DataType >(this->my_deg_);

  
}

/// \details 

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_y(const Coord &pt,
                                    std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  const DataType deg = static_cast< DataType >(this->my_deg_);
  
}

/// \details There are no z derivatives in 2D

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_z(const Coord &pt,
                                    std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

/// \details 

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_xx(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  const DataType deg = static_cast< DataType >(this->my_deg_);
  
}

/// \details 

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_xy(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  const DataType deg = static_cast< DataType >(this->my_deg_);
  
}

/// \details There are no z derivatives in 2D

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_xz(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

/// \details 

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_yy(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  assert(weight.size() == this->get_nb_dof_on_cell());

  const DataType deg = static_cast< DataType >(this->my_deg_);
  
}

/// \details There are no z derivatives in 2D

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_yz(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}

/// \details There are no z derivatives in 2D

template < class DataType, int DIM >
void FEBDMTri< DataType, DIM >::N_zz(const Coord &pt,
                                     std::vector< DataType > &weight) const {
  assert (DIM == 2);
  weight.assign(weight.size(), 0.);
}


} // namespace doffem
} // namespace hiflow
#endif
