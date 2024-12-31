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

#include "fem/cell_trafo/cell_transformation.h"

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow {
namespace doffem {
    
  template class CellTransformation <double, 3>;
    
/*
template<>
void CellTransformation<double,1>::transform ( const Vec<1,double>& coord_ref, Vec<1,double>& coord_mapped ) const 
{
  coord_mapped[0] = this->x ( coord_ref );
}

template<>
void CellTransformation<float,1>::transform ( const Vec<1,float>& coord_ref, Vec<1,float>& coord_mapped ) const 
{
  coord_mapped[0] = this->x ( coord_ref );
}

template<>
void CellTransformation<double,2>::transform ( const Vec<2,double>& coord_ref, Vec<2,double>& coord_mapped ) const 
{
  coord_mapped[0] = this->x ( coord_ref );
  coord_mapped[1] = this->y ( coord_ref );
}

template<>
void CellTransformation<float,2>::transform ( const Vec<2,float>& coord_ref, Vec<2,float>& coord_mapped ) const 
{
  coord_mapped[0] = this->x ( coord_ref );
  coord_mapped[1] = this->y ( coord_ref );
}

template<>
void CellTransformation<double,3>::transform ( const Vec<3,double>& coord_ref, Vec<3,double>& coord_mapped ) const 
{
  coord_mapped[0] = this->x ( coord_ref );
  coord_mapped[1] = this->y ( coord_ref );
  coord_mapped[2] = this->z ( coord_ref );
}

template<>
void CellTransformation<float,3>::transform ( const Vec<3,float>& coord_ref, Vec<3,float>& coord_mapped ) const 
{
  coord_mapped[0] = this->x ( coord_ref );
  coord_mapped[1] = this->y ( coord_ref );
  coord_mapped[2] = this->z ( coord_ref );
}
/*
template<>
void CellTransformation<double,1>::jacobi ( const Vec<1,double>& coord, Mat<1, 1, double>& mat ) const 
{
  mat(0,0) = this->x_x(coord);
}

template<>
void CellTransformation<float,1>::jacobi ( const Vec<1,float>& coord, Mat<1, 1, float>& mat ) const 
{
  mat(0,0) = this->x_x(coord);
}

template<>
void CellTransformation<double,2>::jacobi ( const Vec<2,double>& coord, Mat<2, 2, double>& mat ) const 
{
  mat(0,0) = this->x_x(coord);
  mat(0,1) = this->x_y(coord);
  mat(1,0) = this->y_x(coord);
  mat(1,1) = this->y_y(coord);
}

template<>
void CellTransformation<float,2>::jacobi ( const Vec<2,float>& coord, Mat<2, 2, float>&  mat ) const 
{
  mat(0,0) = this->x_x(coord);
  mat(0,1) = this->x_y(coord);
  mat(1,0) = this->y_x(coord);
  mat(1,1) = this->y_y(coord);
}

template<>
void CellTransformation<double,3>::jacobi ( const Vec<3,double>& coord, Mat<3, 3, double>& mat ) const 
{
  mat(0,0) = this->x_x(coord);
  mat(0,1) = this->x_y(coord);
  mat(0,2) = this->x_z(coord);
  
  mat(1,0) = this->y_x(coord);
  mat(1,1) = this->y_y(coord);
  mat(1,2) = this->y_z(coord);

  mat(2,0) = this->z_x(coord);
  mat(2,1) = this->z_y(coord);
  mat(2,2) = this->z_z(coord);
}

template<>
void CellTransformation<float,3>::jacobi ( const Vec<3,float>& coord, Mat<3, 3, float>& mat) const 
{
  mat(0,0) = this->x_x(coord);
  mat(0,1) = this->x_y(coord);
  mat(0,2) = this->x_z(coord);
  
  mat(1,0) = this->y_x(coord);
  mat(1,1) = this->y_y(coord);
  mat(1,2) = this->y_z(coord);

  mat(2,0) = this->z_x(coord);
  mat(2,1) = this->z_y(coord);
  mat(2,2) = this->z_z(coord);
}
* */
* */
/*
template<>
void CellTransformation<float,1>::hessian ( const Vec<1,float>& coord, size_t d, Mat<1, 1, float>& mat) const 
{
  switch (d)
  {
    case 0:
      mat(0,0) = this->x_xx(coord);
      break;
    default:
      assert(0);
      break;
  }
}

template<>
void CellTransformation<double,1>::hessian ( const Vec<1,double>& coord, size_t d, Mat<1, 1, double>& mat) const 
{
  switch (d)
  {
    case 0:
      mat(0,0) = this->x_xx(coord);
      break;
    default:
      assert(0);
      break;
  }
}

template<>
void CellTransformation<float,2>::hessian ( const Vec<2,float>& coord, size_t d, Mat<2, 2, float>& mat) const 
{
  switch (d)
  {
    case 0:
      mat(0,0) = this->x_xx(coord);
      mat(0,1) = this->x_xy(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->x_yy(coord);
      break;
    case 1:
      mat(0,0) = this->y_xx(coord);
      mat(0,1) = this->y_xy(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->y_yy(coord);
      break;
    default:
      assert(0);
      break;
  }
}

template<>
void CellTransformation<double,2>::hessian ( const Vec<2,double>& coord, size_t d, Mat<2, 2, double>& mat) const 
{
  switch (d)
  {
    case 0:
      mat(0,0) = this->x_xx(coord);
      mat(0,1) = this->x_xy(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->x_yy(coord);
      break;
    case 1:
      mat(0,0) = this->y_xx(coord);
      mat(0,1) = this->y_xy(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->y_yy(coord);
      break;
    default:
      assert(0);
      break;
  }
}

template<>
void CellTransformation<float,3>::hessian ( const Vec<3,float>& coord, size_t d, Mat<3, 3, float>& mat) const 
{
  switch (DIM)
  {
    case 1:
      switch (d)
      {
        case 0:
          mat(0,0) = this->x_xx(coord);
          break;
        default:
          assert(0);
          break;
      }
      break;
    case 2:
      switch (d)
      {
        case 0:
          mat(0,0) = this->x_xx(coord);
          mat(0,1) = this->x_xy(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->x_yy(coord);
          break;
        case 1:
          mat(0,0) = this->y_xx(coord);
          mat(0,1) = this->y_xy(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->y_yy(coord);
          break;
        default:
          assert(0);
          break;
      }
      break;
    case 3:
      switch (d)
      {
        case 0:
          mat(0,0) = this->x_xx(coord);
          mat(0,1) = this->x_xy(coord);
          mat(0,2) = this->x_xz(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->x_yy(coord);
          mat(1,2) = this->x_yz(coord);
          mat(2,0) = mat(0,2);
          mat(2,1) = mat(1,2);
          mat(2,2) = this->x_zz(coord);
          break;
        case 1:
          mat(0,0) = this->y_xx(coord);
          mat(0,1) = this->y_xy(coord);
          mat(0,2) = this->y_xz(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->y_yy(coord);
          mat(1,2) = this->y_yz(coord);
          mat(2,0) = mat(0,2);
          mat(2,1) = mat(1,2);
          mat(2,2) = this->y_zz(coord);
          break;
        case 2:
          mat(0,0) = this->z_xx(coord);
          mat(0,1) = this->z_xy(coord);
          mat(0,2) = this->z_xz(coord);
          mat(1,0) = mat(0,1);
          mat(1,1) = this->z_yy(coord);
          mat(1,2) = this->z_yz(coord);
          mat(2,0) = mat(0,2);
          mat(2,1) = mat(1,2);
          mat(2,2) = this->z_zz(coord);
          break;
        default:
          assert(0);
          break;
      }
      break;
    default:
      assert (0);
      break;
  }
}

template<>
void CellTransformation<double,3>::hessian ( const Vec<3,double>& coord, size_t d, Mat<3, 3, double>& mat) const 
{
  switch (d)
  {
    case 0:
      mat(0,0) = this->x_xx(coord);
      mat(0,1) = this->x_xy(coord);
      mat(0,2) = this->x_xz(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->x_yy(coord);
      mat(1,2) = this->x_yz(coord);
      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);
      mat(2,2) = this->x_zz(coord);
      break;
    case 1:
      mat(0,0) = this->y_xx(coord);
      mat(0,1) = this->y_xy(coord);
      mat(0,2) = this->y_xz(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->y_yy(coord);
      mat(1,2) = this->y_yz(coord);
      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);
      mat(2,2) = this->y_zz(coord);
      break;
    case 2:
      mat(0,0) = this->z_xx(coord);
      mat(0,1) = this->z_xy(coord);
      mat(0,2) = this->z_xz(coord);
      mat(1,0) = mat(0,1);
      mat(1,1) = this->z_yy(coord);
      mat(1,2) = this->z_yz(coord);
      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);
      mat(2,2) = this->z_zz(coord);
      break;
    default:
      assert(0);
      break;
  }
}
*/
/*
template<>
double CellTransformation<double,3>::detJ_x ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  return  this->x_xx(pt) * (this->y_y(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_y(pt))
        + this->x_x(pt)  * (this->y_xy(pt) * this->z_z(pt) + this->y_y(pt) * this->z_xz(pt) 
                          - this->y_xz(pt) * this->z_y(pt) - this->y_z(pt) * this->z_xy(pt))
        - this->x_xy(pt) * (this->y_x(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
        - this->x_y(pt)  * (this->y_xx(pt) * this->z_z(pt) + this->y_x(pt) * this->z_xz(pt) 
                          - this->y_xz(pt) * this->z_x(pt) - this->y_z(pt) * this->z_xx(pt))
        + this->x_xz(pt) * (this->y_x(pt)  * this->z_y(pt) - this->y_y(pt) * this->z_x(pt))
        + this->x_z(pt)  * (this->y_xx(pt) * this->z_y(pt) + this->y_x(pt) * this->z_xy(pt)
                          - this->y_xy(pt) * this->z_x(pt) - this->y_y(pt) * this->z_xx(pt));                                            
}

template<>
double CellTransformation<double,3>::detJ_y ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }

  return  this->x_xy(pt) * (this->y_y(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_y(pt))
        + this->x_x(pt)  * (this->y_yy(pt) * this->z_z(pt) + this->y_y(pt) * this->z_yz(pt) 
                          - this->y_yz(pt) * this->z_y(pt) - this->y_z(pt) * this->z_yy(pt))
        - this->x_yy(pt) * (this->y_x(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
        - this->x_y(pt)  * (this->y_xy(pt) * this->z_z(pt) + this->y_x(pt) * this->z_yz(pt) 
                          - this->y_yz(pt) * this->z_x(pt) - this->y_z(pt) * this->z_xy(pt))
        + this->x_yz(pt) * (this->y_x(pt)  * this->z_y(pt) - this->y_y(pt) * this->z_x(pt))
        + this->x_z(pt)  * (this->y_xy(pt) * this->z_y(pt) + this->y_x(pt) * this->z_yy(pt)
                          - this->y_yy(pt) * this->z_x(pt) - this->y_y(pt) * this->z_xy(pt));                                            
}

template<>
double CellTransformation<double,3>::detJ_z ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }

  return  this->x_xz(pt) * (this->y_y(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_y(pt))
        + this->x_x(pt)  * (this->y_yz(pt) * this->z_z(pt) + this->y_y(pt) * this->z_zz(pt) 
                          - this->y_zz(pt) * this->z_y(pt) - this->y_z(pt) * this->z_yz(pt))
        - this->x_yz(pt) * (this->y_x(pt)  * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
        - this->x_y(pt)  * (this->y_xz(pt) * this->z_z(pt) + this->y_x(pt) * this->z_zz(pt) 
                          - this->y_zz(pt) * this->z_x(pt) - this->y_z(pt) * this->z_xz(pt))
        + this->x_zz(pt) * (this->y_x(pt)  * this->z_y(pt) - this->y_y(pt) * this->z_x(pt))
        + this->x_z(pt)  * (this->y_xz(pt) * this->z_y(pt) + this->y_x(pt) * this->z_yz(pt)
                          - this->y_yz(pt) * this->z_x(pt) - this->y_y(pt) * this->z_xz(pt));                                            
}

template<>
double CellTransformation<double,3>::detJ_xx ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return       this->x_xxx(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
        + 2. * this->x_xx(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_xz(pt) 
                                - this->y_xz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_xy(pt))
        +      this->x_x(pt)   * (this->y_xxy(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_xz(pt) 
                                + this->y_xy(pt)  * this->z_xz(pt) + this->y_y(pt)  * this->z_xxz(pt)
                                - this->y_xxz(pt) * this->z_y(pt)  - this->y_xz(pt) * this->z_xy(pt) 
                                - this->y_xz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xxy(pt))
        -      this->x_xxy(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
        - 2. * this->x_xy(pt)  * (this->y_xx(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_xz(pt) 
                                - this->y_xz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xx(pt))
        -      this->x_y(pt)   * (this->y_xxx(pt) * this->z_z(pt)  + this->y_xx(pt) * this->z_xz(pt) 
                                + this->y_xx(pt)  * this->z_xz(pt) + this->y_x(pt)  * this->z_xxz(pt)  
                                - this->y_xxz(pt) * this->z_x(pt)  - this->y_xz(pt) * this->z_xx(pt) 
                                - this->y_xz(pt)  * this->z_xx(pt) - this->y_z(pt)  * this->z_xxx(pt))
        +      this->x_xxz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
        + 2. * this->x_xz(pt)  * (this->y_xx(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_xy(pt)
                                - this->y_xy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xx(pt))
        +      this->x_z(pt)   * (this->y_xxx(pt) * this->z_y(pt)  + this->y_xx(pt) * this->z_xy(pt)
                                + this->y_xx(pt)  * this->z_xy(pt) + this->y_x(pt)  * this->z_xxy(pt)
                                - this->y_xxy(pt) * this->z_x(pt)  - this->y_xy(pt) * this->z_xx(pt)
                                - this->y_xy(pt)  * this->z_xx(pt) - this->y_y(pt)  * this->z_xxx(pt));
}

template<>
double CellTransformation<double,3>::detJ_xy ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return  this->x_xxy(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
        + this->x_xy(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt)
                           - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
        + this->x_xy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_xz(pt) 
                           - this->y_xz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_xy(pt))
        + this->x_x(pt)   * (this->y_xyy(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_yz(pt) 
                           + this->y_yy(pt)  * this->z_xz(pt) + this->y_y(pt)  * this->z_xyz(pt) 
                           - this->y_xyz(pt) * this->z_y(pt)  - this->y_xz(pt) * this->z_yy(pt) 
                           - this->y_yz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyy(pt))
        - this->x_xyy(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
        - this->x_xy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
        - this->x_yy(pt)  * (this->y_xx(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_xz(pt) 
                           - this->y_xz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xx(pt))        
        - this->x_y(pt)   * (this->y_xxy(pt) * this->z_z(pt)  + this->y_xx(pt) * this->z_yz(pt) 
                           + this->y_xy(pt)  * this->z_xz(pt) + this->y_x(pt)  * this->z_xyz(pt) 
                           - this->y_xyz(pt) * this->z_x(pt)  - this->y_xz(pt) * this->z_xy(pt) 
                           - this->y_yz(pt)  * this->z_xx(pt) - this->y_z(pt)  * this->z_xxy(pt))
        + this->x_xyz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
        + this->x_xz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt) 
                           - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt))
        + this->x_yz(pt)  * (this->y_xx(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_xy(pt)
                           - this->y_xy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xx(pt))
        + this->x_z(pt)   * (this->y_xxy(pt) * this->z_y(pt)  + this->y_xx(pt) * this->z_yy(pt) 
                           + this->y_xy(pt)  * this->z_xy(pt) + this->y_x(pt)  * this->z_xyy(pt)
                           - this->y_xyy(pt) * this->z_x(pt)  - this->y_xy(pt) * this->z_xy(pt) 
                           - this->y_yy(pt)  * this->z_xx(pt) - this->y_y(pt)  * this->z_xxy(pt));                                            
}

template<>
double CellTransformation<double,3>::detJ_xz ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return  this->x_xxz(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
        + this->x_xx(pt)  * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
        + this->x_xz(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_xz(pt) 
                           - this->y_xz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_xy(pt))
        + this->x_x(pt)   * (this->y_xyz(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_zz(pt) 
                           + this->y_yz(pt)  * this->z_xz(pt) + this->y_y(pt)  * this->z_xzz(pt) 
                           - this->y_xzz(pt) * this->z_y(pt)  - this->y_xz(pt) * this->z_yz(pt)
                           - this->y_zz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyz(pt))
        - this->x_xyz(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
        - this->x_xy(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_x(pt)) - this->y_z(pt)  * this->z_xz(pt))
        - this->x_yz(pt)  * (this->y_xx(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_xz(pt) 
                           - this->y_xz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xx(pt))
        - this->x_y(pt)   * (this->y_xxz(pt) * this->z_z(pt)  + this->y_xx(pt) * this->z_zz(pt) 
                           + this->y_xz(pt)  * this->z_xz(pt) + this->y_x(pt)  * this->z_xzz(pt) 
                           - this->y_xzz(pt) * this->z_x(pt)  - this->y_xz(pt) * this->z_xz(pt) 
                           - this->y_zz(pt)  * this->z_xx(pt) - this->y_z(pt)  * this->z_xxz(pt))
        + this->x_xzz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
        + this->x_xz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))
        + this->x_zz(pt)  * (this->y_xx(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_xy(pt)
                           - this->y_xy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xx(pt))                                            
        + this->x_z(pt)   * (this->y_xxz(pt) * this->z_y(pt)  + this->y_xx(pt) * this->z_yz(pt) 
                           + this->y_xz(pt)  * this->z_xy(pt) + this->y_x(pt)  * this->z_xyz(pt)
                           - this->y_xyz(pt) * this->z_x(pt)  - this->y_xy(pt) * this->z_xz(pt) 
                           - this->y_yz(pt)  * this->z_xx(pt) - this->y_y(pt)  * this->z_xxz(pt));                                            
}

template<>
double CellTransformation<double,3>::detJ_yy ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return  this->x_xyy(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
          this->x_xy(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
        + this->x_xy(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
        + this->x_x(pt)   * (this->y_yyy(pt) * this->z_z(pt)  + this->y_yy(pt) * this->z_yz(pt) 
                           + this->y_yy(pt)  * this->z_yz(pt) + this->y_y(pt)  * this->z_yyz(pt) 
                           - this->y_yyz(pt) * this->z_y(pt)  - this->y_yz(pt) * this->z_yy(pt) 
                           - this->y_yz(pt)  * this->z_yy(pt) - this->y_z(pt)  * this->z_yyy(pt))
        - this->x_yyy(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
        - this->x_yy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
        - this->x_yy(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
        - this->x_y(pt)   * (this->y_xyy(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_yz(pt) 
                           + this->y_xy(pt)  * this->z_yz(pt) + this->y_x(pt)  * this->z_yyz(pt) 
                           - this->y_yyz(pt) * this->z_x(pt)  - this->y_yz(pt) * this->z_xy(pt) 
                           - this->y_yz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyy(pt))
        + this->x_yyz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
        + this->x_yz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt) 
                           - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt))
        + this->x_yz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt)
                           - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt)) 
        + this->x_z(pt)   * (this->y_xyy(pt) * this->z_y(pt)  + this->y_xy(pt) * this->z_yy(pt) 
                           + this->y_xy(pt)  * this->z_yy(pt) + this->y_x(pt)  * this->z_yyy(pt)
                           - this->y_yyy(pt) * this->z_x(pt)  - this->y_yy(pt) * this->z_xy(pt) 
                           - this->y_yy(pt)  * this->z_xy(pt) - this->y_y(pt)  * this->z_xyy(pt)); 
}

template<>
double CellTransformation<double,3>::detJ_yz ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }

  return  this->x_xyz(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
        + this->x_xy(pt)  * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
        + this->x_xz(pt)  * (this->y_yy(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yy(pt))
        + this->x_x(pt)   * (this->y_yyz(pt) * this->z_z(pt)  + this->y_yy(pt) * this->z_zz(pt) 
                           + this->y_yz(pt)  * this->z_yz(pt) + this->y_y(pt)  * this->z_yzz(pt) 
                           - this->y_yzz(pt) * this->z_y(pt)  - this->y_yz(pt) * this->z_yz(pt) 
                           - this->y_zz(pt)  * this->z_yy(pt) - this->y_z(pt)  * this->z_yyz(pt))
        - this->x_yyz(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
        - this->x_yy(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
        - this->x_yz(pt)  * (this->y_xy(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xy(pt))
        - this->x_y(pt)   * (this->y_xyz(pt) * this->z_z(pt)  + this->y_xy(pt) * this->z_zz(pt) 
                           + this->y_xz(pt)  * this->z_yz(pt) + this->y_x(pt)  * this->z_yzz(pt) 
                           - this->y_yzz(pt) * this->z_x(pt)  - this->y_yz(pt) * this->z_xz(pt) 
                           - this->y_zz(pt)  * this->z_xy(pt) - this->y_z(pt)  * this->z_xyz(pt))
        + this->x_yzz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
        + this->x_yz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))
        + this->x_zz(pt)  * (this->y_xy(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yy(pt)
                           - this->y_yy(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xy(pt))                                            
        + this->x_z(pt)   * (this->y_xyz(pt) * this->z_y(pt)  + this->y_xy(pt) * this->z_yz(pt) 
                           + this->y_xz(pt)  * this->z_yy(pt) + this->y_x(pt)  * this->z_yyz(pt)
                           - this->y_yyz(pt) * this->z_x(pt)  - this->y_yy(pt) * this->z_xz(pt) 
                           - this->y_yz(pt)  * this->z_xy(pt) - this->y_y(pt)  * this->z_xyz(pt));
}

template<>
double CellTransformation<double,3>::detJ_zz ( const Vec<3,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }

  return  this->x_xzz(pt) * (this->y_y(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_y(pt))
        +  this->x_xz(pt) * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
        + this->x_xz(pt)  * (this->y_yz(pt)  * this->z_z(pt)  + this->y_y(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_y(pt)  - this->y_z(pt)  * this->z_yz(pt))
        + this->x_x(pt)   * (this->y_yzz(pt) * this->z_z(pt)  + this->y_yz(pt) * this->z_zz(pt) 
                           + this->y_yz(pt)  * this->z_zz(pt) + this->y_y(pt)  * this->z_zzz(pt) 
                           - this->y_zzz(pt) * this->z_y(pt)  - this->y_zz(pt) * this->z_yz(pt) 
                           - this->y_zz(pt)  * this->z_yz(pt) - this->y_z(pt)  * this->z_yzz(pt))
        - this->x_yzz(pt) * (this->y_x(pt)   * this->z_z(pt)  - this->y_z(pt)  * this->z_x(pt))
        - this->x_yz(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
        - this->x_yz(pt)  * (this->y_xz(pt)  * this->z_z(pt)  + this->y_x(pt)  * this->z_zz(pt) 
                           - this->y_zz(pt)  * this->z_x(pt)  - this->y_z(pt)  * this->z_xz(pt))
        - this->x_y(pt)   * (this->y_xzz(pt) * this->z_z(pt)  + this->y_xz(pt) * this->z_zz(pt)
                           + this->y_xz(pt)  * this->z_zz(pt) + this->y_x(pt)  * this->z_zzz(pt) 
                           - this->y_zzz(pt) * this->z_x(pt)  - this->y_zz(pt) * this->z_xz(pt) 
                           - this->y_zz(pt)  * this->z_xz(pt) - this->y_z(pt)  * this->z_xzz(pt))
        + this->x_zzz(pt) * (this->y_x(pt)   * this->z_y(pt)  - this->y_y(pt)  * this->z_x(pt))
        + this->x_zz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt) 
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))
        + this->x_zz(pt)  * (this->y_xz(pt)  * this->z_y(pt)  + this->y_x(pt)  * this->z_yz(pt)
                           - this->y_yz(pt)  * this->z_x(pt)  - this->y_y(pt)  * this->z_xz(pt))                                          
        + this->x_z(pt)   * (this->y_xzz(pt) * this->z_y(pt)  + this->y_xz(pt) * this->z_yz(pt) 
                           + this->y_xz(pt)  * this->z_yz(pt) + this->y_x(pt)  * this->z_yzz(pt)
                           - this->y_yzz(pt) * this->z_x(pt)  - this->y_yz(pt) * this->z_xz(pt) 
                           - this->y_yz(pt)  * this->z_xz(pt) - this->y_y(pt)  * this->z_xzz(pt));
}

template<>
double CellTransformation<double,2>::detJ_x ( const Vec<2,double>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  return this->x_xx(pt) * this->y_y(pt) + this->x_x(pt) * this->y_xy(pt)  
       - this->x_xy(pt) * this->y_x(pt) - this->x_y(pt) * this->y_xx(pt);                                          
}

template<>
double CellTransformation<double,2>::detJ_xx ( const Vec<2,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return  this->x_xxx(pt) * this->y_y(pt)  + this->x_xx(pt) * this->y_xy(pt) 
        + this->x_xx(pt)  * this->y_xy(pt) + this->x_x(pt)  * this->y_xxy(pt)  
        - this->x_xxy(pt) * this->y_x(pt)  - this->x_xy(pt) * this->y_xx(pt) 
        - this->x_xy(pt)  * this->y_xx(pt) - this->x_y(pt)  * this->y_xxx(pt) ;                                          
}

template<>
double CellTransformation<double,2>::detJ_xy ( const Vec<2,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return  this->x_xxy(pt) * this->y_y(pt)  + this->x_xx(pt) * this->y_yy(pt) 
        + this->x_xy(pt)  * this->y_xy(pt) + this->x_x(pt)  * this->y_xyy(pt)  
        - this->x_xyy(pt) * this->y_x(pt)  - this->x_xy(pt) * this->y_xy(pt) 
        - this->x_yy(pt)  * this->y_xx(pt) - this->x_y(pt)  * this->y_xxy(pt);                                          
}

template<>
double CellTransformation<double,2>::detJ_xz ( const Vec<2,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,2>::detJ_y ( const Vec<2,double>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  return this->x_xy(pt) * this->y_y(pt) + this->x_x(pt) * this->y_yy(pt)  
       - this->x_yy(pt) * this->y_x(pt) - this->x_y(pt) * this->y_xy(pt);                                          
}

template<>
double CellTransformation<double,2>::detJ_yy ( const Vec<2,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return this->x_xyy(pt) * this->y_y(pt)  + this->x_xy(pt) * this->y_yy(pt) 
       + this->x_xy(pt)  * this->y_yy(pt) + this->x_x(pt)  * this->y_yyy(pt)  
       - this->x_yyy(pt) * this->y_x(pt)  - this->x_yy(pt) * this->y_xy(pt) 
       - this->x_yy(pt)  * this->y_xy(pt) - this->x_y(pt)  * this->y_xyy(pt);                                          
}

template<>
double CellTransformation<double,2>::detJ_yz ( const Vec<2,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,2>::detJ_z ( const Vec<2,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,2>::detJ_zz ( const Vec<2,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,1>::detJ_x ( const Vec<1,double>& pt) const 
{
  if (this->order_ <= 1)
  {
    return 0.;
  }
  return this->x_xx(pt);                                          
}

template<>
double CellTransformation<double,1>::detJ_xx ( const Vec<1,double>& pt) const 
{
  if (this->order_ <= 2)
  {
    return 0.;
  }
  return  this->x_xxx(pt);                                          
}

template<>
double CellTransformation<double,1>::detJ_xy ( const Vec<1,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,1>::detJ_xz ( const Vec<1,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,1>::detJ_y ( const Vec<1,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,1>::detJ_yy ( const Vec<1,double>& pt) const 
{
  return 0.
}

template<>
double CellTransformation<double,1>::detJ_yz ( const Vec<1,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,1>::detJ_z ( const Vec<1,double>& pt) const 
{
  return 0.;
}

template<>
double CellTransformation<double,1>::detJ_zz ( const Vec<1,double>& pt) const 
{
  return 0.;
}
*/


// this->x_x(pt) * this->y_y(pt) - this->x_y(pt) * this->y_x(pt);


/*
         this->x_x(pt) * (this->y_y(pt) * this->z_z(pt) - this->y_z(pt) * this->z_y(pt)) 
       - this->x_y(pt) * (this->y_x(pt) * this->z_z(pt) - this->y_z(pt) * this->z_x(pt))
       + this->x_z(pt) * (this->y_x(pt) * this->z_y(pt) - this->y_y(pt) * this->z_x(pt));
*/

} // namespace doffem
} // namespace hiflow
