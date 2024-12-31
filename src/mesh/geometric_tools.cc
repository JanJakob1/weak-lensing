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

#include "common/bbox.h"
#include "common/data_tools.h"
#include "mesh/geometric_tools.h"
#include "fem/cell_trafo/cell_transformation.h"
#include "fem/cell_trafo/linearlinetransformation.h"
#include "fem/cell_trafo/lineartriangletransformation.h"
#include "fem/cell_trafo/lineartetrahedrontransformation.h"
#include "fem/cell_trafo/bilinearquadtransformation.h"
#include "fem/cell_trafo/linearpyramidtransformation.h"
#include "fem/cell_trafo/trilinearhexahedrontransformation.h"
#include "fem/cell_trafo/alignedhexahedrontransformation.h"
#include "fem/cell_trafo/alignedquadtransformation.h"
#include "mesh/boundary_domain_descriptor.h"
#include "mesh/periodicity_tools.h"
#include "mesh/entity.h"
#include "mesh/geometric_search.h"
#include "mesh/iterator.h"
#include "mesh/periodicity_tools.h"
#include "mesh/refinement.h"
#include "mesh/types.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/foreach.hpp>

namespace hiflow {
namespace mesh {

template <class DataType, int DIM> 
Vec<DIM, DataType> intersect_facet(const Vec<DIM, DataType> &point_a,
                                   const Vec<DIM, DataType> &point_b,
                                   const std::vector< DataType > &vertices,
                                   bool &success) {
  assert(!vertices.empty());
//  assert(!(vertices.size() % point_a.size()));

  Vec<DIM, DataType> intersection;
  const GDim gdim = DIM;
  assert((gdim == 2) || (gdim == 3));

  if (static_cast< int >(vertices.size()) != gdim * gdim) {
    // TODO: implementation for hexahedron facets
    NOT_YET_IMPLEMENTED;
    success = false;
    return intersection;
  }

  // implemented for a facet beeing a line in 2D and a triangle in 3D

  /*
   * 2D line:    3D triangle:    connection from point_a to point_b:
   *   B        C                       F = point_a
   *   *        *                      *
   *   |        |\                    /
   *   |        | \ E                / g
   *   |        |  \                *
   *   *        *---*              G = point_b
   *   A        A    B
   *
   *   E: x = A + x1(B-A) [+ x2(C-A)]  with xi >= 0 and sum xi <= 1
   *   g: x = F + x3(G-F)              with 0 <= x3 <= 1
   *   equating yields GLS:
   *   x1(B-A) [+ x2(C-A)] + x3(F-G) = F-A
   */

  // compute intersection
  Mat<DIM, DIM, DataType> mat;
  Vec<DIM, DataType> vec;
  for (int d = 0; d < DIM; ++d) {
    for (int dir = 0; dir < DIM - 1; ++dir) {
      mat(d,dir) = vertices[(dir + 1) * DIM + d] - vertices[d];
    }
    mat(d,DIM - 1) = point_a[d] - point_b[d];
    vec[d] = point_a[d] - vertices[d];
  }

  // If the system is not solvable, line and facet are parallel
  const bool solved = gauss(mat, vec);
  if (!solved) {
    // check parallelism
    std::vector<Vec<DIM, DataType> > directions (gdim - 1);
    for (int d = 0; d < gdim; ++d) {
      for (int dir = 0; dir < gdim - 1; ++dir) {
        directions[dir][d] =
            vertices[(dir + 1) * gdim + d] - vertices[d];
      }
    }
    Vec<DIM, DataType> facet_normal = normal<DataType, DIM>(directions);
    Vec<DIM, DataType> dir_a_b(point_b);
    dir_a_b.Axpy(point_a, -1);
    if (std::abs(dot(dir_a_b, facet_normal)) < GEOM_TOL) {
      // vectors are parallel
      // TODO: check intersection in this case
      NOT_YET_IMPLEMENTED;
    }
    success = false;
    return intersection;
  }

  // the facet is intersected if
  // 0 <= x3 <= 1
  // xi >= 0 and sum xi <= 1
  if (vec[gdim - 1] < -GEOM_TOL || vec[gdim - 1] > 1 + GEOM_TOL) {
    success = false;
    return intersection;
  }
  DataType sum = 0;
  for (int d = 0; d < gdim - 1; ++d) {
    if (vec[d] < -GEOM_TOL) {
      success = false;
      return intersection;
    }
    sum += vec[d];
  }
  if (sum > 1 + GEOM_TOL) {
    success = false;
    return intersection;
  }

  // fill intersection coordinate vector
  for (int d = 0; d < gdim; ++d) {
    intersection[d] = point_a[d] + vec[gdim - 1] * (point_b[d] - point_a[d]);
  }

  success = true;
  return intersection;
}

template <class DataType, int DIM> 
bool crossed_facet(const Vec<DIM, DataType> &point_a,
                   const Vec<DIM, DataType> &point_b,
                   const std::vector< DataType > &vertices) {

  assert(!vertices.empty());
  assert(!(vertices.size() % point_a.size()));

  bool success;
  intersect_facet<DataType, DIM>(point_a, point_b, vertices, success);
  return success;
}

template <class DataType, int DIM>
bool point_inside_entity(const Vec<DIM, DataType> &point, 
                         const int tdim,
                         const std::vector< DataType > &vertices) {
  // implemented for lines,
  //                 triangles,
  //                 quadrilaterals in one plane,
  //                 tetrahedrons and
  //                 hexahedrons
  // in the dimensions 2D and 3D
  assert(tdim == 1 || tdim == 2 || tdim == 3);
  assert(!vertices.empty());
  assert(!(vertices.size() % DIM));

  GDim gdim = DIM;
  assert(gdim == 2 || gdim == 3);
  bool inside = false;

  switch (tdim) {
    // lines, rectangles and triangles can be handled via distance / volume
    // computation: if vol(ABCD) == vol(PAB) + vol(PBC) + vol(PCD) + vol(PDA)
    // the point lies in the entity
  case 1: {
    assert(gdim >= 1);
    assert(static_cast< int >(vertices.size()) == gdim * 2);
    Vec<DIM, DataType> p_a(vertices, 0);
    Vec<DIM, DataType> p_b(vertices, gdim);


    DataType point_to_a = norm(point - p_a);
    DataType point_to_b = norm(point - p_b);
    DataType a_to_b = norm(p_a - p_b);
    return (point_to_a + point_to_b < a_to_b + GEOM_TOL);
    break;
  }
  case 2: {
    assert(gdim >= 2);

    DataType area_sum = 0;
    const int num_vertices = vertices.size() / gdim;
    // Attention: Here we assume that all vertices lie in one plane!
    assert((gdim < 3 || num_vertices < 4 ||
           vertices_inside_one_hyperplane<DataType, DIM>(vertices, gdim)));
    for (int i = 0; i < num_vertices; ++i) {
      std::vector< DataType > tri_verts(DIM);
      if (point[0] == vertices[i * gdim] &&
          point[1] == vertices[i * gdim + 1]) {
        return true;
      }
      for (int d=0; d<DIM; ++d)
      {
        tri_verts[d] = point[d];
      }

      tri_verts.insert(tri_verts.end(), vertices.begin() + i * gdim,
                       vertices.begin() + (i + 1) * gdim);
      tri_verts.insert(tri_verts.end(),
                       vertices.begin() + ((i + 1) % num_vertices) * gdim,
                       vertices.begin() + ((i + 1) % num_vertices + 1) * gdim);
      area_sum += triangle_area<DataType, DIM>(tri_verts);
    }
    DataType entitity_area = 0;
    for (int i = 1; i < num_vertices - 1; ++i) {
      std::vector< DataType > tri_verts(vertices.begin(),
                                        vertices.begin() + gdim);
      tri_verts.insert(tri_verts.end(), vertices.begin() + i * gdim,
                       vertices.begin() + (i + 2) * gdim);
      entitity_area += triangle_area<DataType, DIM>(tri_verts);
    }
    return (area_sum < entitity_area + GEOM_TOL);
    break;
  }
  case 3: {
    assert(gdim == 3);
    if (static_cast< int >(vertices.size()) == (gdim * (gdim + 1))) {
      // Tetrahedron
      // Parametric equation to check, where the point is:
      // P = A + x0(B-A) [ + x1(C-A) ( + x2(D-A) ) ]
      // This algorithm could also handle lines in 1D and triangles in 2D
      Mat<DIM, DIM, DataType > mat;
      Vec<DIM, DataType > vec;
      for (int i = 0; i < gdim; ++i) {
        for (int j = 0; j < gdim; ++j) {
          mat(i,j) = vertices[(j + 1) * gdim + i] - vertices[i];
        }
        vec[i] = point[i] - vertices[i];
      }
      // solve this linear system of equations
      const bool solved = gauss(mat, vec);
      // if the system is not solvable, the cell is degenerated
      assert(solved);

      // the point lies in the cell, if
      // xi >= 0
      // sum xi <= 1
      DataType sum = 0;
      for (int d = 0; d < gdim; ++d) {
        if (vec[d] < -GEOM_TOL) {
          return false;
        }
        sum += vec[d];
      }
      return (sum < 1 + GEOM_TOL);
    } else if (static_cast< int >(vertices.size()) == (gdim * 8)) {
      // TODO: implementation of a more performant solutions for hexahedrons
      // For an arbitrary hex, the four vertices of one facet may not lie in one
      // plane.
      doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
        = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellHexStd<DataType, DIM> );
        
      hiflow::doffem::TriLinearHexahedronTransformation< DataType, DIM > ht(ref_cell);
      ht.reinit(vertices);
      Vec<DIM, DataType> ref_pt;
      return (ht.contains_physical_point(point, ref_pt));
    } else {
      // Pyramid
      doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
        = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellPyrStd<DataType, DIM> );
        
      hiflow::doffem::LinearPyramidTransformation< DataType, DIM > ht(ref_cell);
      ht.reinit(vertices);
      Vec<DIM, DataType> ref_pt;
      return (ht.contains_physical_point(point, ref_pt));
    }
    break;
  }
  default:
    assert(0);
    return false;
    break;
  }
  // A return should have been called before.
  assert(0);
  return inside;
}

template <class DataType, int DIM>
bool point_inside_cell(const Vec<DIM, DataType> &point, 
                       const std::vector< DataType > &vertices,
                       Vec<DIM, DataType> &ref_point) 
{
  // implemented for lines,
  //                 triangles,
  //                 quadrilaterals,
  //                 tetrahedrons,
  //                 hexahedrons,
  //                 pyramids
  // in the dimensions 2D and 3D
  int tdim = DIM;
  GDim gdim = DIM;

  assert(!vertices.empty());
  assert(!(vertices.size() % DIM));
  const int num_vertices = static_cast< int >(vertices.size() / gdim);

  //assert(gdim == 2 || gdim == 3);
  bool inside = false;
  
  switch (tdim) 
  {
    case 1: 
    {
      // line
      assert(gdim == 1);
      assert(num_vertices == 2);

      doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
        = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellLineStd<DataType, DIM> );
      hiflow::doffem::LinearLineTransformation< DataType, DIM > ht(ref_cell);
      ht.reinit(vertices);
      bool found = (ht.contains_physical_point(point, ref_point));
      return found;  
      break;
    }
    case 2: 
    {
      assert(gdim == 2);
      assert(num_vertices == 3 || num_vertices == 4);

      if (num_vertices == 3)
      {
        // triangle
        doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
          = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellTriStd<DataType, DIM> );
        
        hiflow::doffem::LinearTriangleTransformation< DataType, DIM > ht(ref_cell);
        ht.reinit(vertices);
        bool found = (ht.contains_physical_point(point, ref_point));
        return found;  
      }
      else 
      {
        if (is_parallelogram(vertices))
        {
          // parallelogram
          doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
            = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellQuadStd<DataType, DIM> );
        
          hiflow::doffem::AlignedQuadTransformation< DataType, DIM > ht(ref_cell);
          ht.reinit(vertices);
          bool found = (ht.contains_physical_point(point, ref_point));
          return found;  
        }
        else 
        {
          // general quadrilateral
          doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
            = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellQuadStd<DataType, DIM> );
        
          hiflow::doffem::BiLinearQuadTransformation< DataType, DIM > ht(ref_cell);
          ht.reinit(vertices);
          bool found = (ht.contains_physical_point(point, ref_point));
          return found;  
        }
      }
      break;
    }
    case 3: 
    {
      assert(gdim == 3);
      assert(num_vertices == 4 || num_vertices == 5 || num_vertices == 8);
        
      if (num_vertices == 4)
      {
        // Tetrahedron
        doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
          = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellTetStd<DataType, DIM> );
        
        hiflow::doffem::LinearTetrahedronTransformation< DataType, DIM > ht(ref_cell);
        ht.reinit(vertices);
        bool found = (ht.contains_physical_point(point, ref_point));
        return found;  
      } 
      else if (num_vertices == 8)
      {
        // Hexahedron
        if (is_parallelepiped (vertices))
        {
          doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
            = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellHexStd<DataType, DIM> );
        
          hiflow::doffem::AlignedHexahedronTransformation< DataType, DIM > ht(ref_cell);
          ht.reinit(vertices);
          bool found = (ht.contains_physical_point(point, ref_point));
          return found;  
        }
        else
        {
          doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
            = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellHexStd<DataType, DIM> );
        
          hiflow::doffem::TriLinearHexahedronTransformation< DataType, DIM > ht(ref_cell);
          ht.reinit(vertices);
          bool found = (ht.contains_physical_point(point, ref_point));
          return found;  
        }
        return false;
      } 
      else 
      {
        // Pyramid
        doffem::ConstRefCellPtr<DataType, DIM> ref_cell 
        = doffem::ConstRefCellPtr<DataType, DIM>(new doffem::RefCellPyrStd<DataType, DIM> );
        
        hiflow::doffem::LinearPyramidTransformation< DataType, DIM > ht(ref_cell);
        ht.reinit(vertices);
        bool found = (ht.contains_physical_point(point, ref_point));
        return found;  
      }
      break;
    }
    default:
      assert(0);
      return false;
      break;
  }
  // A return should have been called before.
  assert(0);
  return inside;
}

template <class DataType, int DIM>
bool vertices_inside_one_hyperplane(const std::vector< DataType > &vertices,
                                    const GDim gdim) {
  assert(!vertices.empty());
  assert(!(vertices.size() % gdim));
  int num_points = vertices.size() / gdim;
  assert(num_points > gdim);

  std::vector< Vec<DIM, DataType> > directions(DIM-1);
  for (int i = 0; i < DIM - 1; ++i) {
    for (GDim d = 0; d < DIM; ++d) {
      directions[i][d] = vertices[d] - vertices[(i + 1) * DIM + d];
    }
  }

  Vec<DIM, DataType> origin(vertices, 0);
  Vec<DIM, DataType> plane_normal = normal<DataType, DIM> (directions);
  for (int n = gdim; n < num_points; ++n) {
    Vec<DIM, DataType> test_point(vertices, n * gdim);
    if (!in_plane<DataType, DIM>(test_point, origin, plane_normal, GEOM_TOL)) {
      return false;
    }
  }
  // at this point, it is proved, that all points lie in one plane.
  return true;
}

template <class DataType, int DIM>
DataType triangle_area(const std::vector< DataType > &vertices) {
  assert(!vertices.empty());
  assert(!(vertices.size() % 3));
  const GDim gdim = vertices.size() / 3;
  assert(gdim > 1);

  // 0.5 * |A-B|*dist(C,AB)
  Vec<DIM, DataType> p_a(vertices, 0);
  Vec<DIM, DataType> p_b(vertices, DIM);
  Vec<DIM, DataType> p_c(vertices, 2*DIM);
  Vec<DIM, DataType> dir_a_b(p_a);
  DataType dir_a_b_norm = 0;
  for (int i = 0; i < gdim; ++i) {
    dir_a_b[i] -= p_b[i];
    dir_a_b_norm += dir_a_b[i] * dir_a_b[i];
  }
  dir_a_b_norm = std::sqrt(dir_a_b_norm);

  return 0.5 * dir_a_b_norm * distance_point_line<DataType, DIM>(p_c, p_a, dir_a_b);
}

template <class DataType, int DIM>
DataType facet_area(const std::vector< DataType > &vertices,
                    const GDim gdim) {
  assert(!vertices.empty());
  assert(gdim == 2 || gdim == 3);
  assert(!(vertices.size() % gdim));
  Coordinate area = 0;
  switch (gdim) {
  case (2): {
    // in 2D: Facet = Line
    assert(vertices.size() == 4);
    Vec<DIM, DataType > vert_a(vertices, 0);
    Vec<DIM, DataType > vert_b(vertices, gdim);
    area = norm(vert_a - vert_b);
  } break;
  case (3): {
    // The vertices have to lie in one plane.
    const int num_vertices = vertices.size() / gdim;
    assert( ((num_vertices < 4) || vertices_inside_one_hyperplane<DataType, DIM>(vertices, gdim) ));
    for (int i = 1; i < num_vertices - 1; ++i) {
      std::vector< DataType > tri_verts(vertices.begin(),
                                        vertices.begin() + gdim);
      tri_verts.insert(tri_verts.end(), vertices.begin() + i * gdim,
                       vertices.begin() + (i + 2) * gdim);
      area += triangle_area<DataType, DIM>(tri_verts);
    }
  } break;
  default:
    NOT_YET_IMPLEMENTED;
    break;
  }
  return area;
}

template <class DataType, int DIM>
bool in_plane(const Vec<DIM, DataType> &point,
              const Vec<DIM, DataType> &origin,
              const Vec<DIM, DataType> &normal, const DataType eps) {

  const DataType distance = distance_point_hyperplane<DataType, DIM>(point, origin, normal);
  return (std::abs(distance) < eps);
}

template <class DataType, int DIM>
bool crossed_plane(const Vec<DIM, DataType> &point_a,
                   const Vec<DIM, DataType> &point_b,
                   const Vec<DIM, DataType> &origin,
                   const Vec<DIM, DataType> &normal) {
  
  const DataType distance_a = distance_point_hyperplane<DataType, DIM>(point_a, origin, normal);
  const DataType distance_b = distance_point_hyperplane<DataType, DIM>(point_b, origin, normal);
  return (distance_a * distance_b < 0);
}

template <class DataType, int DIM>
DataType distance_point_hyperplane(const Vec<DIM, DataType> &point,
                                   const Vec<DIM, DataType> &origin,
                                   const Vec<DIM, DataType> &normal) {
  return dot(point, normal) - dot(origin, normal);
}

template <class DataType, int DIM>
DataType distance_point_line(const Vec<DIM, DataType> &point,
                             const Vec<DIM, DataType> &origin,
                             const Vec<DIM, DataType> &direction) {

  Vec<DIM, DataType> foot = foot_point_line<DataType, DIM>(point, origin, direction);
  return norm(point - foot);
}

template <class DataType, int DIM>
Vec<DIM, DataType> foot_point_hyperplane(const Vec<DIM, DataType> &point,
                      const Vec<DIM, DataType> &origin,
                      const Vec<DIM, DataType> &normal) {
  // orthogonal projection
  // f = p - n * (p-o) dot n / (n dot n)
  Vec<DIM, DataType> tmp(point);
  tmp.Axpy(origin, -1);
  DataType factor = dot(tmp, normal) / dot(normal, normal);
  tmp = point;
  tmp.Axpy(normal, -factor);
  return tmp; 
}

template <class DataType, int DIM>
Vec<DIM, DataType> foot_point_line(const Vec<DIM, DataType> &point,
                const Vec<DIM, DataType> &origin,
                const Vec<DIM, DataType> &direction) {

  // foot of point on line
  // foot = o + d * ((p-o) dot d) / (d dot d)
  Vec<DIM, DataType> tmp(point);
  tmp.Axpy(origin, -1);
  DataType factor = dot(tmp, direction) / dot(direction, direction);

  tmp = origin;
  tmp.Axpy(direction, factor);
  return tmp; 
}

template <class DataType, int DIM>
Vec<DIM, DataType> normal(const std::vector< Vec<DIM, DataType> > &directions) {
  assert(!directions.empty());

  Vec<DIM, DataType> result;
  if (DIM == 2) {
    assert(directions.size() == 1);

    DataType dir_norm = norm(directions[0]);
    assert(dir_norm > 0);
    result[0] = directions[0][1] / dir_norm;
    result[1] = -directions[0][0] / dir_norm;
  } 
  else if (DIM == 3) {
    assert(directions.size() == 2);
    result = cross(directions[0], directions[1]);
    DataType res_norm = norm(result);
    assert(res_norm > 0);
    result *= (1. / res_norm);
  } else {
    // directions has the wrong size
    assert(0);
  }
  return result;
}

template <class DataType, int DIM>
DataType distance_point_facet(const Vec<DIM, DataType> &point,
                              const std::vector< DataType > &vertices,
                              Vec<DIM, DataType> &closest_point) {
  assert(!vertices.empty());
  assert(!(vertices.size() % point.size()));

  // Working for 3d and 2d
  GDim gdim = DIM;
  closest_point.Zeros();

  int num_verts = vertices.size() / gdim;

  switch (gdim) {
  case 1: {
    NOT_YET_IMPLEMENTED;
  } break;
  case 2: {
    // 2D: Distance to a finite line.
    // Check first if the closest point is the orthogonal projection on the
    // facet.
    assert(num_verts == 2);
    Vec<DIM, DataType> origin(vertices, 0);
    Vec<DIM, DataType> direction;
    for (int i = 0; i < gdim; i++) {
      direction[i] = vertices[i + gdim] - vertices[i];
    }

    Vec<DIM, DataType> foot_point = foot_point_line<DataType, DIM> (point, origin, direction);
    
    if (point_inside_entity<DataType, DIM>(foot_point, 1, vertices)) {
      closest_point = foot_point;
      return norm(closest_point - point);
    }
  } break;
  case 3: {
    // 3D: first we do a orthogonal projection onto the plane spanned
    // by the facet. The projected point is called foot_point.
    assert(num_verts == 3 || num_verts == 4);
    Vec<DIM, DataType> origin(vertices, 0);

    std::vector<Vec<DIM, DataType> >directions (2);
    for (int i = 0; i < 2; ++i) 
    {
      for (int d=0; d<DIM; ++d) 
      {
        directions[i][d] = vertices[(i+1)*DIM+d] - vertices[d];
      }
    }

    Vec<DIM, DataType> facet_normal = normal<DataType, DIM> (directions);
    // TODO: The vertices of the facet currently have to lie in one plane.
    assert((num_verts < 4 || vertices_inside_one_hyperplane<DataType, DIM>(vertices, gdim)));
    Vec<DIM, DataType> foot_point =
        foot_point_hyperplane<DataType, DIM> (point, origin, facet_normal);
    // if the foot_point is inside the entity we are done and return
    // the foot_point as closest_point
    if (point_inside_entity<DataType, DIM>(foot_point, 2, vertices)) {
      closest_point = foot_point;
      return norm(closest_point - point);
    }

    // else we project our point onto the subentities (lines)
    DataType dist = std::numeric_limits< DataType >::max();
    bool found = false;
    for (int j = 0; j < num_verts; j++) {
      std::vector< DataType > line(2 * DIM);
      Vec<DIM, DataType> direction;
      for (int i = 0; i < DIM; i++) {
        line[i] = vertices[i + j * gdim];
        origin[i] = line[i];
        line[i + gdim] = vertices[i + ((j + 1) % num_verts) * DIM];
        direction[i] = line[i + gdim] - line[i];
      }

      foot_point = foot_point_line<DataType, DIM>(point, origin, direction);
      // if the projected point is inside the entity we are done!
      if (point_inside_entity<DataType, DIM>(foot_point, 1, line)) {
        DataType tmp_dist = norm (foot_point - point);
        if (tmp_dist < dist) {
          dist = tmp_dist;
          closest_point = foot_point;
          found = true;
        }
      }
    }
    if (found) {
      return dist;
    }
  } break;
  }
  // If the closest point is not an orthogonal projection of the point
  // on the facet (or its edges in 3D), one of the vertices is the
  // closest point
  DataType dist = std::numeric_limits< DataType >::max();
  for (int n = 0; n < num_verts; n++) {
    Vec<DIM, DataType> current_vertex (vertices, n * gdim);

    DataType temp_dist = norm (current_vertex - point);
    if (temp_dist < dist) {
      dist = temp_dist;
      closest_point = current_vertex;
    }
  }
  return dist;
}

template <class DataType, int DIM>
Vec<DIM, DataType> project_point(const BoundaryDomainDescriptor<DataType, DIM> &bdd,
                                 const Vec<DIM, DataType> &p,
                                 const MaterialNumber mat) {
  LOG_DEBUG(2, "Starting vector: " << string_from_range(p.begin(), p.end()));
  // starting vector is p
  Vec<DIM, DataType> xk = p;
  int num_it = 0;

  // step_length is used as a stopping criteria. Initially chosen 1.
  // to not fullfill the criteria
  DataType step_length = 1.;
  // if steplength < TOL the algorithm will stop and it is assumed that
  // a solution has been found.
  const DataType TOL = 1e-8;
  // maximal number of iteration
  const int max_it = 1e3;

  // leave early if initial point is already a zero of g.
  {
    const DataType init_g_xk = bdd.eval_func(xk, mat);
    const DataType ABS_TOL = 1e-8;
    if (std::abs(init_g_xk) < ABS_TOL) {
      LOG_DEBUG(2, "Left early because point already on boundary");
      return xk;
    }
  }

  // Following does solve the problem:
  //  min_{x \in R^3} 1/2|| p - x ||^2
  //  s.t.  g(x) = 0
  // Where g is given through the BoundaryDomainDescriptor.
  // The algorithm is based on an SQP approach:
  //
  // x_0 starting point
  // for k=0 to ...
  //   min_{d_k \in R^3} 1/2|| p -(x_k + d_k)||^2
  //   s.t.  g(x_k) + d_k*grad(g)(x_k) = 0
  //
  //   x_{k+1} = x_k + d_k
  // endfor
  //
  // The solution of the minimizing problem is done via the Lagrange
  // function: L(d, l) = 1/2|| p -(x_k + d_k)||^2 + l*(g(x_k)+d*grad(g)(x_k))
  // Need to solve grad(L) = 0.
  // This can be done on an algebraic level to get an "exact" solution.

  while ((step_length > TOL) && (num_it < max_it)) {
    const DataType g_xk = bdd.eval_func(xk, mat);
    const Vec<DIM, DataType> grad_xk = bdd.eval_grad(xk, mat);
    Vec<DIM, DataType> pmxk = p - xk; 

    const DataType grad_g_dot_pmxk = dot(grad_xk, pmxk);

    // lambda_k = (grad(g)(x_k)*(p - x_k) + g(x_k))/||grad(g)(x_k)||^2
    const DataType lambdak = (grad_g_dot_pmxk + g_xk) / dot(grad_xk, grad_xk);

    // d_k = p - x_k - lambda_k*grad(x_k)
    Vec<DIM, DataType> dk(pmxk);
    dk.Axpy(grad_xk, -lambdak);

    // damping?
    DataType damping_factor = 1.;
    // x_{k+1} = x_k + d_k
    xk.Axpy(dk, damping_factor);

    step_length = sqrt(norm(dk));
    ++num_it;

    // Some high level debug information.
    LOG_DEBUG(99, "lambdak: " << lambdak);
//    LOG_DEBUG(99, "dk: " << string_from_range(dk.begin(), dk.end()));
//    LOG_DEBUG(99, "xk(updated): " << string_from_range(xk.begin(), xk.end()));
    LOG_DEBUG(99, "g(xk): " << g_xk);
//    LOG_DEBUG(99, "grad_g(xk): " << string_from_range(grad_xk.begin(),
//                                                      grad_xk.end()));
    LOG_DEBUG(99, "steplength: " << step_length);
  }
  if (num_it == max_it)
  {
    LOG_DEBUG(2, "Stopped Iteration because max Iteration was reached");
  }
//  LOG_DEBUG(2, "Final vector: " << string_from_range(xk.begin(), xk.end()));
  LOG_DEBUG(2, "Final defect: " << bdd.eval_func(xk));
  LOG_DEBUG(2, "Number of Iterations needed: " << num_it);

  return xk;
}

template < class DataType, int DIM >
bool is_point_on_subentity(const Vec<DIM, DataType>& point, const std::vector< Vec<DIM, DataType> > &points) 
{
  int verbatim_mode = 0; // 0 -> no output
  // 1 -> some output

  if (verbatim_mode > 0) {
    std::cout << "Is point (";
    for (size_t k = 0; k < point.size(); ++k) {
      std::cout << point[k] << " ";
    }
    std::cout << ") on subentity?\n";
  }

  assert(!points.empty());

  DataType eps = 1.0e-6;
  DataType delta = 1.0e-6;
  if (typeid(DataType) == typeid(double)) {
    delta = 1.0e-13;
  }

  DataType dist = 1.0; // value to be calculated

  if (points.size() == 1) 
  {
    // single point
    dist = norm(point - points[0]);
  } 
  else if (points.size() == 2) 
  {
    // line

    // distance between point P and even (defined by point A and B)

    Vec<DIM, DataType> AP (point);
    AP -= points[0];

    DataType norm_AP = norm(AP);

    Vec<DIM,DataType> AB (points[1]);
    AB -= points[0];
    
    DataType norm_AB = norm(AB);

    DataType AP_times_AB = dot(AP, AB);

    DataType tmp = norm_AP * norm_AP - (AP_times_AB * AP_times_AB) / (norm_AB * norm_AB);

    tmp = std::abs(tmp);
    assert(tmp >= 0.0 || tmp <= eps);
    dist = std::sqrt(tmp);

    if (verbatim_mode > 0) {
      std::cout << "dist to line = " << dist << "\n";
    }

    if (dist < delta) {

      // further test: point between the two defining points of the line?

      // Compute alpha such that P = A + alpha * AB . Then
      // alpha == 0 => P = A
      // alpha == 1 => P = B
      // 0 < alpha < 1 => P lies between A and B
      // otherwise P lies outside A and B .
      const DataType alpha = dot(AB, AP) / dot(AB, AB);

      if (alpha < -eps || alpha > 1 + eps) {
        return false;
      }
    }

  } else {
    // face

    assert(DIM == 3);

    // 1. calculate normed normal vector of face (p0, p1, p2)

    Vec<DIM,DataType> v(points[1]); // v := points.at(1) - points.at(0)
    v -= points[0];
    
    Vec<DIM,DataType> w(points[2]); // w := points.at(2) - points.at(0)
    w -= points[0];

    Vec<DIM,DataType> normal; // cross product
    normal[0] = (v[1] * w[2] - v[2] * w[1]);
    normal[1] = (v[2] * w[0] - v[0] * w[2]);
    normal[2] = (v[0] * w[1] - v[1] * w[0]);

    // normalize normal vector
    DataType norm_normal = norm(normal);
    
    normal[0] = normal[0] / norm_normal;
    normal[1] = normal[1] / norm_normal;
    normal[2] = normal[2] / norm_normal;

    // 2. calculate parameter d

    DataType d = dot(points[0], normal);

    // 3. calculate value for considered point

    DataType d_point = dot(point, normal);

    // 4. calculate distance point to face

    dist = d - d_point;
    dist = std::abs(dist);

    // 5. in case of more than 3 points

    for (size_t i = 3, e_i = points.size(); i < e_i; ++i) {
      DataType d_temp = dot(points[i], normal);
      assert(std::abs(d - d_temp) < eps);
    }

    // 6. for refined cells

    // until now it has been checked that the DoF point is on the plane
    // defined by the given points, now it is further checked whether
    // the DoF point lies within the convex envelope

    /// Idea of point on face test:
    /// Calculate the normal vectors by sucessive cross products
    /// \f \vec{n}_i(p_{i}-p_{i-1})\times (p-p_{i-1})\quad i=1\cdots n_points \f
    /// and check whether these normal vectors are parallel or zero.

    /// TODO (staffan): an alternative, but similar algorithm that one might try
    /// out would be as follows:
    ///
    /// For each edge compute edge vector e_i = p_i - p_{i-1} and
    /// point vector v_i = p - p_{i-1} (like above). Then compute a
    /// vector u_i = cross(e_i, v_i) and the normal to a plane
    /// containing e_i: n_i = cross(e_i, u_i) . Finally check the sign
    /// of dot(v_i, e_i) to find out what side of the plane that p
    /// lies on. This should remain the same for all edges.
    ///
    /// A special case is when p lies on the edge in question. As in
    /// the current algorithm, this would have to be tested for
    /// separately. Perhaps both algorithms actually work the same?

    if (dist < eps) {
      // normal vector
      Vec< DIM, DataType > vec_normal;
      bool first_time = true;

      for (size_t i = 1, e_i = points.size(); i < e_i; ++i) {
        // set vectors

        Vec< DIM, DataType > vec_edge(points[i]); // p_i - p_{i-1}
        vec_edge -= points[i-1];

        Vec< DIM, DataType > vec_to_point(point); // point - p_{i-1}
        vec_to_point -= points[i-1];

        // cross product
        Vec< DIM, DataType > vec_normal_temp;
        vec_normal_temp = cross(vec_edge, vec_to_point);
        DataType norm = std::abs(hiflow::norm(vec_normal_temp));
        if (verbatim_mode > 0) {
          std::cout << "vec_edge = " << vec_edge << ", "
                    << "\nvec_to_point = " << vec_to_point
                    << ",\nvec_normal_temp = " << vec_normal_temp << "\n";
          std::cout << "norm = " << norm << "\n";
        }
        if (norm > delta) {
          vec_normal_temp =
              vec_normal_temp * (1. / hiflow::norm(vec_normal_temp));
        }
        if (norm > delta) {
          if (first_time) {
            vec_normal = vec_normal_temp;
            first_time = false;
          } else {
            DataType diff = 0;
            for (size_t d = 0; d < DIM; ++d) {
              diff += std::abs(vec_normal[d] - vec_normal_temp[d]);
            }

            if (diff > delta) {
              dist = 1.0; // not within the convex envelope
            }
          }
        }
      }
    } else {
      if (verbatim_mode > 0) {
        std::cout << "dist > eps: dist = " << dist << ", eps = " << eps << "\n";
      }
    }
  }
  return dist < eps;
}

/// find all subentities on which a given point lies
template < class DataType, int DIM > 
void find_subentities_containing_point (const Vec<DIM, DataType>& pt, 
                                        const mesh::CellType *ref_cell,
                                        const std::vector< Vec<DIM, DataType> >& coord_vertices, 
                                        std::vector< std::vector<int> > &dof_on_subentity)
{
  const int verbatim_mode = 0; 
  
  dof_on_subentity.clear();
  dof_on_subentity.resize(ref_cell->tdim()+1);

  // trival case tdim (=cell)
  dof_on_subentity[ref_cell->tdim()].push_back(0);

  // get point coordinates for this entity, including refined entities
  const int gdim = DIM;
  std::vector< double > cv; // coord_vertices_ in other numeration
  for (size_t p = 0, e_p = coord_vertices.size(); p != e_p; ++p) 
  {
    for (size_t c = 0, e_c = gdim; c != e_c; ++c)
    {
      cv.push_back( static_cast<double>(coord_vertices[p][c]) );
    }
  }

  std::vector< double > coord_vertices_refined;
  mesh::compute_refined_vertices(*ref_cell, gdim, cv, coord_vertices_refined);

  if (verbatim_mode > 0) 
  {
    for (size_t i = 0, e_i = coord_vertices_refined.size(); i != e_i; ++i) 
    {
      std::cout << "\t" << coord_vertices_refined[i] << std::endl;
    }
  }

  // insert filtered DoFs
  /// \todo also tdim = ref_cell_->tdim() should be available
  /// \todo coord_ could also be stored by coord_on_subentity_

  for (size_t tdim = 0, e_tdim = ref_cell->tdim(); tdim != e_tdim; ++tdim) 
  {
    // for output purposes
    std::string entity_type;
    if (tdim == 0) 
    {
      entity_type = "point";
    }
    if (tdim == 1) 
    {
      entity_type = "line";
    }
    if (tdim == 2) 
    {
      entity_type = "face";
    }
    if (tdim == 3) 
    {
      entity_type = "volume";
    }

    for (size_t idx_entity = 0, e_idx_entity = ref_cell->num_entities(tdim);
         idx_entity != e_idx_entity; ++idx_entity) 
    {
      if (verbatim_mode > 0) 
      {
        std::cout << "DoF points on " << entity_type << " " << idx_entity << ":" << std::endl;
      }

      // get point coordinates for this entity
      std::vector< Vec<DIM, DataType> > points;
      std::vector< int > vertex = ref_cell->local_vertices_of_entity(tdim, idx_entity);

      // points can't be handled using local_vertices_of_entity()
      if (tdim == 0) 
      {
        vertex.clear();
        vertex.push_back(idx_entity);
      }

      // fill list of points that define the entity
      for (size_t point = 0, e_point = vertex.size(); point != e_point; ++point) 
      {
        Vec<DIM, DataType> temp;
        for (int d = 0; d < gdim; ++d) 
        {
          temp[d] = coord_vertices_refined[vertex[point] * gdim + d];
        }
        points.push_back(temp);
      }

      if (verbatim_mode > 0) 
      {
        std::cout << "defining points: " << std::endl;
        for (size_t p = 0; p < points.size(); ++p) 
        {
          for (size_t d = 0; d < points[0].size(); ++d) 
          {
            std::cout << "\t" << points[p][d];
          }
          std::cout << "\t --- ";
        }
        std::cout << std::endl;
      }

      // filter points that are on subentity
      if (verbatim_mode > 0) 
      {
        std::cout << "Filtering points on subentity (dim = " << tdim << ", idx = " << idx_entity << ")\n";
      }

      Vec<DIM, DataType> dof_coord = pt;

      // filter DoF
      if (static_cast< bool >(mesh::is_point_on_subentity<DataType, DIM> (dof_coord, points))) 
      {
          //coord_on_subentity_[tdim][idx_entity].push_back(dof_coord);
          dof_on_subentity[tdim].push_back(idx_entity);

          if (verbatim_mode > 0) 
          {
            std::cout << "-> ";
            for (size_t d = 0; d < dof_coord.size(); ++d) 
            {
              std::cout << dof_coord[d] << "\t";
            }
            std::cout << std::endl;
          }
        } // if (is_point_on_subentity(dof_coord, points) == true)
      else 
      {
          if (verbatim_mode > 0) 
          {
            std::cout << " dof point is not on subentity.\n";
          }
      }
      if (verbatim_mode > 0) 
      {
        std::cout << "\n\n";
      }
    }   // for (int idx_entity=0; idx_entity<ref_cell_->num_entities(tdim); ++idx_entity)
  } // for (int tdim=0; tdim<ref_cell_->tdim(); ++tdim)
  if (verbatim_mode > 0) 
  {
    std::cout << "\n\n";
  }

  // print summary
  if (verbatim_mode > 0) 
  {
    std::cout << "========" << std::endl;
    for (int tdim = 0; tdim < ref_cell->tdim()+1; ++tdim) 
    {
        std::cout << dof_on_subentity[tdim].size() << " subentities:" << std::endl;
        std::cout << string_from_range(
                         dof_on_subentity[tdim].begin(),
                         dof_on_subentity[tdim].end())
                  << "\n\n";
      for (size_t l = 0; l < dof_on_subentity[tdim].size(); ++l) 
      {
        size_t idx_entity = dof_on_subentity[tdim][l];
        std::cout << "TDim " << tdim << ", " << idx_entity << " ";
        if (idx_entity >= ref_cell->num_regular_entities(tdim)) 
        {
          std::cout << "(refined entity) ";
        } 
        else 
        {
          std::cout << "(regular entity) ";
        }
      }
      std::cout << "========" << std::endl;
    }
    std::cout << std::endl;
  }
}

template<class DataType, int DIM>
bool map_ref_coord_to_other_cell ( const Vec<DIM, DataType> & my_ref_coord,
                                   Vec<DIM, DataType> & other_ref_coord,
                                   doffem::CellTransformation<DataType, DIM> const * my_trans,
                                   doffem::CellTransformation<DataType, DIM> const * other_trans,
                                   std::vector< mesh::MasterSlave > const &period )
{
  const DataType COEFFICIENT_EPS = 1.e3 * std::numeric_limits< DataType >::epsilon();

  const int gdim = DIM;
  
  // Compute coordinates of my cell's DoFs (physical cell).
  Vec<DIM, DataType> phys_coord;
  my_trans->transform(my_ref_coord, phys_coord);

  // Map coordinates back to other reference cell 
  bool found = other_trans->inverse(phys_coord, other_ref_coord);
  if (found) 
  {
    found = other_trans->contains_reference_point(other_ref_coord);
  }
/*
  LOG_DEBUG(2, " ======================================== ");
  LOG_DEBUG(2, "My Cell = " << my_cell.index() << ", other Cell = " << other_cell.index());
  LOG_DEBUG(2, "my reference coords = " << string_from_range(my_ref_coord.begin(), my_ref_coord.end()));
  LOG_DEBUG(2, "physical coords = " << string_from_range(phys_coord.begin(), phys_coord.end()));
  LOG_DEBUG(2, "other reference coords = " << string_from_range(other_ref_coord.begin(), other_ref_coord.end()));
  LOG_DEBUG(2, "Coord trafo successful? : " << found );
*/
  // exception handling when other reference point could not be computed
  if (!found) 
  {
    LOG_DEBUG(2, " - - - - - - - - - - - - - - - - - - - - ");

    // There should be no problems in computing the above inverse if A and B do coincide
    //Sassert(my_cell.index() != other_cell.index());

    // try inverse map again, now take into account periodic boundaries
    // 1. check whether point lies on master or slave  boundary
    std::vector< std::vector< mesh::PeriodicBoundaryType > > per_type =
          mesh::get_periodicity_type<DataType>(phys_coord, DIM, period);

    assert(per_type.size() == 1);
    assert(per_type[0].size() == period.size());

    // TODO_ avoid conversion
    std::vector<DataType> phys_coord_tmp (DIM, 0.);
    for (int d=0; d<DIM; ++d)
    {
      phys_coord_tmp[d] = phys_coord[d];
    }

    // 2. loop over all periods and swap master and slave boundaries until
    // inverse is found
    for (int k = 0; k < period.size(); ++k) 
    {
      if (per_type[0][k] == mesh::NO_PERIODIC_BDY) 
      {
        // point does not lie on periodic boundary
        continue;
      }
      if (per_type[0][k] == mesh::MASTER_PERIODIC_BDY) 
      {
        // dof lies on master periodic boundary -> compute corresponding
        // coordinates on slave boundary

        std::vector<DataType> phys_coord_slave_tmp = mesh::periodify_master_to_slave<DataType>(phys_coord_tmp, gdim, period, k);
        Vec<DIM, DataType> phys_coord_slave ( phys_coord_slave_tmp );

        LOG_DEBUG(2, "Point mapped to slave boundary of period " << k << " : "
                     << string_from_range(phys_coord_slave_tmp.begin(), phys_coord_slave_tmp.end()));

        found = other_trans->inverse(phys_coord_slave, other_ref_coord);
        if (found) 
        {
          found = other_trans->contains_reference_point(other_ref_coord);
        }
      }
      if (per_type[0][k] == mesh::SLAVE_PERIODIC_BDY) 
      {
        // dof lies on slave periodic boundary -> compute corresponding
        // coordinates on master boundary

        std::vector<DataType> phys_coord_master_tmp = mesh::periodify_slave_to_master<DataType>(phys_coord_tmp, gdim, period, k);
        Vec<DIM, DataType> phys_coord_master (phys_coord_master_tmp);

        LOG_DEBUG(2, "Point mapped to master boundary of period " << k << " : "
                     << string_from_range(phys_coord_master_tmp.begin(), phys_coord_master_tmp.end()));

        found = other_trans->inverse(phys_coord_master, other_ref_coord);
        if (found) 
        {
          found = other_trans->contains_reference_point(other_ref_coord);
        }
      }
      if (found) 
      {
        LOG_DEBUG(2, "Found cell after periodifying w.r.t. period " << k);
        LOG_DEBUG(2, "other reference coords = " << string_from_range(other_ref_coord.begin(), other_ref_coord.end()));
        break;
      }
    }
  }
#ifndef NDEBUG
  if (!found) 
  {
    LOG_DEBUG(2, "other cell: did not find inverse to physical point "
                   << string_from_range(phys_coord.begin(), phys_coord.end()));
    LOG_DEBUG(2, " my cell ");
    std::vector< Vec<DIM, DataType> > my_coord_vtx = my_trans->get_coordinates();
    for (int v=0; v<my_coord_vtx.size(); ++v)
    {
      LOG_DEBUG(2, "vertex " << v << " : " << my_coord_vtx[v][0] << ", " << my_coord_vtx[v][1] << ", " << my_coord_vtx[v][2]);
                 
    }
    LOG_DEBUG ( 2," other cell " );
    std::vector< Vec<DIM, DataType> > other_coord_vtx = other_trans->get_coordinates();
    for (int v=0; v<other_coord_vtx.size(); ++v)
    {
      LOG_DEBUG(2, "vertex " << v << " : " << other_coord_vtx[v][0] << ", " << other_coord_vtx[v][1] << ", " << other_coord_vtx[v][2]);
    }
  }
#endif
  return found;
}

template < class DataType, int DIM >
void create_bbox_for_mesh (ConstMeshPtr meshptr, BBox<DataType, DIM>& bbox)
{
  bbox.reset(DIM);
  int tdim = meshptr->tdim();
  
  // check if mesh has a periodic boundary
  std::vector< MasterSlave > period = meshptr->get_period();

  if (period.size() == 0) 
  {
    // no periodic boundary -> simply loop over all vertices in mesh
    for (int vert = 0; vert < meshptr->num_entities(0); ++vert) 
    {
      std::vector<double> tmp_coord = meshptr->get_coordinates(0, vert);
      assert (tmp_coord.size() == DIM);
      Vec<DIM, DataType> pt;
      for (size_t d=0; d<DIM; ++d)
      {
        pt[d] = static_cast<DataType>(tmp_coord[d]);
      }
      bbox.add_point(pt);
    }
  } 
  else 
  {
    // periodic boundary is present -> need to take into account slave boundary
    // points which are not present as vertices in the mesh loop over all cells
    for (mesh::EntityIterator it = meshptr->begin(tdim),
                              end_it = meshptr->end(tdim);
                              it != end_it; ++it) 
    {
      // get coordinates of current cell and unperiodify them
      std::vector< DataType > periodic_coords_on_cell;
      it->get_coordinates(periodic_coords_on_cell);

      std::vector< DataType > unperiodic_coords_on_cell =
          unperiodify(periodic_coords_on_cell, DIM, period);
      
      int num_vert = unperiodic_coords_on_cell.size() / DIM;
      std::vector<Vec<DIM, DataType> > coords_on_cell (num_vert);
      for (int i=0; i<num_vert; ++i)
      { 
        for (int d=0; d<DIM; ++d)
        {
          coords_on_cell[i][d] = unperiodic_coords_on_cell[i*DIM +d];
        } 
      }
      bbox.add_points(coords_on_cell);
    }
  }
  bbox.uniform_extension(10 * GEOM_TOL);
}

template < class DataType, int DIM >
std::vector< DataType > compute_mean_edge_length (ConstMeshPtr meshptr)
{
  const int num_edges = meshptr->num_entities(1);
  const int sqrt_num_edges = std::sqrt(num_edges);
  std::vector< DataType > mean_edge_length(DIM, 0.);
  int count = 0;
  for (int index = 0; index < num_edges; index += sqrt_num_edges, ++count) 
  {
    std::vector< double > coords = meshptr->get_coordinates(1, index);
    assert(static_cast< int >(coords.size()) == 2 * DIM);
    for (int d = 0; d < DIM; ++d) 
    {
      mean_edge_length[d] += std::abs(static_cast<DataType>(coords[d] - coords[d + DIM]));
    }
  }
  for (int d = 0; d < DIM; ++d) 
  {
    mean_edge_length[d] /= static_cast< DataType >(count);
  }
  return mean_edge_length;
}

template < class DataType, int DIM >
void find_adjacent_cells (ConstMeshPtr in_mesh, 
                          ConstMeshPtr out_mesh, 
                          std::map< int, std::set<int> >& cell_map)
{
  assert (in_mesh != nullptr);
  assert (out_mesh != nullptr);
  
  cell_map.clear();
  
  const int tdim = in_mesh->tdim();
  
  // construct common grid search for both meshes
  BBox<DataType, DIM> in_bbox(DIM);
  BBox<DataType, DIM> out_bbox(DIM);
  
  create_bbox_for_mesh (in_mesh, in_bbox);
  create_bbox_for_mesh (out_mesh, out_bbox);
  
  std::vector<DataType> extends(2*DIM, 0.);
  for (size_t d=0; d<DIM; ++d)
  {
    extends[2*d] = std::min(in_bbox.min(d), out_bbox.min(d));
    extends[2*d+1] = std::max(in_bbox.max(d), out_bbox.max(d));
  }
  BBox<DataType, DIM> bbox(extends);
    
  std::vector< DataType > in_mean_edge = compute_mean_edge_length<DataType, DIM> (in_mesh);
  std::vector< DataType > out_mean_edge = compute_mean_edge_length<DataType, DIM> (out_mesh);
  std::vector< DataType > mean_edge (DIM, 0.);
  
  for (int d = 0; d < DIM; ++d) 
  {
    if (in_mean_edge[d] <= 10. * GEOM_TOL) 
    {
      in_mean_edge[d] = in_bbox.max(d) - in_bbox.min(d);
    }
    if (out_mean_edge[d] <= 10. * GEOM_TOL) 
    {
      out_mean_edge[d] = out_bbox.max(d) - out_bbox.min(d);
    }
    mean_edge[d] = std::min(in_mean_edge[d], out_mean_edge[d]);
  }

  std::vector<int> num_intervals (DIM, 0);
  for (int d = 0; d < DIM; ++d) 
  {
    DataType num = std::max(static_cast<DataType>(1.), (bbox.max(d) - bbox.min(d)) / mean_edge[d]);
    num_intervals[d] = static_cast< int >(num);
  }
  
  // create rectangular grid that covers both meshes
  GridGeometricSearch<DataType,DIM> in_search  (in_mesh, bbox, num_intervals);
  GridGeometricSearch<DataType,DIM> out_search (out_mesh, bbox, num_intervals);
                      
  // map: grid index to cell index
  const std::vector< std::list< int > >& in_grid_2_cells = in_search.get_cell_map();
  const std::vector< std::list< int > >& out_grid_2_cells = out_search.get_cell_map();

  // map: cell index to grid index
  const std::vector< std::list< int > >& in_cell_2_grids = in_search.get_inverse_cell_map();
  const std::vector< std::list< int > >& out_cell_2_grids = out_search.get_inverse_cell_map();
  
  std::map< int, std::set<int> > checked_cells;
 
  // loop over all cells of in_mesh
  for (EntityIterator in_cell = in_mesh->begin(tdim); in_cell != in_mesh->end(tdim); ++in_cell) 
  { 
    const int in_cell_index = in_cell->index();
    std::set<int> adj_cells;
    cell_map[in_cell_index] = adj_cells;
    
    std::set<int> cur_checked_cells;
    checked_cells[in_cell_index] = cur_checked_cells;
    
    std::vector<double> tmp_in_vertex_coords = in_mesh->get_coordinates(DIM, in_cell_index); 
    std::vector<DataType> in_vertex_coords(tmp_in_vertex_coords.size());
    for (size_t l=0; l<tmp_in_vertex_coords.size(); ++l)
    {
      in_vertex_coords[l] = static_cast<DataType>(tmp_in_vertex_coords[l]);
    }
        
    std::vector< Vec<DIM, DataType> > in_points;
    interlaced_coord_to_points<DataType, DataType, DIM> (in_vertex_coords, in_points);
           
    const size_t in_num_vert = in_points.size();
    
    // get possible neighbors belonging to out_mesh by using rectangular grid
    // loop over all grid cells that contain current in_cell
    for (auto grid_it = in_cell_2_grids[in_cell->index()].begin();  
         grid_it != in_cell_2_grids[in_cell->index()].end(); ++grid_it) 
    {
      // loop over all out_cells, that are contained in current grid_cell
      for (auto out_cell_index = out_grid_2_cells[*grid_it].begin(); 
           out_cell_index != out_grid_2_cells[*grid_it].end(); ++out_cell_index)
      {
        if (checked_cells[in_cell_index].find(*out_cell_index) != checked_cells[in_cell_index].end())
        {
          continue; 
        }
        
        // check for intersection of in_cell and out_cell
        std::vector<double> tmp_out_vertex_coords = out_mesh->get_coordinates(DIM, *out_cell_index); 
        std::vector<DataType> out_vertex_coords(tmp_out_vertex_coords.size());
        for (size_t l=0; l<tmp_out_vertex_coords.size(); ++l)
        {
          out_vertex_coords[l] = static_cast<DataType>(tmp_out_vertex_coords[l]);
        }
        
        //TODO: maybe adjust the refinement parameter or make it a parameter in this function
        bool intersection = cells_intersect<DataType, DIM>(in_vertex_coords, out_vertex_coords, 10);   
          
        if (intersection)
        {
          cell_map[in_cell_index].insert(*out_cell_index);
        }
        checked_cells[in_cell_index].insert(*out_cell_index);
      }
    }
  }
}

typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double>, false > Polygon;    //false for counterclockwise orientation

template <class DataType, int DIM>
bool cells_intersect (const std::vector < DataType>& in_vertex_coords, 
                      const std::vector < DataType>& out_vertex_coords, 
                      int edge_refinement)
{ 
  int tdim = DIM;
  GDim gdim = DIM;
  bool intersection = false;
  
  if (tdim == 1) {    //two lines intersect each other iff a vertex of one line lies within the other line
    
    std::vector< Vec<DIM, DataType> > in_points;
    interlaced_coord_to_points<DataType, DataType, DIM> (in_vertex_coords, in_points);
    
    Vec<DIM, DataType> tmp_point;
    // check if vertices of in_cell are contained in out_cell
    for (size_t i=0; i<2; ++i)
    {
      bool found = point_inside_cell<DataType, DIM>(in_points[i], out_vertex_coords, tmp_point);
      if (found)
      {
        intersection = true;
        break;
      }
    }
  }
  else if (tdim == 2) 
  {
    // rectangular quads -> check bounding box
    if( is_aligned_rectangular_quad(in_vertex_coords) 
      && is_aligned_rectangular_quad(out_vertex_coords))
    {
      BBox<DataType, DIM> in_box(in_vertex_coords, 4);
      BBox<DataType, DIM> out_box(out_vertex_coords, 4);
      
      return in_box.intersects(out_box, GEOM_TOL);
    }
    
    // TODO: check for other specific cases in order to avoid expensive boost polygon routines
    
    
    // general case
  
    std::string in_vertices_wkt("(");
    std::string out_vertices_wkt("(");
    
    std::string test_str = std::to_string(in_vertex_coords[0]);
    int test_size = test_str.size();
    
    in_vertices_wkt.reserve(in_vertex_coords.size() * (1+test_size) * 2);
    out_vertices_wkt.reserve(out_vertex_coords.size() * (1+test_size) * 2);
    
    //parsing coords in the correct format for the boost function
    for (int i = 0; i < in_vertex_coords.size(); i+=2) {
      in_vertices_wkt += std::to_string(in_vertex_coords[i]) + " " + std::to_string(in_vertex_coords[i+1]) + ","; 
    }
    
    //boost requires start and end point to coincide
    in_vertices_wkt += std::to_string(in_vertex_coords[0]) + " " + std::to_string(in_vertex_coords[1]) + ")";
    //std::cout << in_vertices_wkt << std::endl;
    
    for (int i = 0; i < out_vertex_coords.size(); i+=2) {
      out_vertices_wkt += std::to_string(out_vertex_coords[i]) + " " + std::to_string(out_vertex_coords[i+1]) + ","; 
    }
    
    //boost requires start and end point to coincide
    out_vertices_wkt += std::to_string(out_vertex_coords[0]) + " " + std::to_string(out_vertex_coords[1]) + ")";
    //std::cout << out_vertices_wkt << std::endl;
    
    Polygon in, out;
    std::deque<Polygon> intersection_polygon;
    
    boost::geometry::read_wkt(
        "POLYGON(" + in_vertices_wkt + ")", in);
    boost::geometry::read_wkt(
        "POLYGON(" + out_vertices_wkt + ")", out);
    
    intersection = boost::geometry::intersects(in, out);
    BOOST_FOREACH(Polygon const& p, intersection_polygon)
    {
        //std::cout << "area: " << boost::geometry::area(p) << std::endl;
    }
  }
  else if (tdim == 3) 
  {
    //check if the polyeders are aligned rectangular cuboids by checking if all edges are parallel to unit normal vectors
    
    //std::cout << is_aligned_rectangular_cuboid(in_vertex_coords) << "------------------" << std::endl;
    std::vector <Vec<3, DataType>> in_points;
    std::vector <Vec<3, DataType>> out_points;
    interlaced_coord_to_points<DataType, DataType, 3> (in_vertex_coords, in_points);
    interlaced_coord_to_points<DataType, DataType, 3> (out_vertex_coords, out_points);
    
    if(is_aligned_rectangular_cuboid(in_vertex_coords) 
        && is_aligned_rectangular_cuboid(out_vertex_coords) )
    {
      BBox<DataType, DIM> in_box(in_vertex_coords, 8);
      BBox<DataType, DIM> out_box(out_vertex_coords, 8);
      
      return in_box.intersects(out_box, GEOM_TOL);
    }
    else 
    {
      std::vector<Vec<3, DataType>> in_dir_vectors;
      std::vector<Vec<3, DataType>> in_sup_vectors;
      std::vector<Vec<3, DataType>> out_dir_vectors;
      std::vector<Vec<3, DataType>> out_sup_vectors;
      
      in_dir_vectors.reserve(30);
      out_dir_vectors.reserve(30);
      in_sup_vectors.reserve(30);
      out_sup_vectors.reserve(30);

      
      parametrize_object<DataType> (in_points, in_dir_vectors, in_sup_vectors);
      parametrize_object<DataType> (out_points, out_dir_vectors, out_sup_vectors);
      
      bool found = false;
      for (int i = 0; i< in_dir_vectors.size(); ++i) {
        
        for (int j = 0; j < edge_refinement; ++j) {
          
          Vec<3, DataType> in_pt = in_sup_vectors[i]
            + j / static_cast<DataType>(edge_refinement) * in_dir_vectors[i];

          //std::cout << "Direction vector " << i << std::endl;
          //std::cout << "Part " << j << std::endl;
          //std::cout << pt << std::endl;        
          Vec<3, DataType> in_ref_point;

          assert(  (point_inside_cell<DataType, 3> (in_pt, in_vertex_coords, in_ref_point)) );
          found = point_inside_cell<DataType, 3>(in_pt, out_vertex_coords, in_ref_point);
          
          if (found) 
          {
            intersection = true;
            break;
          }
        }
        if (intersection) break;    
      }
      
      if(!intersection) {   //only check if out_object intersects in_object if no intersection has been found yet
        for (int i = 0; i< out_dir_vectors.size(); ++i) {
          
          for (int j = 0; j < edge_refinement; ++j) {
            
        
            Vec<3, DataType> out_pt = out_sup_vectors[i]
              + j / static_cast<DataType>(edge_refinement) * out_dir_vectors[i];

            Vec<3, DataType> out_ref_point;
           
            assert(  (point_inside_cell<DataType, 3> (out_pt, out_vertex_coords, out_ref_point)) );   
            found = point_inside_cell<DataType, 3>(out_pt, in_vertex_coords, out_ref_point);
            
            if (found) 
            {
              intersection = true;
              break;
            }
          }
          if (intersection) break;  
        }
      }
    }
  }
  return intersection;
}

template <class DataType>
bool is_aligned_rectangular_cuboid (const std::vector < DataType>& vertex_coords) 
{
  //        7 ----------------- 6
  //       /|                /|
  //      / |               / |
  //     /  |z             /  |
  //    4 ----------------- 5 |
  //    |   |             |   |
  //    |   |       y     |   |
  //    |   3 ------------|-- 2
  //    |  /              |  /
  //    | /x              | /
  //    |/                |/
  //    0 ----------------- 1
  
  bool is_aligned_rectangular_cuboid;
  
  if (vertex_coords.size() != 24)
    return false;

  std::vector <Vec<3, DataType>> points;
  interlaced_coord_to_points<DataType, DataType, 3> (vertex_coords, points);
  
  std::vector < Vec<3, DataType > > x_edges(4);
  std::vector < Vec<3, DataType > > y_edges(4);
  std::vector < Vec<3, DataType > > z_edges(4);
  
  x_edges[0] = (points[0] - points[3]);
  x_edges[1] = (points[1] - points[2]);
  x_edges[2] = (points[4] - points[7]);
  x_edges[3] = (points[5] - points[6]);
  
  y_edges[0] = (points[1] - points[0]);
  y_edges[1] = (points[2] - points[3]);
  y_edges[2] = (points[5] - points[4]);
  y_edges[3] = (points[7] - points[6]);
  
  z_edges[0] = (points[4] - points[0]);
  z_edges[1] = (points[5] - points[1]);
  z_edges[2] = (points[7] - points[3]);
  z_edges[3] = (points[6] - points[2]);
  
  //x- and y-edges could be swapped if cuboid is rotated by 90 degrees, so check orientation first
  if (std::abs(x_edges[0][1]) < GEOM_TOL) { //if true, these edges can only be aligned to x-axis
    for (int i = 0; i < 4; i++) {   //you could check fewer edges whether they are aligned to their respective axis
    
      if (std::abs(x_edges[i][1]) > GEOM_TOL || std::abs(x_edges[i][2]) > GEOM_TOL ) return false;
      if (std::abs(y_edges[i][0]) > GEOM_TOL || std::abs(y_edges[i][2]) > GEOM_TOL ) return false;
      if (std::abs(z_edges[i][0]) > GEOM_TOL || std::abs(z_edges[i][1]) > GEOM_TOL ) return false;
    }
    return true;
  }
  else if (std::abs(x_edges[0][0]) < GEOM_TOL) { //if true, these edges can only be aligned to y-axis
    for (int i = 0; i < 4; i++) {   
      if (std::abs(x_edges[i][0]) > GEOM_TOL || std::abs(x_edges[i][2]) > GEOM_TOL ) return false;
      if (std::abs(y_edges[i][1]) > GEOM_TOL || std::abs(y_edges[i][2]) > GEOM_TOL ) return false;
      if (std::abs(z_edges[i][0]) > GEOM_TOL || std::abs(z_edges[i][1]) > GEOM_TOL ) return false;
    }
    return true;
  }
  else return false; //aligned to neither of the axes -> not aligned at all
}

template <class DataType>
bool is_aligned_rectangular_quad(const std::vector < DataType>& vertex_coords) 
{
  //                y    
  //       3 ----------------- 2
  //      /                  /
  //     / x                /
  //    /                  /
  //   0 -----------------1
  
  if (vertex_coords.size()!= 8)
  {
    return false;
  }
  std::vector <Vec<2, DataType>> points;
  interlaced_coord_to_points<DataType, DataType, 2> (vertex_coords, points);
  
  std::vector < Vec<2, DataType > > a_edges(2);
  std::vector < Vec<2, DataType > > b_edges(2);

  
  a_edges[0] = (points[0] - points[3]);
  a_edges[1] = (points[1] - points[2]);
  
  b_edges[0] = (points[1] - points[0]);
  b_edges[1] = (points[2] - points[3]);

  return (   (std::abs(a_edges[0][0]) + std::abs(a_edges[1][0]) < 2.0 * GEOM_TOL)
          && (std::abs(b_edges[0][1]) + std::abs(b_edges[1][1]) < 2.0 * GEOM_TOL) )
      || (   (std::abs(a_edges[0][1]) + std::abs(a_edges[1][1]) < 2.0 * GEOM_TOL)
          && (std::abs(b_edges[0][0]) + std::abs(b_edges[1][0]) < 2.0 * GEOM_TOL) );
}

template <class DataType>
bool is_parallelogram(const std::vector<DataType>& vertex_coords)
{
  //                y
  //       3 ----------------- 2
  //      /                  /
  //     / x                /
  //    /                  /
  //   0 -----------------1

  if (vertex_coords.size() != 8)
    return false;

  std::vector<Vec<2, DataType>> points;
  interlaced_coord_to_points<DataType, DataType, 2>(vertex_coords, points);

  Vec<2, DataType> a1 = points[1] - points[0];
  Vec<2, DataType> a2 = points[2] - points[3];

  Vec<2, DataType> b1 = points[2] - points[1];
  Vec<2, DataType> b2 = points[3] - points[0];

  if (!is_parallel(a1, a2))
  {
    return false;
  }
  return is_parallel(b1, b2);
}

template <class DataType>
bool is_parallelepiped(const std::vector<DataType>& vertex_coords)
{
  //        7 ----------------- 6
  //       /|                /|
  //      / |               / |
  //     /  |z             /  |
  //    4 ----------------- 5 |
  //    |   |             |   |
  //    |   |       y     |   |
  //    |   3 ------------|---- 2
  //    |  /              |  /
  //    | /x              | /
  //    |/                |/
  //    0 ----------------- 1

  if (vertex_coords.size() != 24)
    return false;

  std::vector<Vec<3, DataType>> points;
  interlaced_coord_to_points<DataType, DataType, 3>(vertex_coords, points);

  Vec<3, DataType> edges[3][4];
  edges[0][0] = points[2] - points[1];
  edges[0][1] = points[3] - points[0];
  edges[0][2] = points[4] - points[7];
  edges[0][3] = points[5] - points[6];

  edges[1][0] = points[1] - points[0];
  edges[1][1] = points[3] - points[2];
  edges[1][2] = points[5] - points[4];
  edges[1][3] = points[7] - points[6];

  edges[2][0] = points[4] - points[0];
  edges[2][1] = points[5] - points[1];
  edges[2][2] = points[6] - points[2];
  edges[2][3] = points[7] - points[3];

  for (int i = 0; i < 3; ++i) {
    for (int j = 1; j < 4; ++j) {
      if (!is_parallel(edges[i][j], edges[i][j - 1])) {
        return false;
      }
    }
  }

  return true;
}

template <class DataType >
void parametrize_object (const std::vector<Vec<3, DataType>>& in_points,
                         std::vector <Vec<3, DataType>> &dir_vectors, 
                         std::vector < Vec<3, DataType> >&sup_vectors)
{
  assert((in_points.size() == 4 ) || (in_points.size() == 5) || (in_points.size() == 8 ) );
  //check if first polyhedron is a tetrahedron, pyramid or hexahedron
  
  if (in_points.size() == 4) {    //tetrahedron
    
    //check if all edges intersect the out polyhedron and additionally if the lines connecting the midpoints of the edges with the vertices intersect the polhedron -> 18 lines
    
    
    for (int i = 0; i < 6; i++) {
      sup_vectors.push_back(in_points[0]);
    }
    //direction vectors connecting vertex 0 to the other relevant points
    dir_vectors.push_back(in_points[1]- in_points[0]);
    dir_vectors.push_back(in_points[2]- in_points[0]);
    dir_vectors.push_back(in_points[3]- in_points[0]);
    dir_vectors.push_back(in_points[1] + 0.5*(in_points[3]- in_points[1]) - in_points[0]);
    dir_vectors.push_back(in_points[1] + 0.5*(in_points[2]- in_points[1]) - in_points[0]);
    dir_vectors.push_back(in_points[3] + 0.5*(in_points[2]- in_points[3]) - in_points[0]);
    
    for (int i = 0; i< 5; ++i ) {
      sup_vectors.push_back(in_points[1]);
    }
    dir_vectors.push_back(in_points[2]- in_points[1]);
    dir_vectors.push_back(in_points[3]- in_points[1]);
    dir_vectors.push_back(in_points[0] + 0.5*(in_points[3]- in_points[0]) - in_points[1]);
    dir_vectors.push_back(in_points[0] + 0.5*(in_points[2]- in_points[0]) - in_points[1]);
    dir_vectors.push_back(in_points[3] + 0.5*(in_points[2]- in_points[3]) - in_points[1]);
    
    for (int i = 0; i< 4; ++i ) {
      sup_vectors.push_back(in_points[2]);
    }
    dir_vectors.push_back(in_points[3]- in_points[2]);
    dir_vectors.push_back(in_points[0] + 0.5*(in_points[1]- in_points[0]) - in_points[2]);
    dir_vectors.push_back(in_points[0] + 0.5*(in_points[3]- in_points[0]) - in_points[2]);
    dir_vectors.push_back(in_points[3] + 0.5*(in_points[1]- in_points[3]) - in_points[2]);
    
    for (int i = 0; i< 3; ++i ) {
      sup_vectors.push_back(in_points[3]);
    }

    dir_vectors.push_back(in_points[0] + 0.5*(in_points[1]- in_points[0]) - in_points[3]);
    dir_vectors.push_back(in_points[0] + 0.5*(in_points[2]- in_points[0]) - in_points[3]);
    dir_vectors.push_back(in_points[1] + 0.5*(in_points[2]- in_points[1]) - in_points[3]);
    
    
  }
  
  if (in_points.size() == 5) {    //pyramid
    
    for (int i = 0; i < 4; ++i) {
      sup_vectors.push_back(in_points[0]);
    }
    dir_vectors.push_back(in_points[1]- in_points[0]);
    dir_vectors.push_back(in_points[2]- in_points[0]);
    dir_vectors.push_back(in_points[3]- in_points[0]);
    dir_vectors.push_back(in_points[4]- in_points[0]);
    
    for (int i = 0; i < 3; ++i) {
      sup_vectors.push_back(in_points[1]);
    }
    dir_vectors.push_back(in_points[2]- in_points[1]);
    dir_vectors.push_back(in_points[3]- in_points[1]);
    dir_vectors.push_back(in_points[4]- in_points[1]);
    
    for (int i = 0; i < 2; ++i) {
      sup_vectors.push_back(in_points[2]);
    }

    dir_vectors.push_back(in_points[3]- in_points[2]);
    dir_vectors.push_back(in_points[4]- in_points[2]);
    
    sup_vectors.push_back(in_points[3]);
    dir_vectors.push_back(in_points[4] - in_points[3]);
    
    for (int i = 0; i < 4; ++i) {
      sup_vectors.push_back(in_points[4]);
    }
    dir_vectors.push_back(in_points[0] + 0.5 * (in_points[1]- in_points[0]) - in_points[4]);
    dir_vectors.push_back(in_points[1] + 0.5 * (in_points[2]- in_points[1])- in_points[4]);
    dir_vectors.push_back(in_points[2] + 0.5 * (in_points[3]- in_points[2])- in_points[4]);
    dir_vectors.push_back(in_points[3] + 0.5 * (in_points[0]- in_points[3])- in_points[4]);
    
    //lines connecting midpoints of the ground edges
    
    sup_vectors.push_back(in_points[0] + 0.5 * (in_points[1]- in_points[0]));
    sup_vectors.push_back(in_points[1] + 0.5 * (in_points[2]- in_points[1]));
    
    dir_vectors.push_back(in_points[2] + 0.5 * (in_points[3]- in_points[2]) - (in_points[0] + 0.5 * (in_points[1]- in_points[0] ) ) );
    dir_vectors.push_back(in_points[3] + 0.5 * (in_points[0]- in_points[3])- (in_points[1] + 0.5 * (in_points[2]- in_points[1]) )  );
    
    
    
    
  }
  
  if (in_points.size() == 8) {    //hexahedron
    
    for (int i = 0; i < 7; ++i) {
      sup_vectors.push_back(in_points[0]);
    }
    dir_vectors.push_back(in_points[1]- in_points[0]);
    dir_vectors.push_back(in_points[2]- in_points[0]);
    dir_vectors.push_back(in_points[3]- in_points[0]);
    dir_vectors.push_back(in_points[4]- in_points[0]);
    dir_vectors.push_back(in_points[5]- in_points[0]);
    dir_vectors.push_back(in_points[6]- in_points[0]);
    dir_vectors.push_back(in_points[7]- in_points[0]);
    
    for (int i = 0; i < 6; ++i) {
      sup_vectors.push_back(in_points[1]);
    }
    dir_vectors.push_back(in_points[2]- in_points[1]);
    dir_vectors.push_back(in_points[3]- in_points[1]);
    dir_vectors.push_back(in_points[4]- in_points[1]);
    dir_vectors.push_back(in_points[5]- in_points[1]);
    dir_vectors.push_back(in_points[6]- in_points[1]);
    dir_vectors.push_back(in_points[7]- in_points[1]);
    
    for (int i = 0; i < 5; ++i) {
      sup_vectors.push_back(in_points[2]);
    }

    dir_vectors.push_back(in_points[3]- in_points[2]);
    dir_vectors.push_back(in_points[4]- in_points[2]);
    dir_vectors.push_back(in_points[5]- in_points[2]);
    dir_vectors.push_back(in_points[6]- in_points[2]);
    dir_vectors.push_back(in_points[7]- in_points[2]);
    
    for (int i = 0; i < 4; ++i) {
      sup_vectors.push_back(in_points[3]);
    }

    dir_vectors.push_back(in_points[4]- in_points[3]);
    dir_vectors.push_back(in_points[5]- in_points[3]);
    dir_vectors.push_back(in_points[6]- in_points[3]);
    dir_vectors.push_back(in_points[7]- in_points[3]);
    
    for (int i = 0; i < 3; ++i) {
      sup_vectors.push_back(in_points[4]);
    }

    dir_vectors.push_back(in_points[5]- in_points[4]);
    dir_vectors.push_back(in_points[6]- in_points[4]);
    dir_vectors.push_back(in_points[7]- in_points[4]);
    
    for (int i = 0; i < 2; ++i) {
      sup_vectors.push_back(in_points[5]);
    }

    dir_vectors.push_back(in_points[6]- in_points[5]);
    dir_vectors.push_back(in_points[7]- in_points[5]);
    
    
    sup_vectors.push_back(in_points[6]);
    
    dir_vectors.push_back(in_points[7] - in_points[6]);
    
  }
  
  
}


template Vec<3, float> normal<float, 3>(const std::vector< Vec<3, float> > &);
template Vec<2, float> normal<float, 2>(const std::vector< Vec<2, float> > &);
template Vec<1, float> normal<float, 1>(const std::vector< Vec<1, float> > &);
template Vec<3, double> normal<double, 3>(const std::vector< Vec<3, double> > &);
template Vec<2, double> normal<double, 2>(const std::vector< Vec<2, double> > &);
template Vec<1, double> normal<double, 1>(const std::vector< Vec<1, double> > &);

template float distance_point_hyperplane<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &);
template float distance_point_hyperplane<float, 2>(const Vec<2, float> &point, const Vec<2, float> &, const Vec<2, float> &);
template float distance_point_hyperplane<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &);
template double distance_point_hyperplane<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &);
template double distance_point_hyperplane<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &);
template double distance_point_hyperplane<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &);

template float distance_point_line<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &);
template float distance_point_line<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const Vec<2, float> &);
template float distance_point_line<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &);
template double distance_point_line<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &);
template double distance_point_line<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &);
template double distance_point_line<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &);

template Vec<3, float> foot_point_hyperplane<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &);
template Vec<2, float> foot_point_hyperplane<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const Vec<2, float> &);
template Vec<1, float> foot_point_hyperplane<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &);
template Vec<3, double> foot_point_hyperplane<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &);
template Vec<2, double> foot_point_hyperplane<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &);
template Vec<1, double> foot_point_hyperplane<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &);

template Vec<3, float> foot_point_line<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &);
template Vec<2, float> foot_point_line<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const Vec<2, float> &);
template Vec<1, float> foot_point_line<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &);
template Vec<3, double> foot_point_line<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &);
template Vec<2, double> foot_point_line<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &);
template Vec<1, double> foot_point_line<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &);

template float triangle_area<float, 3>(const std::vector< float > &);
template float triangle_area<float, 2>(const std::vector< float > &);
template float triangle_area<float, 1>(const std::vector< float > &);
template double triangle_area<double, 3>(const std::vector< double > &);
template double triangle_area<double, 2>(const std::vector< double > &);
template double triangle_area<double, 1>(const std::vector< double > &);

template float facet_area<float, 3>(const std::vector< float > &, const GDim );
template float facet_area<float, 2>(const std::vector< float > &, const GDim );
template float facet_area<float, 1>(const std::vector< float > &, const GDim );
template double facet_area<double, 3>(const std::vector< double > &, const GDim );
template double facet_area<double, 2>(const std::vector< double > &, const GDim );
template double facet_area<double, 1>(const std::vector< double > &, const GDim );

template bool in_plane<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &, const float ) ;
template bool in_plane<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const Vec<2, float> &, const float ) ;
template bool in_plane<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &, const float ) ;
template bool in_plane<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &, const double ) ;
template bool in_plane<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &, const double ) ;
template bool in_plane<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &, const double ) ;

template bool crossed_plane<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &, const Vec<3, float> &);
template bool crossed_plane<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const Vec<2, float> &, const Vec<2, float> &);
template bool crossed_plane<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &, const Vec<1, float> &);
template bool crossed_plane<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &, const Vec<3, double> &);
template bool crossed_plane<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &, const Vec<2, double> &);
template bool crossed_plane<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &, const Vec<1, double> &);

template bool crossed_facet<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const std::vector< float > &);
template bool crossed_facet<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const std::vector< float > &);
template bool crossed_facet<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const std::vector< float > &);
template bool crossed_facet<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const std::vector< double > &);
template bool crossed_facet<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const std::vector< double > &);
template bool crossed_facet<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const std::vector< double > &);

template Vec<3, float> intersect_facet<float, 3>(const Vec<3, float> &, const Vec<3, float> &, const std::vector< float > &, bool &);
template Vec<2, float> intersect_facet<float, 2>(const Vec<2, float> &, const Vec<2, float> &, const std::vector< float > &, bool &);
template Vec<1, float> intersect_facet<float, 1>(const Vec<1, float> &, const Vec<1, float> &, const std::vector< float > &, bool &);
template Vec<3, double> intersect_facet<double, 3>(const Vec<3, double> &, const Vec<3, double> &, const std::vector< double > &, bool &);
template Vec<2, double> intersect_facet<double, 2>(const Vec<2, double> &, const Vec<2, double> &, const std::vector< double > &, bool &);
template Vec<1, double> intersect_facet<double, 1>(const Vec<1, double> &, const Vec<1, double> &, const std::vector< double > &, bool &);

template float distance_point_facet<float, 3>(const Vec<3, float> &, const std::vector< float > &, Vec<3, float> &);
template float distance_point_facet<float, 2>(const Vec<2, float> &, const std::vector< float > &, Vec<2, float> &);
template float distance_point_facet<float, 1>(const Vec<1, float> &, const std::vector< float > &, Vec<1, float> &);
template double distance_point_facet<double, 3>(const Vec<3, double> &, const std::vector< double > &, Vec<3, double> &);
template double distance_point_facet<double, 2>(const Vec<2, double> &, const std::vector< double > &, Vec<2, double> &);
template double distance_point_facet<double, 1>(const Vec<1, double> &, const std::vector< double > &, Vec<1, double> &);

template bool point_inside_entity<float, 3>(const Vec<3, float> &, const TDim , const std::vector< float > &);
template bool point_inside_entity<float, 2>(const Vec<2, float> &, const TDim , const std::vector< float > &);
template bool point_inside_entity<float, 1>(const Vec<1, float> &, const TDim , const std::vector< float > &);
template bool point_inside_entity<double, 3>(const Vec<3, double> &, const TDim , const std::vector< double > &);
template bool point_inside_entity<double, 2>(const Vec<2, double> &, const TDim , const std::vector< double > &);
template bool point_inside_entity<double, 1>(const Vec<1, double> &, const TDim , const std::vector< double > &);

template bool point_inside_cell<float, 3>(const Vec<3, float> &, const std::vector< float > &, Vec<3, float> &);
template bool point_inside_cell<float, 2>(const Vec<2, float> &, const std::vector< float > &, Vec<2, float> &);
template bool point_inside_cell<float, 1>(const Vec<1, float> &, const std::vector< float > &, Vec<1, float> &);
template bool point_inside_cell<double, 3>(const Vec<3, double> &, const std::vector< double > &, Vec<3, double> &);
template bool point_inside_cell<double, 2>(const Vec<2, double> &, const std::vector< double > &, Vec<2, double> &);
template bool point_inside_cell<double, 1>(const Vec<1, double> &, const std::vector< double > &, Vec<1, double> &);

template bool vertices_inside_one_hyperplane<float, 3>(const std::vector< float > &, const GDim );
template bool vertices_inside_one_hyperplane<float, 2>(const std::vector< float > &, const GDim );
template bool vertices_inside_one_hyperplane<float, 1>(const std::vector< float > &, const GDim );
template bool vertices_inside_one_hyperplane<double, 3>(const std::vector< double > &, const GDim );
template bool vertices_inside_one_hyperplane<double, 2>(const std::vector< double > &, const GDim );
template bool vertices_inside_one_hyperplane<double, 1>(const std::vector< double > &, const GDim );

template Vec<3, float> project_point<float, 3>(const BoundaryDomainDescriptor<float, 3> &, const Vec<3, float> &, const MaterialNumber );
template Vec<2, float> project_point<float, 2>(const BoundaryDomainDescriptor<float, 2> &, const Vec<2, float> &, const MaterialNumber );
template Vec<1, float> project_point<float, 1>(const BoundaryDomainDescriptor<float, 1> &, const Vec<1, float> &, const MaterialNumber );
template Vec<3, double> project_point<double, 3>(const BoundaryDomainDescriptor<double, 3> &, const Vec<3, double> &, const MaterialNumber );
template Vec<2, double> project_point<double, 2>(const BoundaryDomainDescriptor<double, 2> &, const Vec<2, double> &, const MaterialNumber );
template Vec<1, double> project_point<double, 1>(const BoundaryDomainDescriptor<double, 1> &, const Vec<1, double> &, const MaterialNumber );

template bool map_ref_coord_to_other_cell<float, 3> ( const Vec<3, float> & ,
                                                      Vec<3, float> & ,
                                                      doffem::CellTransformation<float, 3> const * ,
                                                      doffem::CellTransformation<float, 3> const * ,
                                                      std::vector< mesh::MasterSlave > const & );
template bool map_ref_coord_to_other_cell<float, 2> ( const Vec<2, float> &,
                                                      Vec<2, float> &,
                                                      doffem::CellTransformation<float, 2> const *,
                                                      doffem::CellTransformation<float, 2> const *,
                                                      std::vector< mesh::MasterSlave > const & );
template bool map_ref_coord_to_other_cell<float, 1> ( const Vec<1, float> &,
                                                      Vec<1, float> &,
                                                      doffem::CellTransformation<float, 1> const *,
                                                      doffem::CellTransformation<float, 1> const *,
                                                      std::vector< mesh::MasterSlave > const & );
template bool map_ref_coord_to_other_cell<double, 3> ( const Vec<3, double> &,
                                                      Vec<3, double> &,
                                                      doffem::CellTransformation<double, 3> const *,
                                                      doffem::CellTransformation<double, 3> const *,
                                                      std::vector< mesh::MasterSlave > const & );
template bool map_ref_coord_to_other_cell<double, 2> ( const Vec<2, double> &,
                                                      Vec<2, double> &,
                                                      doffem::CellTransformation<double, 2> const *,
                                                      doffem::CellTransformation<double, 2> const *,
                                                      std::vector< mesh::MasterSlave > const & );
template bool map_ref_coord_to_other_cell<double, 1> ( const Vec<1, double> &,
                                                      Vec<1, double> &,
                                                      doffem::CellTransformation<double, 1> const *,
                                                      doffem::CellTransformation<double, 1> const *,
                                                      std::vector< mesh::MasterSlave > const & );

template void find_subentities_containing_point<float, 3> (const Vec<3, float> &, 
                                                           const mesh::CellType *,
                                                           const std::vector< Vec<3, float> >& , 
                                                           std::vector< std::vector<int> > &);
template void find_subentities_containing_point<float, 2> (const Vec<2, float>&, 
                                                           const mesh::CellType *,
                                                           const std::vector< Vec<2, float> >& , 
                                                           std::vector< std::vector<int> > &);
template void find_subentities_containing_point<float, 1> (const Vec<1, float> &, 
                                                           const mesh::CellType *,
                                                           const std::vector< Vec<1, float> >& , 
                                                           std::vector< std::vector<int> > &);
template void find_subentities_containing_point<double, 3> (const Vec<3, double> &, 
                                                            const mesh::CellType *,
                                                            const std::vector< Vec<3, double> >& , 
                                                            std::vector< std::vector<int> > &);
template void find_subentities_containing_point<double, 2> (const Vec<2, double> &, 
                                                            const mesh::CellType *,
                                                            const std::vector< Vec<2, double> >& , 
                                                            std::vector< std::vector<int> > &);
template void find_subentities_containing_point<double, 1> (const Vec<1, double>&, 
                                                            const mesh::CellType *,
                                                            const std::vector< Vec<1, double> >& , 
                                                            std::vector< std::vector<int> > &);


template bool is_point_on_subentity<float, 3> (const Vec<3, float>& , const std::vector< Vec<3, float> > &);
template bool is_point_on_subentity<float, 2> (const Vec<2, float>& , const std::vector< Vec<2, float> > &);
template bool is_point_on_subentity<float, 1> (const Vec<1, float>& , const std::vector< Vec<1, float> > &);
template bool is_point_on_subentity<double, 3> (const Vec<3, double>& , const std::vector< Vec<3, double> > &);
template bool is_point_on_subentity<double, 2> (const Vec<2, double>& , const std::vector< Vec<2, double> > &);
template bool is_point_on_subentity<double, 1> (const Vec<1, double>& , const std::vector< Vec<1, double> > &);

template void create_bbox_for_mesh<float, 3>  (ConstMeshPtr meshptr, BBox<float, 3>& bbox);
template void create_bbox_for_mesh<float, 2>  (ConstMeshPtr meshptr, BBox<float, 2>& bbox);
template void create_bbox_for_mesh<float, 1>  (ConstMeshPtr meshptr, BBox<float, 1>& bbox);
template void create_bbox_for_mesh<double, 3> (ConstMeshPtr meshptr, BBox<double, 3>& bbox);
template void create_bbox_for_mesh<double, 2> (ConstMeshPtr meshptr, BBox<double, 2>& bbox);
template void create_bbox_for_mesh<double, 1> (ConstMeshPtr meshptr, BBox<double, 1>& bbox);

template std::vector< float >  compute_mean_edge_length<float, 3>  (ConstMeshPtr meshptr);
template std::vector< float >  compute_mean_edge_length<float, 2>  (ConstMeshPtr meshptr);
template std::vector< float >  compute_mean_edge_length<float, 1>  (ConstMeshPtr meshptr);
template std::vector< double > compute_mean_edge_length<double, 3> (ConstMeshPtr meshptr);
template std::vector< double > compute_mean_edge_length<double, 2> (ConstMeshPtr meshptr);
template std::vector< double > compute_mean_edge_length<double, 1> (ConstMeshPtr meshptr);

template void find_adjacent_cells <float, 3> (ConstMeshPtr in_mesh, ConstMeshPtr out_mesh, std::map< int, std::set<int> >& adjacent_map);
template void find_adjacent_cells <float, 2> (ConstMeshPtr in_mesh, ConstMeshPtr out_mesh, std::map< int, std::set<int> >& adjacent_map);
template void find_adjacent_cells <float, 1> (ConstMeshPtr in_mesh, ConstMeshPtr out_mesh, std::map< int, std::set<int> >& adjacent_map);
template void find_adjacent_cells <double, 3> (ConstMeshPtr in_mesh, ConstMeshPtr out_mesh, std::map< int, std::set<int> >& adjacent_map);
template void find_adjacent_cells <double, 2> (ConstMeshPtr in_mesh, ConstMeshPtr out_mesh, std::map< int, std::set<int> >& adjacent_map);
template void find_adjacent_cells <double, 1> (ConstMeshPtr in_mesh, ConstMeshPtr out_mesh, std::map< int, std::set<int> >& adjacent_map);

template bool cells_intersect<double, 1>(const std::vector < double>& in_vertex_coords,const  std::vector < double>& out_vertex_coords, int edge_refinement);
template bool cells_intersect<double, 2>(const std::vector < double>& in_vertex_coords,const  std::vector < double>& out_vertex_coords, int edge_refinement);
template bool cells_intersect<double, 3>(const std::vector < double>& in_vertex_coords,const  std::vector < double>& out_vertex_coords, int edge_refinement);
template bool cells_intersect<float, 1> (const std::vector < float>& in_vertex_coords, const  std::vector < float>& out_vertex_coords, int edge_refinement);
template bool cells_intersect<float, 2> (const std::vector < float>& in_vertex_coords, const  std::vector < float>& out_vertex_coords, int edge_refinement);
template bool cells_intersect<float, 3> (const std::vector < float>& in_vertex_coords, const  std::vector < float>& out_vertex_coords, int edge_refinement);

template bool is_aligned_rectangular_cuboid<float>  (const std::vector < float>& in_vertex_coords);
template bool is_aligned_rectangular_cuboid<double> (const std::vector < double>& in_vertex_coords);

template bool is_aligned_rectangular_quad<float>  (const std::vector < float>& in_vertex_coords);
template bool is_aligned_rectangular_quad<double> (const std::vector < double>& in_vertex_coords);

template bool is_parallelogram<float> (const std::vector<float>& in_vertex_coords);
template bool is_parallelogram<double>(const std::vector<double>& in_vertex_coords);

template bool is_parallelepiped<float> (const std::vector<float>& in_vertex_coords);
template bool is_parallelepiped<double>(const std::vector<double>& in_vertex_coords);

template void parametrize_object<double>(const std::vector<Vec<3, double>>& in_points, std::vector <Vec<3, double >> &dir_vectors, std::vector < Vec<3,double> >&sup_vectors);
template void parametrize_object<float> (const std::vector<Vec<3, float>>& in_points, std::vector <Vec<3, float >> &dir_vectors, std::vector < Vec<3,float> >&sup_vectors);

} // namespace mesh
} // namespace hiflow
