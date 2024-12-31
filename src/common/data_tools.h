// Copyright (C) 2011-2020 Vincent Heuveline
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

#ifndef HIFLOW_DATA_TOOLS_H
#define HIFLOW_DATA_TOOLS_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip> 
#include <vector>
#include <math.h> 
#include "common/vector_algebra.h"

/// @brief T

/// @author Philipp Gerstner

namespace {

} // namespace

namespace hiflow {

template<class DataType>
void number_range (DataType first, DataType inc, size_t length, std::vector<DataType>& result)
{
  result.resize(length, static_cast<DataType>(0));
  for (size_t l=0; l<length; ++l)
  {
    result[l] = first + static_cast<DataType>(l) * inc;
  }
}

template<class DataType>
void log_2d_array(const std::vector< std::vector< DataType> >& vals,
                  std::ostream &s, 
                  int precision = 6) 
{
  s << std::setprecision(precision);
  for (size_t l=0; l<vals.size(); ++l)
  {
    for (size_t k=0; k<vals[l].size(); ++k)
    {
      s << vals[l][k] << " ";
    }
    s << "\n";
  }
}  
  
template <typename T>
void set_to_value(const size_t size, const T& val, std::vector<T> & array)
{
  const size_t old_size = std::min(size, array.size());
  array.resize(size, val);
  
  for (size_t i=0; i<old_size; ++i)
  {
    array[i] = val; 
  }
}

template <typename T>
void set_to_zero(const size_t size, std::vector<T> & array)
{
  const size_t old_size = std::min(size, array.size());
  array.resize(size, static_cast<T>(0));
  
  for (size_t i=0; i<old_size; ++i)
  {
    array[i] = static_cast<T>(0); 
  }
}

template <class DataType>
void double_2_datatype (const std::vector<double>& in_array, 
                        std::vector<DataType>& out_array)
{
  const size_t size = in_array.size();
  out_array.resize(size);
  
  for (size_t i=0; i<size; ++i)
  {
    out_array[i] = static_cast<DataType>(in_array[i]);
  }
}

template <class DataType, class T, int DIM>
void interlaced_coord_to_points (const std::vector<DataType>& inter_coords,   
                                 std::vector<Vec<DIM, T> >& points)
{
  const size_t num_vert = static_cast<size_t> (inter_coords.size() / DIM);
  points.resize(num_vert);
  
  for (size_t i=0; i<num_vert; ++i)
  {
    for (size_t d=0; d<DIM; ++d)
    {
      points[i][d] = static_cast<T>(inter_coords[i * DIM + d]);
    }
  }
}

template <typename T, int DIM>
void flatten_2d_vec_array(const std::vector< std::vector< Vec<DIM,T> > >& input, 
                          std::vector<T> & output, 
                          std::vector<size_t>& sizes) 
{
  if (input.size() == 0)
  {
    output.clear();
    return;
  }

  sizes.clear();
  size_t total_size = 0;
  
  for (size_t d=0; d<input.size(); ++d)
  {
    sizes.push_back(input[d].size());
    total_size += input[d].size();
  }
  
  output.clear();
  output.resize(total_size*DIM);
  
  size_t offset = 0;
  for (size_t d=0; d<input.size(); ++d)
  {
    for (size_t i=0; i<sizes[d]; ++i)
    {
      for (size_t l=0; l<DIM; ++l)
      {
        output[offset+l] = input[d][i][l];
      }
      offset += DIM;
    }
  }
}

template <typename T>
void flatten_2d_array(const std::vector< std::vector<T> >& input, 
                      std::vector<T> & output, 
                      std::vector<size_t>& sizes) 
{
  if (input.size() == 0)
  {
    output.clear();
    return;
  }

  sizes.clear();
  size_t total_size = 0;
  
  for (size_t d=0; d<input.size(); ++d)
  {
    sizes.push_back(input[d].size());
    total_size += input[d].size();
  }
  
  output.clear();
  output.resize(total_size);
  
  size_t offset = 0;
  for (size_t d=0; d<input.size(); ++d)
  {
    for (size_t i=0; i<sizes[d]; ++i)
    {
      output[offset+i] = input[d][i];
    }
    offset += sizes[d];
  }
}

template <typename T>
void expand_to_2d_array(const std::vector<T> & input,
                        const std::vector<size_t>& sizes,
                        std::vector< std::vector<T> >& output)
{
  if (input.size() == 0)
  {
    output.clear();
    return;
  }
  
  const size_t n = sizes.size();
  if (output.size() != n)
  {
    output.resize(n);
  }
  
  size_t offset = 0;
  for (size_t d=0; d<n; ++d)
  {
    const size_t m = sizes[d];
    if (output[d].size() != m)
    {
      output[d].resize(m);
    }
    
    for (size_t i=0; i<m; ++i)
    {
      output[d][i] = input[offset+i];
    }
    offset += m;
  }
}

template <typename T, int DIM>
void expand_to_2d_vec_array(const std::vector<T> & input,
                            const std::vector<size_t>& sizes,
                            std::vector< std::vector<Vec<DIM,T> > >& output)
{
  if (input.size() == 0)
  {
    output.clear();
    return;
  }
  
  const size_t n = sizes.size();
  if (output.size() != n)
  {
    output.resize(n);
  }
  
  size_t offset = 0;
  for (size_t d=0; d<n; ++d)
  {
    const size_t m = sizes[d];
    if (output[d].size() != m)
    {
      output[d].resize(m);
    }
    
    for (size_t i=0; i<m; ++i)
    {
      for (size_t l=0; l<DIM; ++l)
      {
        output[d][i][l] = input[offset+l];
      }
      offset += DIM; 
    }
  }
}

template <typename T>
bool vectors_are_equal(const std::vector<T>& arg1,
                       const std::vector<T>& arg2,
                       T eps)
{
  if (arg1.size() != arg2.size())
  {
    return false;
  }
  if (arg1.size() == 0)
  {
    return true;
  }
  
  const size_t len = arg1.size();
  auto it1 = arg1.begin();
  auto it2 = arg2.begin();
   
  for (size_t i=0; i!=len; ++i)
  {
    if (std::abs(*it1 - *it2) > eps)
    {
      return false;
    }
    ++it1;
    ++it2;
  }
  return true;
}

template <typename T, int DIM>
bool vectors_are_equal(const std::vector<Vec<DIM,T> >& arg1,
                       const std::vector<Vec<DIM,T> >& arg2,
                       T eps)
{
  if (arg1.size() != arg2.size())
  {
    return false;
  }
  if (arg1.size() == 0)
  {
    return true;
  }
  
  const size_t len = arg1.size();
  auto it1 = arg1.begin();
  auto it2 = arg2.begin();
   
  for (size_t i=0; i!=len; ++i)
  {
    if (it1->neq(*it2, eps))
    {
      return false;
    }
    ++it1;
    ++it2;
  }
  return true;
}

template <typename T>
bool contains_nan (const std::vector<T>& arg)
{
  for (size_t i=0, e=arg.size(); i!=e; ++i)
  {
    if (std::isnan(arg[i]))
    {
      return true;
    }
  }
  return false;
}

} // namespace hiflow

#endif /* _DATA_TOOLS_H_ */
