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

#include "fem/fe_mapping.h"
#include "fem/fe_reference.h"

namespace hiflow {
namespace doffem {

template < class DataType, int DIM >
std::vector<DataType> MappedRefFunctionEval<DataType, DIM>::evaluate (const mesh::Entity& entity, 
                                                                      const std::vector<Vec<DIM, DataType> >& pts) const
{
  assert (this->function_ != nullptr);
  assert (this->cell_trafo_ != nullptr);
  assert (this->fe_trafo_ != nullptr);
  const size_t nb_comp = this->nb_comp();
  const size_t nb_func = this->nb_func();
    
  std::vector<DataType> mapped_vals (nb_comp * nb_func * pts.size(), 0.);
  size_t offset = 0;
  for (size_t l = 0; l < pts.size(); ++l)
  {
    Coord ref_pt;
    if (this->cell_trafo_->inverse(pts[l], ref_pt))
    {
      // if physical point lies not on cell associated to cell_trafo ->
      // nothing to do since ansatz functions are extended by 0
    
      // evaluate reference shape function at ref_pt
      std::vector<DataType> shape_vals_pt (this->function_->weight_size(), 0.);
      std::vector<DataType> mapped_vals_pt (this->function_->weight_size(), 0.);
      this->function_->N(ref_pt, shape_vals_pt);

      // map shape function values to element living on physical cell
      this->fe_trafo_->map_shape_function_values (*this->cell_trafo_, ref_pt,
                                                  0, nb_func, nb_comp,
                                                  shape_vals_pt, mapped_vals_pt);
      // insert mapped_vals_pt into mapped_vals
      // loop over shape functions
      //for (size_t i = 0; i < nb_func; ++i)
      for (size_t c = 0; c < nb_comp; ++c)
      {
        // loop over components
        //for (size_t c = 0; c < nb_comp; ++c)
        for (size_t i = 0; i < nb_func; ++i)
        {
          //mapped_vals[offset + i * nb_comp + c] = mapped_vals_pt[this->ref_fe_->iv2ind(i, c)];
          mapped_vals[offset + this->iv2ind(i,c)/*c * nb_func + i*/] = mapped_vals_pt[this->function_->iv2ind(i, c)];
        }
      }
    }
    offset += nb_func * nb_comp;
  }
  return mapped_vals;
}

template class MappedRefFunctionEval<float,1>;
template class MappedRefFunctionEval<float,2>;
template class MappedRefFunctionEval<float,3>;
template class MappedRefFunctionEval<double,1>;
template class MappedRefFunctionEval<double,2>;
template class MappedRefFunctionEval<double,3>;
}
}
