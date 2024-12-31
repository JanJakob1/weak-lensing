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

#include "tools/point_evaluator.h"

namespace hiflow {

template class PointEvaluator <double, 2, double>;

#if 0
template <>
bool PointEvaluator< double, double >::mpi_allreduce (double* send_data,
                                                      double* recv_data,
                                                      int count,
                                                      MPI_Op op,
                                                      MPI_Comm communicator ) const
{
  MPI_Allreduce(send_data, recv_data, count, MPI_DOUBLE, op, communicator);
}

template <>
bool PointEvaluator< double, Vec<2, double> >::mpi_allreduce(Vec<2, double>* send_data,
                                                             Vec<2, double>* recv_data,
                                                             int count,
                                                             MPI_Op op,
                                                             MPI_Comm communicator ) const
{
  std::vector<double> send_data_array;
  std::vector<double> recv_data_array(2*count);
  
  send_data_array.reserve(2*count);
  
  for (size_t l=0; l<count; ++l)
  {
    for (size_t d=0; d<2; ++d)
    {
      send_data_array.push_back(send_data[l][d]);
    }
  }

  int status = MPI_Allreduce(&send_data_array[0], &recv_data_array[0], count * 2, MPI_DOUBLE, op, communicator);
  
  for (size_t l=0; l<count; ++l)
  {
    for (size_t d = 0; d < 2; ++d)
    {
      recv_data[l][0] = recv_data_array[l];
    }
  }
}

template class PointEvaluator <double, double>;
template class PointEvaluator <double, Vec<2,double> >;
#endif

} // namespace hiflow
