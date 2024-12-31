console_info = False

###-------------------------------------------------------------------------###
### Import python modules
###-------------------------------------------------------------------------###

###-------------------------------------------------------------------------###
### Import UQ scripts
###-------------------------------------------------------------------------###

import IO_vtk
import qoi_postprocessing as qoipp

def create_surrogate_from_vtk ( work_folder, vtk_filename, qoi_names,
                                basis, norms, nodes, weights ):
  """
  Create surrogate model from vtk node files
  Combines read_vtk_nodes and create_surrogate_from_qoi
  input:
    work_folder: Folder with collocation node folders
    vtk_filename: vtk file within a respective node folder
    qoi_names: names of the QOIs in the vtk files
    basis: CP basis
    norms: CP basis norms
    nodes: Quadrature nodes
    weights: Quadrature weights
  return:
    surr: surrogate model
    offsets: offsets for indexing original data arrays
  """
  ### read node qoi data
  nn = len( weights )
  nodes_data, offsets = \
    IO_vtk.read_vtk_nodes ( work_folder, vtk_filename, qoi_names, nn )

  ### create cp surroagte
  surr = qoipp.create_surrogate_from_qoi ( basis, norms, nodes, weights, nodes_data )

  return surr, offsets

