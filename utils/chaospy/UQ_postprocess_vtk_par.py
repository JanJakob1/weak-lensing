#!/usr/bin/env python3

console_info = True

###-------------------------------------------------------------------------###
### Import python modules
###-------------------------------------------------------------------------###

import os
import sys
import chaospy as cp

###-------------------------------------------------------------------------###
### Import UQ scripts
###-------------------------------------------------------------------------###

sys.path.append(os.path.abspath("scripts"))
import integration as integ
import qoi_postprocessing as qoipp
import PC_collocation as pcc
import vtk_postprocessing as vtkpp
import IO_vtk as iovtk
import parallel as par

###-------------------------------------------------------------------------###
### Input and parameters 
###-------------------------------------------------------------------------###

### Get template file names

### Get template file names

#work_folder = '/home/philipp/RemoteDir/bwforDev_home/TEHD/Projekt/AK5_7K_0kV_3D_stat_meso_finite_dT_uniform_1'
#pp_folder = '/home/philipp/RemoteDir/bwforDev_home/TEHD/Projekt/AK5_7K_0kV_3D_stat_meso_finite_dT_uniform_1'
#node_2_pvtu_path = 'pvtu/'
#node_2_vtu_path = 'primal_vtu/'
#pvtu_2_vtu_path = "../" + node_2_vtu_path
#vtk_output_name = 'TEHD_primal.0.0000.pvtu'
#qoi_names = [ "F_phi", "F_rad", "F_z", "T", "P", "Phi", "vel_phi", "vel_rad", "vel_x", "vel_y", "vel_z", "vort_phi", "vort_rad", "vort_z"]
#PC_degree = 3
#dist_list = [cp.Uniform(6.5, 7.5)]
#quadrature_rule = 'g'
#sparse_flag = False

# mounted work directory
#work_folder = '/home/philipp/Programme/Hiflow/hiflow_dev/build_gcc_relwithdebinfo/examples/poisson_chaospy'
work_folder = '/home/kratzke/gitlab/hiflow_build/examples/poisson_chaospy'

# location where to write the postprocessed output
pp_folder = '/home/philipp/Programme/Hiflow/hiflow_dev/build_gcc_relwithdebinfo/examples/poisson_chaospy'

# path to pvtk and vtk within node folder
# use './' if there is no additional folder
node_2_pvtu_path = './'
node_2_vtu_path = './'
pvtu_2_vtu_path = ""

# pvtu filename
vtk_output_name = 'solution4.pvtu'

###
qoi_names = ["u"]

### Define UQ parameters
PC_degree = 2
quadrature_rule = 'g'
sparse_flag = False
dist_list = [cp.Uniform(0.1, 0.2), cp.Uniform(0.7, 0.8)]

basis, basis_norms, nodes, weights, dist, num_nodes = pcc.create_PC_collocation (PC_degree, quadrature_rule, sparse_flag, dist_list)


###-------------------------------------------------------------------------###
### Create surrogate model from node results 
###-------------------------------------------------------------------------###

pvtu_dummy = os.path.join( work_folder, "node_0", node_2_pvtu_path, vtk_output_name )

node_2_vtu_filenames, pvtu_2_vtu_filenames = iovtk.read_pvtu_structure( pvtu_dummy, node_2_pvtu_path )

print (node_2_vtu_filenames)
print (pvtu_2_vtu_filenames)

num_tasks = len( node_2_vtu_filenames )

loc_chunk = par.local_chunk( num_tasks )

if node_2_vtu_path != "./":
  os.system("mkdir " + os.path.join( work_folder, node_2_vtu_path ))

if node_2_pvtu_path != "./":
  os.system("mkdir " + os.path.join( work_folder, node_2_pvtu_path ))
  
for rank_task in loc_chunk:

  print ( "Task:\n", rank_task )

  node_2_vtu_filename = node_2_vtu_filenames[rank_task]

  vtk_surr, vtk_offsets = \
     vtkpp.create_surrogate_from_vtk( work_folder, node_2_vtu_filename, qoi_names, 
                                      basis, basis_norms, nodes, weights )

###-------------------------------------------------------------------------###
### Calculate statistics with chaospy 
###-------------------------------------------------------------------------###

### VTK postprocessing
### mean
  mean_vtk = integ.mean( vtk_surr, dist )
  
### Standard deviations
  std_dev_vtk = integ.stddev( vtk_surr, dist )

### Main Sobol indices
  main1_vtk, main2_vtk = integ.sens_m( vtk_surr, dist )

### Total Sobol indices
  total1_vtk, total2_vtk = integ.sens_t( vtk_surr, dist )

### Percentiles
  perc_vtk = qoipp.compute_percentiles( vtk_surr, dist, [10,90] )

###-------------------------------------------------------------------------###
### Write results 
###-------------------------------------------------------------------------###
    
  dummy_file = os.path.join( work_folder, "node_0", node_2_vtu_filename )
  mean_file = os.path.join( work_folder, node_2_vtu_path, "mean_" + str(rank_task) + ".vtu" )
  stddev_file = os.path.join( work_folder, node_2_vtu_path, "std_dev_" + str(rank_task) + ".vtu" )
  main1_file = os.path.join( work_folder, node_2_vtu_path, "main1_" + str(rank_task) + ".vtu" )
  main2_file = os.path.join( work_folder, node_2_vtu_path, "main2_" + str(rank_task) + ".vtu" )
  total1_file = os.path.join( work_folder, node_2_vtu_path, "total1_" + str(rank_task) + ".vtu" )
  total2_file = os.path.join( work_folder, node_2_vtu_path, "total2_" + str(rank_task) + ".vtu" )
  p10_file = os.path.join( work_folder, node_2_vtu_path, "perc10_" + str(rank_task) + ".vtu" )
  p90_file = os.path.join( work_folder, node_2_vtu_path, "perc90_" + str(rank_task) + ".vtu" )
  iovtk.write_vtu (dummy_file, mean_file, mean_vtk, vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, stddev_file, std_dev_vtk, vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, main1_file, main1_vtk, vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, main2_file, main2_vtk, vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, total1_file, total1_vtk, vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, total2_file, total2_vtk, vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, p10_file, perc_vtk[0], vtk_offsets, qoi_names)
  iovtk.write_vtu (dummy_file, p90_file, perc_vtk[1], vtk_offsets, qoi_names)


if ( par.rank() == 0 ):
  print ("Rank 0 writes pvtu.")
  mean_file_list = [pvtu_2_vtu_path + "mean_" + str( p ) + ".vtu" for p in range( num_tasks )]
  stddev_file_list = [pvtu_2_vtu_path + "std_dev_" + str( p ) + ".vtu" for p in range( num_tasks )]
  main1_file_list = [pvtu_2_vtu_path + "main1_" + str( p ) + ".vtu" for p in range( num_tasks )]
  main2_file_list = [pvtu_2_vtu_path + "main2_" + str( p ) + ".vtu" for p in range( num_tasks )]
  total1_file_list = [pvtu_2_vtu_path + "total1_" + str( p ) + ".vtu" for p in range( num_tasks )]
  total2_file_list = [pvtu_2_vtu_path + "total2_" + str( p ) + ".vtu" for p in range( num_tasks )]
  p10_file_list = [pvtu_2_vtu_path + "perc10_" + str( p ) + ".vtu" for p in range( num_tasks )]
  p90_file_list = [pvtu_2_vtu_path + "perc90_" + str( p ) + ".vtu" for p in range( num_tasks )]
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "mean.pvtu" ),
                    pvtu_2_vtu_filenames, mean_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "std_dev.pvtu" ),
                    pvtu_2_vtu_filenames, stddev_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "main1.pvtu" ),
                    pvtu_2_vtu_filenames, main1_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "main2.pvtu" ),
                    pvtu_2_vtu_filenames, main2_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "total1.pvtu" ),
                    pvtu_2_vtu_filenames, total1_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "total2.pvtu" ),
                    pvtu_2_vtu_filenames, total2_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "perc10.pvtu" ),
                    pvtu_2_vtu_filenames, p10_file_list )
  iovtk.adapt_pvtu( pvtu_dummy, os.path.join( work_folder, node_2_pvtu_path, "perc90.pvtu" ),
                    pvtu_2_vtu_filenames, p90_file_list )

