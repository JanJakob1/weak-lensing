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

import solver_input as si
import PC_collocation as pcc

###-------------------------------------------------------------------------###
### Input and parameters 
###-------------------------------------------------------------------------###

# work directory
#work_folder = '/home/philipp/Programme/Hiflow/hiflow_dev/build_gcc_relwithdebinfo/examples/poisson_chaospy'
work_folder = '/home/kratzke/gitlab/hiflow_build/examples/poisson_chaospy'

# script that should be executed for each node
per_node_submit_script = 'run_uq.sh'

# parameters for per_node_submit_script
per_node_submit_args = work_folder

# mounted template directory 
template_folder = os.path.join(work_folder,'template')

num_threads = 2;

### Define UQ parameters

PC_degree = 2
quadrature_rule = 'g'
sparse_flag = False
dist_list = [cp.Uniform(0.1, 0.2), cp.Uniform(0.7, 0.8)]

basis, basis_norms, nodes, weights, dist, num_nodes = pcc.create_PC_collocation (PC_degree, quadrature_rule, sparse_flag, dist_list)

###-------------------------------------------------------------------------###
### Adapt application input files 
###-------------------------------------------------------------------------###

### copy template dir

state = si.create_node_folders( template_folder, work_folder, nodes )

if console_info:
  print ("Created", num_nodes, "node input folders:", state)

###-------------------------------------------------------------------------###
### Test submission script creation 
###-------------------------------------------------------------------------###

si.create_job_submission_script( os.path.join( work_folder, "submission.sh" ), work_folder, num_nodes, per_node_submit_script, per_node_submit_args)

###-------------------------------------------------------------------------###
### Run application 
###-------------------------------------------------------------------------###

si.run_jobs_locally (work_folder, per_node_submit_script, per_node_submit_args, num_threads, num_nodes)

