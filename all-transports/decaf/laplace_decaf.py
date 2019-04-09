# a small 2-node example, just a producer and consumer
# modified from the linear-2nodes-topo1.py

# --- include the following 4 lines each time ---

import networkx as nx
import os
import imp
import sys
import argparse
wf = imp.load_source('workflow', os.environ['DECAF_PREFIX'] + '/python/decaf.py')
# --- set your options here ---

# path to .so module for dataflow callback functions
mod_path = os.environ['BUILD_DIR'] + '/lib/libmod_laplace_decaf_shared.so'

# define workflow graph
# 2-node workflow
#
#    prod (4 procs) -> con (2 procs)
#
#  entire workflow takes 8 procs (2 dataflow procs between producer and consumer)
#  dataflow can be overlapped, but currently all disjoint procs (simplest case)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
wf.initParserForTopology(parser)
args = parser.parse_args()

procs_prod = int(os.environ['procs_prod'])
procs_con = int(os.environ['procs_con'])
procs_link = int(os.environ['procs_link'])

#Dan comment
# Creating the topology
#topo = wf.topologyFromArgs(args)
#topoProd = topo.subTopology("prod", procs_prod, 0)
#topoDflow = topo.subTopology("dflow", procs_link, procs_prod)
#topoCon = topo.subTopology("con", procs_con, procs_prod+procs_link)

# Creating the graph
#w = nx.DiGraph()
#w.add_node("prod", topology=topoProd, func='prod', cmdline=os.environ['BUILD_DIR'] + '/bin/lammps_decaf')
#w.add_node("con", topology=topoCon, func='con', cmdline=os.environ['BUILD_DIR'] + '/bin/lammps_decaf')
#w.add_edge("prod", "con", topology=topoDflow, func='dflow', path=mod_path,
#          prod_dflow_redist='count', dflow_con_redist='count',cmdline=os.environ['BUILD_DIR'] +'/bin/lammps_decaf')

# --- Graph definition ---
# --- Graph definition ---
laplace   = wf.Node("prod", start_proc=0, nprocs=procs_prod, func='prod', cmdline=os.environ['BUILD_DIR'] + '/bin/laplace_decaf')
outPort = laplace.addOutputPort("out")

#print1   = wf.Node("dflow", start_proc=procs_prod, nprocs=procs_link, func='print', cmdline='./lammps')
#inPort1  = print1.addInputPort("in")

con   = wf.Node("con", start_proc= procs_prod+procs_link, nprocs=procs_con, func='con', cmdline=os.environ['BUILD_DIR'] + '/bin/laplace_decaf') 
inPort  = con.addInputPort("in")


link    = wf.Edge(laplace.getOutputPort("out"), con.getInputPort("in"), start_proc=procs_prod, nprocs=procs_link, func='dflow',  path=mod_path, prod_dflow_redist='count', dflow_con_redist='count', cmdline=os.environ['BUILD_DIR'] +'/bin/laplace_decaf')

#lammps_detect = wf.Edge(lammps.getOutputPort("out"), detect.getInputPort("in"),
 #       start_proc=0,  nprocs=0, func='', prod_dflow_redist='count', dflow_con_redist='count')

#link2    = wf.Edge(lammps.getOutputPort("out"), print2.getInputPort("in"), start_proc=6, nprocs=1, func='dflow',
#           path=mod_path, prod_dflow_redist='count', dflow_con_redist='count', cmdline='./lammps')

#link3    = wf.Edge(print2.getOutputPort("out"), print3.getInputPort("in"), start_proc=8, nprocs=1, func='dflow',
#           path=mod_path, prod_dflow_redist='count', dflow_con_redist='count', cmdline='./lammps')


# --- convert the nx graph into a workflow data structure and run the workflow ---
#wf.processGraph(w, "lammps_decaf")
print args.custom_args
wf.processGraph("laplace_decaf", args.custom_args)

#wf.processGraph("lammps_decaf", args.custom_args)
