from simtk.openmm import app, unit
from simtk import openmm
import subprocess
import mdtraj as md

import argparse
parser = argparse.ArgumentParser(description="Save the last frame of a given trajectory")
parser.add_argument('-d','--dir', type=str,help='the directory containing snapshots.pbd and trajectory.dcd.')
parser.add_argument('-n','--name', type=str,help='the name of the directory where the last frame will be output.')
args = parser.parse_args()

pdbname = 'nanoparticle.pdb'
traj = md.load(args.dir+'trajectory.dcd',top=args.dir+'snapshots.pdb')
traj[-1::].save_pdb(args.name+'/'+pdbname)

cmd = "obabel -ipdb {0} -opdb -O {1}".format(args.name+'/'+pdbname,args.name+'/'+pdbname)
subprocess.call(cmd,shell=True)
