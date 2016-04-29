# Large box
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import math

prmtop = AmberPrmtopFile('restart_out.prmtop')
inpcrd = AmberInpcrdFile('restart_out.inpcrd')

bxsize = quantity.Quantity(((16.1,0.0,0.0),(0.0,16.1,0.0),(0.0,0.0,16.1)),unit=nanometer)
prmtop.topology.setPeriodicBoxVectors(bxsize)

system = prmtop.createSystem(nonbondedMethod=CutoffPeriodic,nonbondedCutoff=1*nanometer,constraints=HBonds,implicitSolvent=OBC1)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Set GPU processor:
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

# Very quick minimisation:
simulation.minimizeEnergy(maxIterations=10)

# Equilibrate
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(500)

# Production for 117.6 nanoseconds
simulation.reporters.append(DCDReporter('restarted_trajectory.dcd', 150000))
simulation.reporters.append(StateDataReporter('restarted_simdata.txt', 150000, step=True,potentialEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('restarted_snapshots.pdb', 10000000))
simulation.step(58800000)
