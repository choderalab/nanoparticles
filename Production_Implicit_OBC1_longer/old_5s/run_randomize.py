# Large box
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import math

prmtop = AmberPrmtopFile('out.prmtop')
inpcrd = AmberInpcrdFile('out.inpcrd')

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

# Randomise and equilibrate
# Simulated annealing
Ti = 1000.
Tf = 300.
N = 100
simulation.reporters.append(PDBReporter('annealing.pdb', 500))
for t in range(N+1):
    temp = Ti*math.exp(-((1.0*t)/N)*math.log(Ti/Tf))*kelvin
    simulation.context.setVelocitiesToTemperature(temp)
    simulation.step(500)
simulation.reporters.pop(0)

# Production for 300 nanoseconds
simulation.reporters.append(DCDReporter('trajectory.dcd', 150000))
simulation.reporters.append(StateDataReporter('simdata.txt', 150000, step=True,potentialEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('snapshots.pdb', 10000000))
simulation.step(150000000)
