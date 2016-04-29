# Large box
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

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

# Equilibrate:
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(500)

# Production for 100 nanoseconds
simulation.reporters.append(app.DCDReporter('trajectory.dcd', 50000))
simulation.reporters.append(StateDataReporter('simdata.txt', 50000, step=True,potentialEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('snapshots.pdb', 10000000))
simulation.step(50000000)





