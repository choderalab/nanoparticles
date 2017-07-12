from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import time
from openmoltools.forcefield_generators import gaffTemplateGenerator

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run an openmm simulation with salt exchange moves.")
    parser.add_argument('-f','--file',type=str,help="the pdb file containing the nanoparticle,default=Thermalized_bigger.pdb",default="Thermalized_bigger.pdb")
    parser.add_argument('-l','--length',type=float, help="the length of the simulation in nanometers, default = 100",default = 100)
    parser.add_argument("--gpu",action='store_true',help="whether the simulation will be run on a GPU, default=False",default=False)
    args = parser.parse_args()

# Initialsing the forcefield so that it can generate antechamber parameters
forcefield = ForceField('gaff.xml','tip3p.xml','amber99sbildn.xml')
forcefield.registerTemplateGenerator(gaffTemplateGenerator)

# Loading the starting structure of the simulations
pdb = PDBFile(args.file)

# Creating everything necessary for openmm simulation:
system = forcefield.createSystem(pdb.topology,nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=None)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
system.addForce(MonteCarloBarostat(1*atmospheres, 300*kelvin, 25))

if args.gpu:
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
else:
    simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Minimisation:
#t = time.time()
#simulation.minimizeEnergy()
#print 'System minimised in {0} seconds.'.format(time.time()-t)


# Equilibrate. Gradually heating:
t = time.time()
temperature = range(1,300,9)
#simulation.context.setVelocitiesToTemperature(300*kelvin)
#simulation.step(500)

simulation.context.setVelocitiesToTemperature(1*kelvin)
for temp in temperature:
    integrator.setTemperature(temp*kelvin)
    simulation.step(100)
print 'System thermalised in {0} seconds.'.format(time.time()-t)

positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
PDBFile.writeFile(simulation.topology,positions, open('equil.pdb', 'w'))

# Running simulation
nsteps = int(1000*args.length/0.002)   # Converting desired length of simulation to number of steps
nsnapshots = int(nsteps/1000.0)
ndata  = int(nsteps/100.0)
npdb = int(nsteps/5.0)

t = time.time()
simulation.reporters.append(app.DCDReporter('trajectory.dcd', nsnapshots))
simulation.reporters.append(StateDataReporter('simdata.txt', ndata, step=True,potentialEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('snapshots.pdb',npdb ))
simulation.step(nsteps)
print 'Simulation finished in {0} seconds.'.format(time.time()-t)
