from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import time
from openmoltools.forcefield_generators import gaffTemplateGenerator

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run an openmm simulation with salt exchange moves.")
    parser.add_argument('-f','--file',type=str,help="the pdb file containing the nanoparticle,default=nanoparticle.pdb",default="nanoparticle.pdb")
    parser.add_argument("--gpu",action='store_true',help="whether the simulation will be run on a GPU, default=False",default=False)
    args = parser.parse_args()

# Initialsing the forcefield so that it can generate antechamber parameters
forcefield = ForceField('gaff.xml','tip3p.xml','amber99sbildn.xml')
forcefield.registerTemplateGenerator(gaffTemplateGenerator)

# Loading the starting structure of the simulations
pdb = PDBFile(args.file)

# Solvating, neutralising and adding salt:
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield,padding=1*nanometer,neutralize=True,ionicStrength=0.0001*molar)
app.PDBFile.writeFile(modeller.topology,modeller.positions,open('SolvatedParticle.pdb','w'))

# Creating everything necessary for openmm simulation:
system = forcefield.createSystem(modeller.topology,nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
if args.gpu:
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    simulation = Simulation(modeller.topology, system, integrator, platform, properties)
else:
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

t = time.time()
# Minimisation:
simulation.minimizeEnergy()
print 'System minimised in {0} seconds.'.format(time.time()-t)

t = time.time()
# Equilibrate:
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(500)
print 'System thermalised in {0} seconds.'.format(time.time()-t)

positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
PDBFile.writeFile(simulation.topology,positions, open('equil.pdb', 'w'))
