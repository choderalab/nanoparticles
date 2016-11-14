from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import time

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run an openmm simulation with salt exchange moves.")
    parser.add_argument('-b','--boxlength',type=float,help="the length of the simulation box in Angstroms, default=160.0",default=160.0)
    parser.add_argument('-n','--nmols',type=int,nargs = 2, help="the number of molecules in the mixture")
    parser.add_argument('-m','--molnames',type=str,nargs = 2, help="the name of the mol2 file names")
    parser.add_argument('-l','--length',type=float, help="the length of the simulation in nanometers, default = 100",default = 100)
    parser.add_argument("--implicit",action='store_true',help="whether the an implicit solvent model, default=False",default=False)
    parser.add_argument("--gpu",action='store_true',help="whether the simulation will be run on a GPU, default=False",default=False)
    parser.add_argument("--prep",action='store_true',help="whether only the preparation is done, default=False",default=False)
    parser.add_argument("--run",action='store_true',help="whether only the preparation is done, default=False",default=False)
    args = parser.parse_args()


molnames = args.molnames
nmols = args.nmols
boxsize =  args.boxlength

if args.prep:
    from MixtureSetup_leap import nanopipeline
    if args.implicit:
        nanopipeline(molnames[0],molnames[1],ndrugs=nmols[0],ndye=nmols[1],boxsize=boxsize,neutralize=False,implicit=True)
    else:
        nanopipeline(molnames[0],molnames[1],ndrugs=nmols[0],ndye=nmols[1],boxsize=boxsize,neutralize=True,implicit=False)

if args.run:
    prmtop = app.AmberPrmtopFile('out.prmtop')
    inpcrd = app.AmberInpcrdFile('out.inpcrd')

    # Converting into nanometers
    bxsize = quantity.Quantity(((boxsize/10.0,0.0,0.0),(0.0,boxsize/10.0,0.0),(0.0,0.0,boxsize/10.0)),unit=nanometer)
    prmtop.topology.setPeriodicBoxVectors(bxsize)

    if args.implicit:
        system = prmtop.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=1*nanometer,constraints=HBonds,implicitSolvent=OBC1)
        # Restraining the system in a spherical well
        fstring = '100*max(0, r-{0})^2; r=sqrt((x-{0})^2+(y-{0})^2+(z-{0})^2)'.format(boxsize/10.0/3.0)
        force = CustomExternalForce(fstring)
        system.addForce(force)
        for i in range(system.getNumParticles()):
            force.addParticle(i, [])
    else:
        system = prmtop.createSystem(nonbondedMethod=CutoffPeriodic,nonbondedCutoff=1*nanometer,constraints=HBonds)
        system.addForce(MonteCarloBarostat(1*atmospheres, 300*kelvin, 25))

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    if args.gpu:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    else:
        simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

    # Very quick minimisation:
    simulation.minimizeEnergy()

    # Equilibrate:
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    simulation.step(500)

    positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    PDBFile.writeFile(simulation.topology,positions, open('equil.pdb', 'w'))

    # Running simulation
    nsteps = int(1000*args.length/0.002)   # Converting desired length of simulation to number of steps
    nsnapshots = int(nsteps/1000.0)
    ndata  = int(nsteps/100.0)
    npdb = int(nsteps/5.0)

    simulation.reporters.append(app.DCDReporter('trajectory.dcd', nsnapshots))
    simulation.reporters.append(StateDataReporter('simdata.txt', ndata, step=True,potentialEnergy=True, temperature=True))
    simulation.reporters.append(PDBReporter('snapshots.pdb',npdb ))
    simulation.step(nsteps)