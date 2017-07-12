## Implicit Assembly of nanoparticle models

The entire procedure to create a nanoparticle from Mol2 files is automated with the run_nanomix.py script, which uses functions from the MixtureSetup_leap.py. The full details can be found within the scripts, but the creation of the nanoparticles models was as follows:
* Copies of drug and dye molecules were added uniformaly to a large box with packmol.
* The molecules were then parametrized with leap and antechamber
* Using the AMBER topology files, the systems were prepared to simulate in implicit solvent with OpenMM.
* To speed up the formation of the nanoparticles, a harmonic constraint was added to keep all the molecules within a central spherical region. Even without constraints, all molecules were found to form nanoparticle-like aggregates in implicit solvent, with the rate of formation being limited by diffusion times.
* The systems were then simulated with OpenMM for up to 300ns.

The folder `NanoStructures` contains the final structures of the simulations.
