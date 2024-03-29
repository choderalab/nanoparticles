## Explicit solvent simulations of models of nanoparticles

The folders `1e`, `2m`, ..., `10v` contain the run scripts and structures from the simulations. Within each folder 
* `Thermalized_bigger.pdb` was created with `run_equil.py` using structures from `../ImplicitAssembly/Nanoparticles`.
* `submit` is the submission script for Hal that calls `run_explicit.py`.
* `snapshots.pdb` are frames taken from the simulation.
* `nanoparticle.pdb` is the centered nanoparticle from the final frame of the simulation with all water molecules and ions stripped away.

### Manifest

```
run_equil.py:         The script used to solvate and thermalize the model nanoparticle from the implicit solvent simulations.
```

```
run_explicit.py:      The script used for the explicit solvent simulations, based on the output from run_equil.py.
```

```
GetNanoParticles.py:  Script to extract final frame of nanoparticle simulations. 
```

```
Analysis              Folder use to analyze nanoparticle simulations based on the scripts in ../AnalysisTools.
```
