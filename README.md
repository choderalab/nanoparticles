# nanoparticles

## Manifest

'OriginalFiles/'

```
    The original mol2 files of drugs and dye molecules sent to Gregory Ross by Yosef Shanmay.
```

`CleanedFiles/`

```
    Mol2 files from OriginalFiles/ given cleaner names and ligand residue names given.
```

`ParamFiles/`

```
    Mol2 files from CleanedFiles parameterised with antechamber and parmchk.
    Imprecision in Ambermini's sqm code results in molecules with (slightly) non unit charge.
    This procedure will be replaced with parameteriztion with OpenEye from Openmoltools.
```

`Production_Implicit_OBC1_longer/`

```
    Packmol and leap scripts to prepare mixtures for 10 molecules. As well as run*.py scripts to perform Openmm
    implicit solvation simulations upto 300ns. Each simulation comprises of 50 drug molecules and 10
    dye molecules in a box of 161 by 1 Angstroms 
```



