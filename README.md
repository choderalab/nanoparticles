# nanoparticles

Scripts and files for nanoparticle simulations.

## Manifest


```
name_list.txt:   the naming scheme of the molecules and simulations.
```

`DLSdata`
```
Data and analysis for high throughput constant pH experiment conducted by Mehtap. 
```

`ImplicitAssembly`
```
Scripts and structures from the assembly of nanoparticles in implicit solvent.
```

`ExplicitProduction`
```
Unbiased molecular dynamics of the model nanoparticles that were created in implicit solvent
```

`AnalysisTools`
```
Python scripts for nanoparticle analysis.
```

`OriginalFiles/`

```
    The original mol2 files of drugs and dye molecules sent to Gregory Ross by Yosef Shamay.
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

`GrantFigures`
```
Figures produced for the collaborative grant with the Heller Lab.
```
