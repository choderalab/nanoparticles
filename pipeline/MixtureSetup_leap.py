from openmoltools import openeye, amber         # for Choderalab wrappers of openeye software
from openmoltools.utils import enter_temp_directory
import openeye.oechem as oechem
from openmoltools import packmol
import parmed as pmd
import subprocess
import os
import shutil

def oemol_to_antechamber(m, gaff_mol2_filename, frcmod_filename, residue_name="MOL", strictStereo=False):
    """
        Build a molecule from a mol2 file and run antechamber,
        generating GAFF mol2 and frcmod files from a smiles string.  Charges
        will be generated using the OpenEye QuacPac AM1-BCC implementation.
        
        Created by hacking openmoltools/openeye.py
        
        Parameters
        ----------
        m : oechem molecule object
        Molecule to construct and charge
        gaff_mol2_filename : str
        Filename of mol2 file output of antechamber, with charges
        created from openeye
        frcmod_filename : str
        Filename of frcmod file output of antechamber.  Most likely
        this file will be almost empty, at least for typical molecules.
        residue_name : str, optional, default="MOL"
        
        OpenEye writes mol2 files with <0> as the residue / ligand name.
        This chokes many mol2 parsers, so we replace it with a string of
        your choosing.  This might be useful for downstream applications
        if the residue names are required to be unique.
        strictStereo : bool, optional, default=False
        If False, permits smiles strings with unspecified stereochemistry.
        See https://docs.eyesopen.com/omega/usage.html
    """
    #oechem = import_("openeye.oechem")
    #if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for oechem!"))
    
    # Get the absolute path so we can find these filenames from inside a temporary directory.
    gaff_mol2_filename = os.path.abspath(gaff_mol2_filename)
    frcmod_filename = os.path.abspath(frcmod_filename)
    
    m = openeye.get_charges(m, strictStereo=strictStereo, keep_confs=1)
    
    with enter_temp_directory():  # Avoid dumping 50 antechamber files in local directory.
        _unused = openeye.molecule_to_mol2(m, "./tmp.mol2", residue_name=residue_name)
        net_charge = oechem.OENetCharge(m)
        tmp_gaff_mol2_filename, tmp_frcmod_filename = amber.run_antechamber("tmp", "./tmp.mol2", charge_method=None, net_charge=net_charge)  # USE OE AM1BCC charges!
        shutil.copy(tmp_gaff_mol2_filename, gaff_mol2_filename)
        shutil.copy(tmp_frcmod_filename, frcmod_filename)


def parametrize_mol(inname,resname):
    """
        Parameterizes molecule with antechamber and adds charges using Open Eye's am1-bcc method from a mol2 file.
        
        Parameters
        ----------
        inname  : str
        the names of the mol2 file of the molecule
        resname : str
        the name you wish to assign as the residue name
        
        Returns
        -------
        nothing
        
        BUT 'resname'_gaff.mol2 and 'resname'.frcmod are created in the working directory.
    """
    mol = oechem.OEGraphMol()           # initialising molecule object
    ifs = oechem.oemolistream()         # initialising input stream for reading in data
    ifs.SetFormat(oechem.OEFormat_MOL2) # specifying mol2 format input
    ifs.open(inname)
    if oechem.OEReadMolecule(ifs,mol):  # this function automatically returns True or False, to help spot for errors.
        pass
    else:
        print "Problem loading molecule!"
    mol = oechem.OEMol(mol)
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)
    oechem.OEAddExplicitHydrogens(mol)
    mol.SetTitle(resname)
    if any([atom.GetName() == '' for atom in mol.GetAtoms()]):
        oechem.OETriposAtomNames(mol)
    oemol_to_antechamber(mol,resname+'_gaff.mol2',resname+'.frcmod',residue_name=resname)

def mol2pdb(inname):
    outname = inname.split('.')[0] + '.pdb'
    parm = pmd.load_file(inname,structure=False)
    parm.save(outname,overwrite=True)

def oemol2pdb(mol2name,pdbname):
    mol = oechem.OEGraphMol()           # initialising molecule object
    ifs = oechem.oemolistream()         # initialising INPUT stream for reading in data
    ofs = oechem.oemolostream()         # initialising the OUTPUT stream for writing data
    # Reading:
    ifs.SetFormat(oechem.OEFormat_MOL2) # specifying mol2 format input
    ifs.open(mol2name)
    if oechem.OEReadMolecule(ifs,mol):
        pass
    else:
        print "Problem loading molecule!"
    # Writing:
    ofs.SetFormat(oechem.OEFormat_PDB)      # specifying that I want a pdb format
    ofs.open(pdbname)
    oechem.OEWriteMolecule(ofs,mol)


def tleap_nanosim(drugname,dyename,boxname,solvate=True):
    tleap_input = """
        source leaprc.ff99SB
        source leaprc.gaff
        {0} = loadmol2 {0}_gaff.mol2
        {1} = loadmol2 {1}_gaff.mol2
        loadamberparams {0}.frcmod
        loadamberparams {1}.frcmod
        {2} = loadpdb {2}.pdb
        solvateBox {2} TIP3PBOX 0.5
        addIons2 box Na+ 0
        saveamberparm {2} out.prmtop out.inpcrd
        quit
        """.format(drugname,dyename,boxname)
    file_handle = open('tleap_commands', 'w')
    file_handle.writelines(tleap_input)
    file_handle.close()
    
    cmd = "tleap -f {0}".format(file_handle.name)
    subprocess.call(cmd,shell=True)

def tleap_mixture(resnames,boxname,solvate=True,implicit=False,run=True):
    """
        Creates and runs tleap commands for setting up mixture simulations of small molecules.
        The mixture can be solvated with water.
        
        Parameters
        ----------
        resnames  : list of str
        the names of the molecules in the mixture. Must match residue name and antechamber mol2 file.
        i.e. resname[1]= 'LIG', implies the residue name is LIG and there is an antechamber mol2 file 'LIG_gaff.mol2'.
        boxname : str
        the name of pdb file (without extension .pdb) file that contains the mixture of molecules in resnames
        solvate = bool
        whether to fill the free space in boxname.pdb with TIP3P water.
        The box is neutralized assuming there is a net negative charge.
        run: bool
        whether to run tleap
        
        Returns
        -------
        tleap_input : str
        the input file used for tleap
        """
    tleap_input = "source leaprc.ff99SB\nsource leaprc.gaff\n"
    for name in resnames:
        tleap_input += "{0} = loadmol2 {0}_gaff.mol2\n".format(name)
        tleap_input += "loadamberparams {0}.frcmod\n".format(name)
    tleap_input += "{0} = loadpdb {0}.pdb\n".format(boxname)
    tleap_input += "setbox {0} centers\n".format(boxname)
    if solvate == True and implicit == False:
        tleap_input += "solvateBox box TIP3PBOX 0.5\n"
    #tleap_input += "addIons2 box Na+ 0\n"       # Not necessary if ions are already added by packmol
    elif solvate == True and implicit == True:
        tleap_input += "set default PBRadii mbondi2\n"  # Set for IGB=5 radii
    elif implicit == True:
        tleap_input += "set default PBRadii mbondi2\n"  # Set for IGB=5 radii
    tleap_input += "saveamberparm {0} out.prmtop out.inpcrd\nquit".format(boxname)
    if run == True:
        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleap_input)
        file_handle.close()
        cmd = "tleap -f {0}".format(file_handle.name)
        subprocess.call(cmd,shell=True)
    return tleap_input

def nanopipeline(drugmol2,dyemol2,drugname='LIG',dyename='DYE',ndrugs=10,ndye=20,boxsize=50,neutralize=True,implicit=False):
    """
        Starting with the mol2 files of drug and dye, the molecules are parametrised, assembled uniformly in
        cubic box, and AMBER prmtop and inpcrd files are created for the system.
        
        Parameters
        ----------
        drugmol2, dyemol2  : str, str
        the mol2 file name of the drug/dye, with or without hydrogens
        drugname, dyename  : str, str
        the residue name you want for the drug/dye residue name and output file naming schemes
        ndrugs/ndye: int
        the number of drug/dye molecules you want in the box
        boxsize: int
        the length of the cubic box in Angstroms
        
        Returns
        -------
        out.inpcrd and out.out.prmtop are created in the working directory, as well gaff output files
        """
    
    # Parmetrizing with antechamber via openeye:
    parametrize_mol(drugmol2,drugname)
    parametrize_mol(dyemol2,dyename)
    # Creating box of drugs and dye. Packmol requires PDB files
    mol2pdb(drugname+'_gaff.mol2')
    mol2pdb(dyename+'_gaff.mol2')
    if neutralize == True:
        box = packmol.pack_box([drugname+'_gaff.pdb',dyename+'_gaff.pdb','Na.pdb'],[ndrugs,ndye,ndye], tolerance=2.0, box_size=boxsize)
    else:
        box = packmol.pack_box([drugname+'_gaff.pdb',dyename+'_gaff.pdb'],[ndrugs,ndye], tolerance=2.0, box_size=boxsize)
    box.save('box.pdb')
    #os.remove(drugname+'_gaff.pdb')
    #os.remove(dyename+'_gaff.pdb')
    # Running tleap to create out.prmtop out.inpcrd for use with openmm
    #tleap_nanosim(drugname,dyename,'box')
    if implicit == False:
        tleap_mixture([drugname,dyename],'box',solvate=True,implicit=False,run=True)
    else:
        tleap_mixture([drugname,dyename],'box',solvate=False,implicit=True,run=True)

def mixture_setup(molfiles,resnames,nmols,boxsize,solvate=True):
    """
        Function that generalises the function 'nanopipeline' for general mixtures.
        Starting with the mol2 files of drug and dye, the molecules are parametrised, assembled uniformly in
        cubic box, and AMBER prmtop and inpcrd files are created for the system.
        
        Parameters
        ----------
        molfiles : list of str
        the mol2 file name of the compounds, with or without hydrogens
        resnames : list of str
        the residue names you want for the molecules  residue name and output file naming schemes
        nmols : list of int
        the number of molecules you want in the box, excluding water and ions
        boxsize: int
        the length of the cubic box in Angstroms
        
        Returns
        -------
        out.inpcrd and out.out.prmtop are created in the working directory, as well gaff output files
        """
    
    # Parmetrizing with antechamber via openeye:
    for i in range(len(molfiles)): parametrize_mol(molfiles[i],resnames[i])
    # Creating box of drugs and dye. Packmol requires PDB files
    for i in range(len(molfiles)): mol2pdb(resnames[i]+'_gaff.mol2')
    molfiles_pdb = [mol + '_gaff.pdb' for mol in resnames]                   
    box = packmol.pack_box(molfiles_pdb, nmols, tolerance=2.0, box_size=boxsize)
    box.save('box.pdb')
    # Running tleap to create out.prmtop out.inpcrd for use with openmm
    tleap_input = tleap_mixture(resnames,'box',solvate=solvate,run=True)