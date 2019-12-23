# AtomicRenamer

Given a molecule label the atoms (names/labels) according to a reference ligand from the PDB.

## Why?

RDKit is not good with atom names, while macromolecular tools like Rosetta or even PyMOL rely heavily on them. This is aimed at fixing that. An atom has several properties, here I will mention:

* symbol (or element name), say C for calcium
* index (an integer specifying the order it is in which often comes from the order they appeared in a SMILES string)
* labels (or atom name), say CA for C-&alpha;

Open babel if given a sdf of an amino acid and asked to convert to a mol2 will label them as CA etc. But most other times and most programs don't.

This stems also from the fact that mol (sdf) files do not specify atom labels in the main block. Mol2 and PDB do however.

## Script
This short script given a molecule (_e.g._ `mol = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=N2)N')`) and reference PDB ligand code (_e.g._ `ATP`) will label in place the molecule (adding the property `AtomLabel`) and return a list of atom names (with indices matching the atomic indices obvious).

    >>> from rdkit import Chem
    >>> from atomic_renamer import AtomicNamer
    >>> mol = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=N2)N')
    >>> labels = AtomicNamer('ATP').name(mol)
    >>> labels
    ['C8', 'N9', 'C4', 'C5', 'N7', 'C6', 'OX1', 'N6', 'C2', 'N3', 'N1']
    
The atom name/label is assigned to the prop `AtomLabel` (following https://www.rdkit.org/docs/RDKit_Book.html).

Note, that while there is a bound method called `.display(mol)`, I have not finished it as I been able to get it to work. 
These labels can be saved as `mol2` in a convoluted way, becuase the mol2 writer in Rdkit is a bit tempramental. So using open babel is better and using the bound method `AtomicNamer.fix` to fix these.

    >>> mol.UpdatePropertyCache() # I might have changed some atoms around
    >>> mol = Chem.AddHs(mol) #protonate explicitly
    >>> Chem.GetSSSR(mol) #not communists, but resonance fixing
    >>> AllChem.EmbedMolecule(mol) #initialise for 3d.
    >>> AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    >>> Chem.MolToMolFile(mol, 'guanine.mol')
    >>> os.system(f"obabel -i mol guanine.mol -o mol2 -O guanine.mol2")
    >>> AtomicNamer.fix('guanine.mol2', 'guanine.better.mol2', labels)
    >>> os.system(f"obabel -i mol2 guanine2.better.mol2 -o mol2 -O guanine.conf.mol2 --conformer --nconf 30 --writeconformers")

Once this is done, the mol2 can be used. If using Rosetta and are about to parametrise it, why not check out my [2to3 port of mol_to_params.py](https://github.com/matteoferla/mol_to_params.py)?

Also, for more stuff, see [my blog post about Rdkit](https://blog.matteoferla.com/2019/10/rdkit-for-rosetta-plp-ligand-space-as.html).

## Under the hood

The attributes `.ref` and `.reflabels` contain the RDKit `Chem.rdchem.Mol` object and the list of atom names. So if you want to use something that isn't a PDB ligand code you can.

The reason for using the PDB ligand code is that if you change the name of a residue (in TextEdit or PyMOL) and run the structure through Rosetta will change it. This is handy for post translation modifications —for more see [my blog post about PTMs and Rosetta](https://blog.matteoferla.com/2019/01/phosphorylated-pdb-files.html).

