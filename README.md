# AtomicRenamer
Given a molecule label the atoms (names/labels) according to a reference ligand from the PDB.

RDKit is not good with atom names, while tools like Rosetta rely heavily on them. This is aimed at fixing that. An atom has several properties, here I will mention:

* symbol (or element name), say C for calcium
* index (an integer specifying the order it is in which often comes from the order they appeared in a SMILES string)
* labels (or atom name), say CA for C-&alpha;

Open babel if given a sdf of an amino acid and asked to convert to a mol2 will label them as CA etc. But most other times and most programs don't.

This short script given a molecule (_e.g._ `mol = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=N2)N')`) and reference PDB ligand code (_e.g._ `ATP`) will label in place the molecule (adding the property `AtomLabel`) and return a list of atom names (with indices matching the atomic indices obvious).


  from rdkit import Chem
  from atomic_renamer import AtomicNamer
  mol = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=N2)N')
  AtomicNamer('ATP').name(mol)
  \['C8', 'N9', 'C4', 'C5', 'N7', 'C6', 'OX1', 'N6', 'C2', 'N3', 'N1']
   
The atom name/label is assigned to the prop `AtomLabel` (following https://www.rdkit.org/docs/RDKit_Book.html).

Note, that while there is a bound method called `.display(mol)`, I have not finished it as I been able to get it to work.   
