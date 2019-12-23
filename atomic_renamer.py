__author__ = "Matteo Ferla"
__doc__ = """
'AtomLabel
"""

import requests, copy
from rdkit import Chem
from rdkit.Chem import rdFMCS, Draw
from IPython.display import display
from typing import List, Dict


class AtomicNamer:
    """
    Given a molecule label it according to a reference ligand from the PDB.

    >>> atpnamer = AtomicNamer('ATP')
    >>> mol = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=N2)N')
    >>> labels = atpnamer.name(mol)


    """

    def __init__(self, ligand_code: str, mark_upmatched: bool = True):
        """
        Initialise the namer with the reference ligand

        :param ligand_code: PDB 3letter code for ligand. 4letter is fine too. But has to exist on RCSB PDB.
        :type ligand_code: str
        :param mark_upmatched: Add an X between the symbol and the number?
        :type mark_upmatched: bool
        """
        self.ligand_code = ligand_code  #: PDB 3/4letter code for ligand
        self.mark_upmatched = mark_upmatched  #: Add an X between the symbol and the number?
        ref = self.mol_from_3letter(ligand_code)
        reflabels = self.atomlabels_from_3letter(ligand_code)
        self.label(ref, reflabels)
        self.ref = ref  #: Mol of the ref (from PDB ligand)
        self.reflabels = reflabels  #: atom labels for ref

    def name(self, mol: Chem.rdchem.Mol) -> List:
        """
        Given a molecule labelled it according to the reference ligand.

        :param mol: the mol to label
        :type mol: Chem.rdchem.Mol
        :return: atom labels
        :rtype: List[str]
        """
        common = Chem.MolFromSmarts(rdFMCS.FindMCS([mol, self.ref]).smartsString)
        commonlabels = [self.reflabels[a] for a in self.ref.GetSubstructMatches(common)[0]]
        self.label(common, commonlabels)
        mollabels_dict = {a: commonlabels[i] for i, a in enumerate(mol.GetSubstructMatches(common)[0])}
        mollabels = self.complete_labels(mol, mollabels_dict, self.mark_upmatched)
        self.label(mol, mollabels)
        return mollabels

    @staticmethod
    def _get_from_PDB(ligand_code: str, frmt: str) -> str:
        """
        Get the ligand ``ligand_code`` in the format ``frmt``

        :param ligand_code: RCSB PDB code
        :type  ligand_code: str
        :param frmt: ``_model.sdf``, ``_ideal.sdf`` or ``.cif``
        :type  frmt: str
        :return: the text of the molecule
        :rtype" str
        """
        return requests.get(f'http://files.rcsb.org/ligands/view/{ligand_code.lower()}{frmt}').text

    def mol_from_3letter(self, ligand_code: str) -> Chem.rdchem.Mol:
        """
        Get the ligand ``ligand_code`` as an sdf (mol) block.

        :param ligand_code: RCSB PDB code
        :type  ligand_code: str
        :return: the molecule
        :rtype" Chem.rdchem.Mol
        """
        sdf = self._get_from_PDB(ligand_code, '_model.sdf')
        return Chem.MolFromMolBlock(sdf)

    def atomlabels_from_3letter(self, ligand_code: str, dehydrogenated: bool = True) -> List:
        """
        The index of the list is the index of atoms.

        :param ligand_code: RCSB PDB code
        :type  ligand_code: str
        :param dehydrogenated: remove the hydrogens?
        :type dehydrogenated: bool
        :return: list of atomlabels
        """
        cif = self._get_from_PDB(ligand_code, '.cif')
        namer = lambda row: row[3:9].strip().replace('\"', "")
        ref_atoms = [namer(row) for row in cif.split('\n') if len(row.split()) > 15]
        if dehydrogenated:
            ref_atoms = self.dehydrogenate(ref_atoms)
        return ref_atoms

    @staticmethod
    def dehydrogenate(atomlabels: List) -> List:
        """
        The sdf will lack expicit hydrogens. So the list needs dehydrogenating.
        dodgy way of doing it as all etas will be lost.

        :param atomlabels: atom labels
        :type atomlabels: List[str]
        :return: the atomlabels with no hydrogens.
        :rtype: List[str]
        """
        return [a for a in atomlabels if 'H' not in a]

    @staticmethod
    def label(mol: Chem.rdchem.Mol, atomlabels: List) -> None:  # -> mol inplace.
        """
        Assign the prop ``AtomLabel``... https://www.rdkit.org/docs/RDKit_Book.html

        :param mol: the molecule to be labelled _in place_.
        :type mol: Chem.rdchem.Mol
        :param atomlabels: atom labels
        :type atomlabels: List[str]
        :return: None
        """
        assert len(
            atomlabels) == mol.GetNumAtoms(), 'the number of atoms in mol has to be the same as atomlabels. Hydrogens? dehydrogenate!'
        for idx in range(mol.GetNumAtoms()):
            mol.GetAtomWithIdx(idx).SetProp('AtomLabel', atomlabels[idx])
        return None

    @staticmethod
    def complete_labels(mol: Chem.rdchem.Mol, mollabels_dict: Dict, mark_upmatched: bool = True) -> List:
        """
        Complete the gaps in the atom labels dictionary (normally a list), by given names like CX1.

        :param mol: the molecule to be labelled _in place_.
        :type mol: Chem.rdchem.Mol
        :param mollabels_dict: key is index (int) and value is name like for a normal atomlabels (but with gaps)
        :type mollabels_dict: Dict
        :param mark_upmatched: Add an X between the symbol and the number
        :type mark_upmatched: bool
        :return: atom labels
        :rtype: List[str]
        """
        mollabels = []
        counters = {}
        for i in range(mol.GetNumAtoms()):
            if i in mollabels_dict:
                mollabels.append(mollabels_dict[i])
            else:
                el = mol.GetAtomWithIdx(i).GetSymbol().upper()
                if el in counters:
                    counters[el] += 1
                else:
                    counters[el] = 1
                if mark_upmatched:
                    mollabels.append(f'{el}X{counters[el]}')
                else:
                    mollabels.append(el + str(counters[el]))
        return mollabels

    @staticmethod
    def display(mol: Chem.rdchem.Mol, show='name'):
        # show = 'index' | 'name'
        if show:
            atoms = mol.GetNumAtoms()
            mol = copy.deepcopy(mol)
            for idx in range(atoms):
                if show == 'index':
                    mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(idx))
                elif show == 'name':
                    raise NotImplementedError(
                        'I need to figure out what property is needed as molAtomMapNumber is an str(int)')
                    mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',
                                                    str(mol.GetAtomWithIdx(idx).GetProp('AtomLabel')))
                else:
                    raise ValueError
        display(Draw.MolToImage(mol))
        return None

def test():
    atpnamer = AtomicNamer('ATP')
    mol = Chem.MolFromSmiles('C1=NC2=C(N1)C(=O)NC(=N2)N')
    labels = atpnamer.name(mol)
    assert labels == ['C8', 'N9', 'C4', 'C5', 'N7', 'C6', 'OX1', 'N6', 'C2', 'N3', 'N1']
    print(labels)
