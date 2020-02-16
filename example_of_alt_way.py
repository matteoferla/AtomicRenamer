__doc__ = """
This is an example of the normal approach.
"""

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign
import pymol2
from IPython.display import display

with pymol2.PyMOL() as pymol:
    pymol.cmd.fetch('3SRG')
    pymol.cmd.remove('solvent')
    pymol.cmd.remove('segi E+F+G+H+I+J')
    pymol.cmd.save('3SRG.clean.pdb')
    pymol.cmd.save('3SRG.clean.lig_only.mol', 'resn OCH')
    pymol.cmd.save('3SRG.clean.lig_only.pdb', 'resn OCH')

#https://www.rdkit.org/docs/Cookbook.html
# The reference molecule
# OCH
ref = Chem.MolFromSmiles('O=c1ccc2ccccc2[nH]1')
# The PDB conformations
target = Chem.MolFromPDBFile('3SRG.clean.lig_only.pdb')#, sanitize=False, strictParsing=False)
#mol1 = Chem.MolFromMolFile('3SRG.clean.lig_only.mol')#, sanitize=False, strictParsing=False)
target = AllChem.AssignBondOrdersFromTemplate(ref, target)
print('from structure')
display(target)
#mol2 = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1') #'phenylacetate'
probe = Chem.MolFromSmiles('Oc1(O)ccc2ccccc2[nH]1')
Chem.AddHs(probe)
AllChem.EmbedMolecule(probe)
AllChem.UFFOptimizeMolecule(probe, maxIters=2000)
Chem.rdPartialCharges.ComputeGasteigerCharges(probe)
print('new')
display(probe)

### find what is common
res=Chem.rdFMCS.FindMCS([probe, target],
                        #matchValences=True,
                        atomCompare=Chem.rdFMCS.AtomCompare.CompareElements,
                        bondCompare=Chem.rdFMCS.BondCompare.CompareOrder
                       )
common = Chem.MolFromSmarts(res.smartsString)
print('Common')
display(common)

### Align them
overlap_target = target.GetSubstructMatch(common)
overlap_probe = probe.GetSubstructMatch(common)
atomMap = [(probe_at,target_at) for probe_at,target_at in zip(overlap_probe, overlap_target)]
print(atomMap)
rms = rdMolAlign.AlignMol(probe, target, atomMap=atomMap, maxIters=500)
print(rms)
Chem.MolToMolFile(probe, 'inter.aligned.mol')