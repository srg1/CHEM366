from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

largest_Fragment = rdMolStandardize.LargestFragmentChooser()

suppl = Chem.SupplierFromFilename("Take 3 with Names.sdf")

with Chem.SDWriter('Take 3 with Names.sdf') as w:
    for i,mol in enumerate(suppl):
        ID = mol.GetProp("Catalog ID")
        mol.SetProp("_Name",ID)
        mol = largest_Fragment.choose(mol)
        mol = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        try:
            mol.GetProp("_Name")
        except:
            mol.SetProp("_Name",ID)
        w.write(mol, confId=0)