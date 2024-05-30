from IPython.display import clear_output
from io import StringIO
import ase
from ase.io import read
from ase.calculators.mopac import MOPAC
from ase.optimize import BFGS

try: from openbabel import pybel
except: pass

try: from rdkit.Chem import AllChem
except: pass

try: from rdkit import Chem
except: pass

try: from xtb_ase import XTB
except: pass

try: import torchani
except: pass

try: from aimnet2ase import aimnet2_optimize
except: pass

try: from tblite.ase import TBLite
except: pass

def smiles2mol(smiles:str, N_conformers:int=50)->str:
  """
    get mol format block from smiles
  """
  # smiles to mol block
  mol = Chem.MolFromSmiles(smiles)
  mol = Chem.AddHs(mol)

  # sanitize
  Chem.SanitizeMol(mol)

  # generate 3D embedded conformers
  num_conformers = N_conformers
  params = AllChem.ETKDGv3() # embedding algorithm
  params.randomSeed = 42

  conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)

  energies = []
  for conf_id in conformer_ids:
    # MMFF94 힘장 설정
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)

    # force field optimize
    ff.Minimize()
    energy = ff.CalcEnergy()
    energies.append((conf_id, energy))

  # get lowest conformer
  lowest_energy_idx, _ = min(energies, key=lambda x: x[-1])

  # mol format block
  mol_block = Chem.MolToMolBlock(mol, confId=lowest_energy_idx)

  return mol_block



