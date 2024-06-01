from IPython.display import clear_output
from io import StringIO
import ase
from ase.io import read
from ase.calculators.mopac import MOPAC
from ase.optimize import BFGS
from rdkit import Chem
from rdkit.Chem import AllChem

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

def rdkit_optimize(smiles_mol:str, charge:int, forcefield:str, mmffVariant=None, smiles_input=True)->str:
  """
  Description
  -----------
  rdkit geometry optimize using force field

  Parameters
  ----------
    - smiles_mol(str) : smiles or mol string
    - charge(int) : molecular total charge
    - forcefield(str) : supported force field : ['MMFF', 'UFF']
    - mmffVariant(str) : MMFF forcefield variants : ['MMFF94', 'MMFF94s']
    - smiles_input(bool) : whether smiles_mol is smiles or mol format

  Returns
  -------
    - optimized xyz string(str)
  """
  # forcefield parameter check
  if forcefield == "MMFF" and mmffVariant not in ["MMFF94", "MMFF94s"]: raise ValueError(f'mmffVariant expected "MMFF94" or "MMFF94s", but got {mmffVariant}')

  # read smiles
  if smiles_input:
    mol_block = smiles2mol(smiles_mol)
  else:
    mol_block = smiles_mol
  mol = Chem.MolFromMolBlock(mol_block, removeHs=False)

  # charge check
  Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
  partial_chargies = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms()]
  charge_from_molecule = round(sum(partial_chargies))
  charge_from_input = charge

  if charge_from_molecule != charge_from_input:
    print("Incosistency of charge")
    return None

  # optimize
  if forcefield == 'UFF':
      status = AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
  elif forcefield == 'MMFF':
      status = AllChem.MMFFOptimizeMolecule(mol, mmffVariant, maxIters=2000)
  else:
      raise ValueError(f"{forcefield} is supported. Use 'UFF' or 'MMFF'.")

  # optimize status check
  if status == 1:
     print("Optimizer : Not converged in 2000 iteractions.")
     return None
  elif status == -1:
     print("Optimizer : Fail to optimize, status code : -1")
     return None

  optimized_xyz = Chem.MolToXYZBlock(mol)
  return optimized_xyz

def openbabel_optimize(smiles_mol:str, charge:int, forcefield:str, smiles_input:bool=True)->str:
  """
  Description
  -----------
  openbabel geometry optimize using force field

  Parameters
  ----------
    - smiles_mol(str) : smiles or mol string
    - charge(int) : molecular total charge
    - forcefield(str) : supported force field : ['uff', 'mmff94', 'mmff94s', 'ghemical', 'gaff']
    - smiles_input(bool) : whether smiles_mol is smiles or mol format

  Returns
  -------
    - optimized xyz string(str)
  """
  # read smiles
  if smiles_input:
    mol_block = smiles2mol(smiles_mol)
  else:
    mol_block = smiles_mol
  molecule = pybel.readstring(format="mol", string=mol_block)

  # set charge
  charge_from_molecule = molecule.charge
  charge_from_input = charge

  if charge_from_molecule != charge_from_input:
    print("Incosistency of charge")
    return None

  # optmize using FF
  molecule.localopt(forcefield=forcefield)
  optimized_xyz = molecule.write("xyz")
  return optimized_xyz

def mopac_optimize(xyz_string:str, charge:int, method:str, clear_log=True)->str:
  """
  Description
  -----------
  xtb geometry optimize 함수

  Parameters
  ----------
  - xyz_string (str) : xyz format string
  - charge (int) : molecular total chage
  - clear_log (bool) : clear optimization logging
  - method (str) : supported semi-empirical method
    - supported methods = ['AM1', 'MNDO', 'MNDOD', 'PM3', 'PM6', 'PM6-D3',
                           'PM6-DH+', 'PM6-DH2', 'PM6-DH2X', 'PM6-D3H4',
                           'PM6-D3H4X', 'PMEP', 'PM7', 'PM7-TS', 'RM1']

  Returns
  -------
  - optimized xyz (str)
  """

  # convert xyz format --> ase.Atoms
  mol = ase.io.read(StringIO(xyz_string), format="xyz")

  # set XTB calculator
  mol.calc = MOPAC(method=method, charge=charge )

  # geometry optimize
  optimizer = BFGS(mol)
  optimizer.run()

  # get xyz format string
  with StringIO() as output:
    ase.io.write(output, mol, format="xyz")
    opt_xyz = output.getvalue()

  # clear cell output
  if clear_log:
    clear_output()

  return opt_xyz

def torchani_optimize(xyz_string:str, charge:int, model:str, clear_log=True)->str:
  """
  Description
  -----------
  torchani geometry optimize 함수

  Parameters
  ----------
  - xyz_string (str) : xyz format string
  - charge (int) : molecular total chage
  - model (str) : torchani ML potential model
  - clear_log (bool) : clear optimization logging

  Supported models
  ----------------
  - ani1ccx :  CCSD(T)*/CBS (DPLNO-CCSD(T)) | HCNO
  - ani1x   :  wB97X/6-31G(d)               | HCNO
  - ani2x   :  wB97X/6-31G(d)               | HCNOFSCl

  Returns
  -------
  - optimized xyz (str)
  """
  # convert xyz format --> ase.Atoms
  mol = ase.io.read(StringIO(xyz_string), format="xyz")

  # check compatibility between the torchani model and the elements in the molecule
  model = model.lower()
  supporting_elements = {
      "ani1ccx" : {'H', 'C', 'N', 'O'},
      "ani1x"   : {'H', 'C', 'N', 'O'},
      "ani2x"   : {'H', 'C', 'N', 'O', 'F', 'S', 'Cl'}}

  elements_in_mol = set(mol.get_chemical_symbols())
  unsupported_elements = elements_in_mol - supporting_elements[model]
  if len(unsupported_elements) != 0:
    print(f"{unsupported_elements} are not compatible with the {model} model")
    return None

  # torchani supports only neutral compounds
  if charge != 0:
    print("torchani supports only neutral molecule")
    return None

  # set torchani calculator
  if model == "ani1ccx":
    mol.calc = torchani.models.ANI1ccx().ase()
  elif model == "ani1x":
    mol.calc = torchani.models.ANI1x().ase()
  elif model == "ani2x":
    mol.calc = torchani.models.ANI2x().ase()
  else:
    raise ValueError("Not available model")

  # geometry optimize
  optimizer = BFGS(mol)
  optimizer.run()

  # get xyz format string
  with StringIO() as output:
    ase.io.write(output, mol, format="xyz")
    opt_xyz = output.getvalue()

  # clear cell output
  if clear_log:
    clear_output()

  return opt_xyz

def xtb_optimize(xyz_string:str, charge:int, method:str,
                 spinpol:bool|None=None, uhf:int=0, clear_log=True)->str:
  """
  Description
  -----------
  xtb geometry optimize 함수

  Parameters
  ----------
  - xyz_string (str) : xyz format string
  - charge (int) : molecular total chage
  - method (str) : xtb method
    - supported methods : ["gfn1-xtb", "gfn2-xTB", "gfn-ff"]
  - spinpol ( bool | None ) : number of unpaired electrons
  - uhf (int) : whether to use spin-polarized xTB ( see :  https://github.com/Andrew-S-Rosen/xtb_ase/blob/main/src/xtb_ase/calculator.py )
  - clear_log (bool) : clear optimization logging

  Returns
  -------
  - optimized xyz (str)
  """

  # convert xyz format --> ase.Atoms
  mol = ase.io.read(StringIO(xyz_string), format="xyz")

  # set XTB calculator
  mol.calc = XTB(method=method, charge=charge, spinpol=spinpol, uhf=uhf)

  # geometry optimize
  optimizer = BFGS(mol)
  optimizer.run()

  # get xyz format string
  with StringIO() as output:
    ase.io.write(output, mol, format="xyz")
    opt_xyz = output.getvalue()

  # clear cell output
  if clear_log:
    clear_output()

  return opt_xyz

def tblite_optimize(xyz_string:str, charge:int, method:str, clear_log=True)->str:
  """
  Description
  -----------
  TBLite geometry optimize 함수

  Parameters
  ----------
  - xyz_string (str) : xyz format string
  - charge (int) : molecular total chage
  - method (str) : xtb method
    - supported methods : ["GFN1-xTB", "GFN2-xTB", "IPEA1-xTB"]
  - clear_log (bool) : clear optimization logging

  Returns
  -------
  - optimized xyz (str)
  """

  # convert xyz format --> ase.Atoms
  mol = ase.io.read(StringIO(xyz_string), format="xyz")

  # set TBLite calculator
  mol.calc = TBLite(method=method, charge=charge)

  # geometry optimize
  optimizer = BFGS(mol)
  optimizer.run()

  # get xyz format string
  with StringIO() as output:
    ase.io.write(output, mol, format="xyz")
    opt_xyz = output.getvalue()

  # clear cell output
  if clear_log:
    clear_output()

  return opt_xyz



