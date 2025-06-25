
# === Standard Library Imports ===
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from PIL import Image  # Needed to create blank image
import pandas as pd
import os
from typing import Optional

from app.config import Config 
import random
import base64
from io import BytesIO

   

class MoleculeInfo:

    def __init__(self, role='placeholder', ID=000, smiles='*'):
        self.img_dir = Config.IMG_DIR
        self.role = role
        self.ID = ID
        self.smiles = smiles
        self.mol = self._safe_mol_from_smiles(smiles)
        self.img_name = None
        self.image = None

        if self.mol:
            self.img_name, self.image = self._generate_and_save_image()
        else:
            self.img_name, self.image = self._generate_blank_image()

    def __repr__(self):
        return f"MoleculeInfo(role='{self.role}', ID='{self.ID}', smiles='{self.smiles}')"

    def _safe_mol_from_smiles(self, smiles):
        try:
            return Chem.MolFromSmiles(smiles)
        except Exception as e:
            print(f"Invalid SMILES for {self.ID}: {e}")
            return None

    def update_role(self, new_role):
        """Update the molecule's role."""
        self.role = new_role

    def update_role_with_bridge_position(self, bridge_position):
        """Update the molecule's role."""
        self.role = f"{self.role} (Bridge Position: {bridge_position})"

    def _generate_and_save_image(self):
        try:
            os.makedirs(self.img_dir, exist_ok=True)

            unique_id = f"{random.randint(0, 999999):06d}" 
            filename = f"{unique_id}_{self.ID}.png"
            filepath = os.path.join(self.img_dir, filename)

            img = Draw.MolToImage(self.mol, size=(300, 300))
            img.save(filepath)

            return filename, img

        except Exception as e:
            print(f"Failed to generate image for {self.ID}: {e}")
            return self._generate_blank_image()

    def _generate_blank_image(self):
        img = Image.new('RGB', (300, 300), color='white')
        filename = f"{self.ID}_blank.png"
        filepath = os.path.join(self.img_dir, filename)

        try:
            os.makedirs(self.img_dir, exist_ok=True)
            img.save(filepath)
        except Exception as e:
            print(f"Failed to save blank image: {e}")
            filename = None

        return filename, img

    def _generate_image_base64(self):
        img = Draw.MolToImage(self.mol, size=(200, 200))
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        encoded = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return None, f"data:image/png;base64,{encoded}"


    def update_with_new_smiles(self, new_smiles):
        self.smiles = new_smiles
        self.mol = self._safe_mol_from_smiles(new_smiles)
        if self.mol:
            self.img_name, self.image = self._generate_and_save_image()
        else:
            self.img_name, self.image = self._generate_blank_image()


def get_random_molecule_info(
    df: pd.DataFrame,
    synthon_number: int,
    role: str,
    id_col: str = 'synton_id',
    smiles_col: str = 'SMILES',
    synthon_col: str = 'synton#'
    ) -> Optional[MoleculeInfo]:
    """
    Returns a MoleculeInfo object with given role and id based on a random SMILES.
    """
    if synthon_col not in df.columns or smiles_col not in df.columns or id_col not in df.columns:
        return MoleculeInfo()
        
    filtered = df[df[synthon_col] == synthon_number]
    if filtered.empty:
        return MoleculeInfo()

    random_row = filtered.sample(n=1).iloc[0]
    smiles = random_row[smiles_col]
    mol_id = random_row[id_col]
    
    return MoleculeInfo(role=role, ID=mol_id, smiles=smiles)


def generate_product(
    rxn_smarts: str,
    mol_infos: list[MoleculeInfo],
    reaction_id: str
) -> Optional[MoleculeInfo]:
    """
    Generate product molecule from reaction SMARTS and input synthon MoleculeInfo list.
    Returns a MoleculeInfo object for the product or None if failed.
    """
    # Extract RDKit Mol objects from MoleculeInfo, skipping any None

    if mol_infos[2].smiles == '*':
        reactants = (mol_infos[0].mol, mol_infos[1].mol)
    else:
        reactants = (mol_infos[0].mol, mol_infos[1].mol, mol_infos[2].mol, )

    rxn_mol = AllChem.ReactionFromSmarts(rxn_smarts)
    try:
        product_mols = rxn_mol.RunReactants(reactants)
        if not product_mols:
            return None
        product = product_mols[0][0]
        Chem.SanitizeMol(product)
        ids = [m.ID for m in mol_infos]
        product_id = f"{reaction_id}_{'_'.join(map(str, ids))}"
        product_smiles = Chem.MolToSmiles(product)

        return MoleculeInfo(role="product", ID=product_id, smiles=product_smiles)

    except Exception as e:
        print(f"Failed to generate product for {reaction_id}: {e}")
        return MoleculeInfo(role="product", ID=None, smiles=None)

