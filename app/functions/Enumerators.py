 # === Standard Library Imports ===
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, AllChem, Draw
from PIL import Image
import pandas as pd
import re
import json
from pathlib import Path
import io
import random
from rdkit import Chem
from typing import Optional
import os
import io

from app.config import Config 
from app.functions.utills import extract_single_value
from app.functions.Molecules import get_random_molecule_info, generate_product, MoleculeInfo


def enumerator_2_component(object_1, object_2, rxn_mol, synthon_position, reaction_id):

    product_tuples = rxn_mol.RunReactants([object_1.mol, object_2.mol])
    product = product_tuples[0][0]
    product_smiles = Chem.MolToSmiles(product)

            
    if synthon_position == 1:
        ids= f"{reaction_id}-{object_1.ID}_s2_0"
    else:
        ids= f"{reaction_id}-s1_{object_2.ID}_0"

    return MoleculeInfo(role=f"product_{synthon_position}", ID=ids, smiles=product_smiles)


def enumerator_3_component(object_1, object_2, object_3, rxn_mol, synthon_position, reaction_id):
    """
    Enumerates product from 3 input MoleculeInfo objects using RDKit reaction and returns a MoleculeInfo for the product.
    """
    # Apply the reaction
    product_tuples = rxn_mol.RunReactants([object_1.mol, object_2.mol, object_3.mol])
    product = product_tuples[0][0]
    product_smiles = Chem.MolToSmiles(product)

    # Determine which synthon to use for the ID
    if synthon_position == 1:
        synthon_id = object_1.ID
        product_id = f"{reaction_id}-{synthon_id}_s2_s3"
    elif synthon_position == 2:
        synthon_id = object_2.ID
        product_id = f"{reaction_id}-s1_{synthon_id}_s3"
    else:  # synthon_position == 3
        synthon_id = object_3.ID
        product_id = f"{reaction_id}-s1_s2_{synthon_id}"

    return MoleculeInfo(role=f"product_{synthon_position}", ID=product_id, smiles=product_smiles)


def sub_enumeration(current_syntons, current_minimal_cap, synthon_position, reaction_id):

        def replace_labels(smiles):
            SYNTHON_LABELS = ['U', 'Np', 'Pu', 'Am']
            if smiles is not None:
                for label in SYNTHON_LABELS:
                    smiles = smiles.replace(label, 'U')
            return smiles

        ultimate_reaction_smarts = '[U]-[*:1].[U]-[*:2]>>[*:1]-[*:2]'
        rxn_mol = AllChem.ReactionFromSmarts(ultimate_reaction_smarts)

    
        reag_1 = Chem.MolFromSmiles(replace_labels(current_syntons.smiles))
        reag_2 = Chem.MolFromSmiles(replace_labels(current_minimal_cap.smiles))


        product_tuples = rxn_mol.RunReactants((reag_1, reag_2))
        product = product_tuples[0][0]
        product_smiles = Chem.MolToSmiles(product)


        return MoleculeInfo(role=f"product_{synthon_position}", 
                ID=generate_mel_synthon_id(
                    reaction_id=reaction_id,
                    synthon_position=synthon_position,
                    x= current_syntons.ID), 
                smiles=product_smiles)


def generate_mel_synthon_id(reaction_id: str, synthon_position: int, x: str) -> str:
    """
    Generate MEL synthon ID based on reaction ID, synthon position, and x value.

    Args:
        reaction_id (str): The reaction ID.
        synthon_position (int): The position of the synthon (1, 2, or 3).
        x (str): The dynamic part of the ID.

    Returns:
        str: The generated MEL synthon ID.
    """
    if synthon_position == 1:
        return f"{reaction_id}-{x}_s2_s3"
    elif synthon_position == 2:
        return f"{reaction_id}-s1_{x}_s3"
    elif synthon_position == 3:
        return f"{reaction_id}-s1_s2_{x}"
    else:
        raise ValueError(f"Unsupported synthon position: {synthon_position}")
