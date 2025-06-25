# === Standard Library Imports ===
from rdkit import Chem
from rdkit.Chem import AllChem

import pandas as pd
import random
from typing import Optional, Literal

from app.functions.utills import extract_single_value
from app.config import Config
from app.functions.Molecules import get_random_molecule_info, generate_product

import os

def get_minimal_caps(reaction_id: str, mode: Literal["generation", "enumeration"] = "generation"):
    """
    Load minimal cap molecules for a given reaction ID and mode.

    Args:
        reaction_id (str): The ID of the reaction.
        mode (str): One of ["generation", "enumeration"]. Determines subfolder to load caps from.

    Returns:
        mol_infos_dict (dict): Dictionary of MoleculeInfo objects.
        mel_type (str): The mel_type associated with the reaction.
    """

    data_base_path = Config.BASE_DATA_PATH
    if not data_base_path:
        raise EnvironmentError("Config.BASE_DATA_PATH is not set.")

    # Load reaction metadata
    reaction_df = pd.read_csv(os.path.join(data_base_path, 'REACTION_file_for_MEL.tsv'), sep='\t')
    reaction_row = reaction_df[reaction_df['reaction_id'] == reaction_id]
    if reaction_row.empty:
        raise ValueError(f"Reaction ID {reaction_id} not found in REACTION_file_for_MEL.tsv.")

    rxn_smarts = extract_single_value(reaction_row, "Reaction")
    mel_type = extract_single_value(reaction_row, 'mel_type')

    # Determine file path based on mode
    mode_folder = {
        "generation": "Minimal_Caps_by_reaction_id_for_generation_MEL",
        "enumeration": "Minimal_Caps_by_reaction_id_for_enumeration_MEL"
    }.get(mode)

    if not mode_folder:
        raise ValueError("Mode must be 'generation' or 'enumeration'.")

    min_cap_path = os.path.join(data_base_path, mode_folder, f'{reaction_id}.txt')
    if not os.path.isfile(min_cap_path):
        raise FileNotFoundError(f"Minimal cap file not found: {min_cap_path}")

    min_cap_df = pd.read_csv(min_cap_path, sep='\t')
    mol_infos_dict = {
        f"min_cap_s{i}": get_random_molecule_info(min_cap_df, i, role=f"min_cap_s{i}")
        for i in [1, 2, 3]
    }

    return mol_infos_dict, mel_type


def get_random_synthons_for_reaction_id(reaction_id: str):

    data_base_path = Config.BASE_DATA_PATH
    if not data_base_path:
        raise EnvironmentError("Config.MEL_DATA_PATH is not set.")

    # Reaction Dataframe
    reaction_df = pd.read_csv(f'{data_base_path}/REACTION_file_for_MEL.tsv', sep='\t')
    reaction_row = reaction_df[reaction_df['reaction_id'] == reaction_id]

    if reaction_row.empty:
        raise ValueError(f"Reaction ID {reaction_id} not found in reaction_df.")

    rxn_smarts = extract_single_value(reaction_row, "Reaction")
    mel_type = extract_single_value(reaction_row, 'mel_type')

    # Synthons DataFrame
    synthon_path = f'{data_base_path}/SYNTHONS_by_reaction_id_modified_for_MEL/{reaction_id}.txt'
    synthon_df = pd.read_csv(synthon_path, sep='\t')

    s1= get_random_molecule_info(synthon_df, 1, role=f"synthon_1") 
    s2 = get_random_molecule_info(synthon_df, 2, role=f"synthon_2") 
    s3 = get_random_molecule_info(synthon_df, 3, role=f"synthon_3") 
    product = generate_product(rxn_smarts, [s1, s2, s3], reaction_id)

   # Assemble MoleculeInfo dictionary
    mol_infos_dict = {
        "synthon_1": s1,
        "synthon_2": s2,
        "synthon_3": s3,
        "product": product,
    }
    return mol_infos_dict,  mel_type

