 # === Standard Library Imports ===
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, AllChem, Draw
import pandas as pd
import re
import json
import io
import random
import os
from typing import Optional

from app.config import Config 
from app.functions.utills import extract_single_value
from app.functions.Molecules import get_random_molecule_info, generate_product, MoleculeInfo
from app.functions.Enumerators import enumerator_2_component, enumerator_3_component, sub_enumeration

def get_generation(reaction_id: str):

    data_base_path = Config.BASE_DATA_PATH
    if not data_base_path:
        raise EnvironmentError("Config.BASE_DATA_PATH is not set.")

    # Load and validate reaction row
    reaction_df = pd.read_csv(f'{data_base_path}/REACTION_file_for_MEL.tsv', sep='\t')
    reaction_row = reaction_df[reaction_df['reaction_id'] == reaction_id]
    if reaction_row.empty:
        raise ValueError(f"Reaction ID {reaction_id} not found in reaction_df.")

    # Extract metadata
    rxn_smarts = extract_single_value(reaction_row, "Reaction")
    rxn_mol = AllChem.ReactionFromSmarts(rxn_smarts)
    mel_type = extract_single_value(reaction_row, 'mel_type')

    # File paths
    synthon_path = f'{data_base_path}/SYNTHONS_by_reaction_id_modified_for_MEL/{reaction_id}.txt'
    synthon_df = pd.read_csv(synthon_path, sep='\t')
    min_cap_path = f'{data_base_path}/Minimal_Caps_by_reaction_id_for_generation_MEL/{reaction_id}.txt'
    min_cap_df = pd.read_csv(min_cap_path, sep='\t')

    # Usage with MoleculeInfo objects
    s1= get_random_molecule_info(synthon_df, 1, role=f"synthon_1") 
    s2 = get_random_molecule_info(synthon_df, 2, role=f"synthon_2") 
    s3 = get_random_molecule_info(synthon_df, 3, role=f"synthon_3") 
    min_mol_info_1, min_mol_info_2, min_mol_info_3 = [ get_random_molecule_info(min_cap_df, i, role=f"min_cap_s{i}")for i in [1, 2, 3]]

  
    if mel_type == '2_component':

        product_1 = enumerator_2_component(s1, min_mol_info_2, 
                                                rxn_mol, synthon_position = 1, reaction_id=reaction_id)
        product_2 = enumerator_2_component(min_mol_info_2,  s2, 
                                                rxn_mol, synthon_position = 2, reaction_id=reaction_id)
 
        mol_infos_dict = {
            "reagent_1": s1,
            "reagent_2":min_mol_info_2,
            "product_1": product_1,

            "reagent_3": min_mol_info_1,
            "reagent_4": s2,
            "product_2": product_2,

        }
        for key , item in mol_infos_dict.items() :

            print(key, ':',item.smiles, item.ID)

        return mol_infos_dict,  mel_type

    elif mel_type == '3_component':
        result_dict = {}

        result_dict[1] = enumerator_3_component(
            s1,
            min_mol_info_2,
            min_mol_info_3,
            rxn_mol,
            synthon_position=1, reaction_id=reaction_id
        )
        result_dict[2] = enumerator_3_component(
            min_mol_info_1,
            s2,
            min_mol_info_3,
            rxn_mol,
            synthon_position=2,
            reaction_id = reaction_id
        )
        result_dict[3] = enumerator_3_component(
            min_mol_info_1,
            min_mol_info_2,
            s3,
            rxn_mol,
            synthon_position=3,
            reaction_id=reaction_id
        )


        # Build visual mol_infos_dict
        mol_infos_dict = {
            "reagent_1": s1,
            "reagent_2": min_mol_info_2,
            "reagent_3": min_mol_info_3,
            "product_1": result_dict[1],  # first molecule from route 1

            "reagent_4": min_mol_info_1,
            "reagent_5": s2,
            "reagent_6": min_mol_info_3,
            "product_2": result_dict[2],  # route 2

            "reagent_7": min_mol_info_1,
            "reagent_8": min_mol_info_2,
            "reagent_9": s3,
            "product_3": result_dict[3],  # route 3
        }
        return mol_infos_dict, mel_type

    elif mel_type == 'bridge':

        bridge_position = int(reaction_row.iloc[0]['bridge_position'])


        if bridge_position == 1:
            # Subenumerate: min_s2 + s2
            product_1 = sub_enumeration(
                current_syntons=s2,  # synthon object
                current_minimal_cap=min_mol_info_2,  # mincap object
                synthon_position=2,
                reaction_id=reaction_id,

            )
            # Subenumerate: min_s3 + s3
            product_2 = sub_enumeration(
                current_syntons=s3,
                current_minimal_cap=min_mol_info_3,
                synthon_position=3,
                reaction_id=reaction_id,

            )

            # Bridge is between S2 and S3
            mol_infos_dict = {
                "reagent_1": s2,
                "reagent_2": min_mol_info_2,
                "product_1": product_1,  # min_s2 + s2

                "reagent_3": s3,
                "reagent_4": min_mol_info_3,
                "product_2": product_2,  # min_s3 + s3
            }

        elif bridge_position == 2:
            # Subenumerate: min_s1 + s1
            product_1= sub_enumeration(
                current_syntons=s1,
                current_minimal_cap=min_mol_info_1,
                synthon_position=1,
                reaction_id=reaction_id,

            )
            # Subenumerate: min_s3 + s3
            product_2 = sub_enumeration(
                current_syntons=s3,
                current_minimal_cap=min_mol_info_3,
                synthon_position=3,
                reaction_id=reaction_id,

            )
            # Bridge is between S1 and S3
            mol_infos_dict = {
                "reagent_1": s1,
                "reagent_2": min_mol_info_1,
                "product_1": product_1,

                "reagent_3": s3,
                "reagent_4": min_mol_info_3,
                "product_2": product_2,
            }


        elif bridge_position == 3:
            # Subenumerate: min_s1 + s1
            product_1 = sub_enumeration(
                current_syntons=s1,
                current_minimal_cap=min_mol_info_1,
                synthon_position=1,
                reaction_id=reaction_id,

            )
            # Subenumerate: min_s2 + s2
            product_2 = sub_enumeration(
                current_syntons=s2,
                current_minimal_cap=min_mol_info_2,
                synthon_position=2,
                reaction_id=reaction_id,

            )
        
            mol_infos_dict = {
                "reagent_1": s1,
                "reagent_2": min_mol_info_1,
                "product_1": product_1,

                "reagent_3": s2,
                "reagent_4": min_mol_info_2,
                "product_2": product_2,
            }

        return mol_infos_dict, mel_type
         
    return {},  mel_type

 

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
