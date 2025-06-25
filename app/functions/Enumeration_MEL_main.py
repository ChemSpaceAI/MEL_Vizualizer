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

# Import custom functions from your project
from app.functions.Chemistry_logic import (
    add_mols,
    add_charge,
)
from app.functions.Molecules import MoleculeInfo


from app.config import Config 
from app.functions.utills import extract_single_value
from app.functions.Molecules import get_random_molecule_info, generate_product, MoleculeInfo
from app.functions.Enumerators import enumerator_2_component, enumerator_3_component, sub_enumeration
from app.functions.Enumeration_MEL import *

def get_enumeration(reaction_id: str):

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
    s1, s2, s3 = [ get_random_molecule_info(synthon_df, i, role=f"synthon_s{i}")for i in [1, 2, 3]]

    min_cap_path = f'{data_base_path}/Minimal_Caps_by_reaction_id_for_generation_MEL/{reaction_id}.txt'
    min_cap_df = pd.read_csv(min_cap_path, sep='\t')
    min_mol_info_1, min_mol_info_2, min_mol_info_3 = [ get_random_molecule_info(min_cap_df, i, role=f"min_cap_s{i}")for i in [1, 2, 3]]

    min_cap_path_enum = f'{data_base_path}/Minimal_Caps_by_reaction_id_for_enumeration_MEL/{reaction_id}.txt'
    min_cap_df_enum = pd.read_csv(min_cap_path_enum, sep='\t')
    min_mol_enum_1, min_mol_enum_2, min_mol_enum_3 = [ get_random_molecule_info(min_cap_df_enum, i, role=f"min_cap_s{i}")for i in [1, 2, 3]]

    len_s1 = synthon_df[(synthon_df['reaction_id'] == reaction_id) & (synthon_df['synton#'] == 1)].shape[0]
    len_s2 = synthon_df[(synthon_df['reaction_id'] == reaction_id) & (synthon_df['synton#'] == 2)].shape[0]
    len_s3 = synthon_df[(synthon_df['reaction_id'] == reaction_id) & (synthon_df['synton#'] == 3)].shape[0]

    if mel_type == '2_component':

        product_1a = enumerator_2_component(s1, min_mol_info_2, 
                                                rxn_mol, synthon_position = 1, reaction_id=reaction_id)

        s1a, s2a, product_2a = enumerate_after_generation_2_comp(product_1a.ID, rxn_mol, synthon_df, reaction_id)

        product_1b = enumerator_2_component(min_mol_info_2,  s2, 
                                                rxn_mol, synthon_position = 2, reaction_id=reaction_id)
                                                
        s1b, s2b, product_2b = enumerate_after_generation_2_comp(product_1b.ID, rxn_mol, synthon_df, reaction_id)

        mol_infos_dict = {
            "query1a": product_1a,
            "reagent_1a": s1a,
            "reagent_2a": s2a,
            "product_2a": product_2a,

            "query1b": product_1b,
            "reagent_1b": s1b,
            "reagent_2b": s2b,
            "product_2b": product_2b,
        }
        
        tables_info = {"len_s1": f"{len_s1:,}",
                     "len_s2": f"{len_s2:,}",
                     "subspace": f"{len_s1 * len_s2:,}"}

        return mol_infos_dict, tables_info, mel_type


    elif mel_type == '3_component':
        result_dict = {}
        dict_of_mol_infos = {}

        synthon_combinations = {
            1: (s1, min_mol_info_2, min_mol_info_3),
            2: (min_mol_info_1, s2, min_mol_info_3),
            3: (min_mol_info_1, min_mol_info_2, s3),
        }

        for position, (syn1, syn2, syn3) in synthon_combinations.items():
            try:
                # Step 1: First generation
                product_query = enumerator_3_component(
                    syn1, syn2, syn3, rxn_mol,
                    synthon_position=position,
                    reaction_id=reaction_id
                )

                if not product_query:
                    raise ValueError("Empty product_query")

                # Step 2: First iteration - two branches
                result_1, result_2 = enumerate_after_generation_3_comp_first_iter(
                    product_query.ID, rxn_mol, synthon_df, reaction_id, min_cap_df_enum
                )
                s1a, s2a, s3a, product_iter_1a = result_1
                s1b, s2b, s3b, product_iter_1b = result_2

                # Step 3: Second iteration for both branches
                s1aa, s2aa, s3aa, product_iter_2a = enumerate_after_generation_3_comp_second_iter(
                    product_iter_1a.ID, rxn_mol, synthon_df, reaction_id
                )
                s1bb, s2bb, s3bb, product_iter_2b = enumerate_after_generation_3_comp_second_iter(
                    product_iter_1b.ID, rxn_mol, synthon_df, reaction_id
                )

            except Exception as e:
                print(e)
                # Use empty MoleculeInfo fallback on error
                product_query = product_iter_1a = product_iter_1b = MoleculeInfo()
                product_iter_2a = product_iter_2b = MoleculeInfo()
                s1a = s2a = s3a = s1b = s2b = s3b = MoleculeInfo()
                s1aa = s2aa = s3aa = s1bb = s2bb = s3bb = MoleculeInfo()

            # Store everything in the dictionary
            dict_of_mol_infos[f"case_{position}"] = {
                "query_1": product_query,

                "reagent_1a": s1a,
                "reagent_2a": s2a,
                "reagent_3a": s3a,
                "product_1a": product_iter_1a,

                "reagent_1b": s1b,
                "reagent_2b": s2b,
                "reagent_3b": s3b,
                "product_1b": product_iter_1b,

                "query_2a": product_iter_1a,
                "reagent_4a": s1aa,
                "reagent_5a": s2aa,
                "reagent_6a": s3aa,
                "product_2a": product_iter_2a,

                "query_2b": product_iter_1b,
                "reagent_4b": s1bb,
                "reagent_5b": s2bb,
                "reagent_6b": s3bb,
                "product_2b": product_iter_2b,
            }


        tables_info = {"len_s1": f"{len_s1:,}",
                    "len_s2": f"{len_s2:,}",
                    "len_s3": f"{len_s3:,}",
                    "subspace": f"{len_s1 * len_s2 * len_s3:,}"}


        return dict_of_mol_infos,  mel_type,tables_info

    elif mel_type == 'bridge':

        bridge_position = int(reaction_row.iloc[0]['bridge_position'])
        dict_of_mol_infos = {}

        if bridge_position == 1:
            # Subenumerate: min_s2 + s2
            
            mel_fragment_1a = sub_enumeration(
                current_syntons=s2,  # synthon object
                current_minimal_cap=min_mol_info_2,  # mincap object
                synthon_position=2,
                reaction_id=reaction_id,

            )
            mel_fragment_1b = sub_enumeration(
                current_syntons=s3,
                current_minimal_cap=min_mol_info_3,
                synthon_position=3,
                reaction_id=reaction_id,
            )
            
            dict_of_mol_infos = {
                "bridge": s1,
                "subbridge_a": min_mol_info_2,
                "subbridge_b": min_mol_info_3,
            }

        elif bridge_position == 2:
            # Subenumerate: min_s1 + s1
            mel_fragment_1a= sub_enumeration(
                current_syntons=s1,
                current_minimal_cap=min_mol_info_1,
                synthon_position=1,
                reaction_id=reaction_id,

            )
            # Subenumerate: min_s3 + s3
            mel_fragment_1b = sub_enumeration(
                current_syntons=s3,
                current_minimal_cap=min_mol_info_3,
                synthon_position=3,
                reaction_id=reaction_id,

            )
            dict_of_mol_infos = {
                "bridge": s2,
                "subbridge_a": min_mol_info_1,
                "subbridge_b": min_mol_info_3,
            }


        elif bridge_position == 3:
            # Subenumerate: min_s1 + s1
            mel_fragment_1a = sub_enumeration(
                current_syntons=s1,
                current_minimal_cap=min_mol_info_1,
                synthon_position=1,
                reaction_id=reaction_id,

            )
            # Subenumerate: min_s2 + s2
            mel_fragment_1b = sub_enumeration(
                current_syntons=s2,
                current_minimal_cap=min_mol_info_2,
                synthon_position=2,
                reaction_id=reaction_id,

            )
            dict_of_mol_infos = {
                "bridge": s3,
                "subbridge_a": min_mol_info_1,
                "subbridge_b": min_mol_info_2,
            }

        s1a, s2a, s3a, product_iter_1a = enumerate_after_generation_bridge_first_iter(
                mel_fragment_1a.ID, rxn_mol, synthon_df, reaction_id, min_cap_df_enum, bridge_position
        )
        s1aa, s2aa, s3aa, product_iter_2a = enumerate_after_generation_3_comp_second_iter(
                product_iter_1a.ID, rxn_mol, synthon_df, reaction_id
            )


        s1b, s2b, s3b, product_iter_1b = enumerate_after_generation_bridge_first_iter(
                mel_fragment_1b.ID, rxn_mol, synthon_df, reaction_id, min_cap_df_enum, bridge_position
        )
        s1bb, s2bb, s3bb, product_iter_2b = enumerate_after_generation_3_comp_second_iter(
                product_iter_1b.ID, rxn_mol, synthon_df, reaction_id
            )

        if bridge_position == 1:
            mel_fragment_1a.update_role(bridge_position)
            mel_fragment_1b.update_role(bridge_position)
        elif bridge_position == 2:
            mel_fragment_1a.update_role_with_bridge_position(bridge_position)
            mel_fragment_1b.update_role_with_bridge_position(bridge_position)
        elif bridge_position == 3:
            mel_fragment_1a.update_role_with_bridge_position(bridge_position)
            mel_fragment_1b.update_role_with_bridge_position(bridge_position)

        dict_of_mol_infos['case_1'] = {
            "query_1a": mel_fragment_1a,
            "reagent_1a": s1a,
            "reagent_2a": s2a,
            "reagent_3a": s3a,
            "product_1a": product_iter_1a,

            "query_2a": product_iter_1a,
            "reagent_4a": s1aa,
            "reagent_5a": s2aa,
            "reagent_6a": s3aa,
            "product_2a": product_iter_2a,

        }
        dict_of_mol_infos['case_2'] = {
            "query_1a": mel_fragment_1a,
            "reagent_1a": s1b,
            "reagent_2a": s2b,
            "reagent_3a": s3b,
            "product_1a": product_iter_1b,

            "query_2a": product_iter_1b,
            "reagent_4a": s1bb,
            "reagent_5a": s2bb,
            "reagent_6a": s3bb,
            "product_2a": product_iter_2b,

        }

        tables_info = {"len_s1": f"{len_s1:,}",
                    "len_s2": f"{len_s2:,}",
                    "len_s3": f"{len_s3:,}",
                    "subspace": f"{len_s1 * len_s2 * len_s3:,}",
                    "bridge_position": bridge_position}

        return dict_of_mol_infos, mel_type, tables_info

    return {},  mel_type, {}
