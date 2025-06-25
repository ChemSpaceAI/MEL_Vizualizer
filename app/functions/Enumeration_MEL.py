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




def enumerate_after_generation_bridge_first_iter(synthons_formula_ID, rxn_mol, synthons_df, reaction_id, minimal_caps_df, bridge_position ):

    # Parse synthon IDs from the formula

    _, synthons_formula_ID = synthons_formula_ID.split('-')
    s1_ID, s2_ID, s3_ID = synthons_formula_ID.split('_')

    #s1_s2_26213476 s1_23795974_s3 24703580_s2_s3
    
    if ("s2" == s2_ID) and ('s3' == s3_ID):

        #Enumeration 1: synthon_1 (anchor)  + synthon_2 + synthon_3

        anchor_synthon_ID = int(s1_ID)

        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID] 

        if bridge_position == 2:

            variable_synthons = synthons_df[synthons_df['synton#'] == 2]
            cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 3]
            
            result = enumerator_engine_for_3_component(
                synthons_1= anchor_synthon_row,
                synthons_2= variable_synthons,
                synthons_3= cap_synthons,
                rxn_mol=rxn_mol,
                reaction_id=reaction_id,
                cap_position=3
            )
            s1, s2, s3, product = result
            s1.update_role('Cap')
            s2.update_role('Anchor')
            s3.update_role('Variable (Bridge)')
            result = s1, s2, s3, product

            
        elif bridge_position == 3: 

            variable_synthons = synthons_df[synthons_df['synton#'] == 3]
            cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 2]
            
            result = enumerator_engine_for_3_component(
                synthons_1= anchor_synthon_row,
                synthons_2= cap_synthons,
                synthons_3= variable_synthons,
                rxn_mol=rxn_mol,
                reaction_id=reaction_id,
                cap_position=2
            )
            s1, s2, s3, product = result
            s1.update_role('Anchor')
            s2.update_role('Cap')
            s3.update_role('Variable (Bridge)')
            result = s1, s2, s3, product

    elif ("s1" == s1_ID) and ('s2' == s2_ID):

        #Enumeration 2:  synthon_1 + synthon_2 + synthon_3 (anchor)
        
        anchor_synthon_ID = int(s3_ID)

        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID] 

        if bridge_position == 1:

            variable_synthons = synthons_df[synthons_df['synton#'] == 1]
            cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 2]
            
            result = enumerator_engine_for_3_component(
                synthons_1= variable_synthons,
                synthons_2= cap_synthons,
                synthons_3= anchor_synthon_row,
                rxn_mol=rxn_mol,
                reaction_id=reaction_id,
                cap_position=2
            )

            s1, s2, s3, product = result
            s1.update_role('Variable (Bridge)')
            s2.update_role('Cap')
            s3.update_role('Anchor')
            result = s1, s2, s3, product

        elif bridge_position == 2:

            variable_synthons = synthons_df[synthons_df['synton#'] == 2]
            cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 1]
            
            result = enumerator_engine_for_3_component(
                synthons_1= cap_synthons,
                synthons_2= variable_synthons,
                synthons_3= anchor_synthon_row,
                rxn_mol=rxn_mol,
                reaction_id=reaction_id,
                cap_position=1
            )

            s1, s2, s3, product = result
            result =  s1.update_role('Cap'), s2.update_role('Variable (Bridge)'), s3.update_role('Anchor'), product

    elif ("s1" == s1_ID) and ('s3' == s3_ID):
        #Enumeration 3: synthon_1 + synthon_2 (anchor) + synthon_3 

        anchor_synthon_ID = int(s2_ID)
        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID] 

        if bridge_position == 1:

            variable_synthons = synthons_df[synthons_df['synton#'] == 1]
            cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 3]
            
            result = enumerator_engine_for_3_component(
                synthons_1= variable_synthons,
                synthons_2= anchor_synthon_row,
                synthons_3= cap_synthons,
                rxn_mol=rxn_mol,
                reaction_id=reaction_id,
                cap_position=3
            )

            s1, s2, s3, product = result
            s1.update_role('Variable (Bridge)')
            s2.update_role('Anchor')
            s3.update_role('Cap')
            result = s1, s2, s3, product

        elif bridge_position == 3:

            variable_synthons = synthons_df[synthons_df['synton#'] == 3]
            cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 1]
            
            result = enumerator_engine_for_3_component(
                synthons_1= cap_synthons,
                synthons_2= anchor_synthon_row,
                synthons_3= variable_synthons,
                rxn_mol=rxn_mol,
                reaction_id=reaction_id, 
                cap_position=1
            )
            s1, s2, s3, product = result
            s1.update_role('Cap')
            s2.update_role('Anchor')
            s3.update_role('Variable (Bridge)')
            result = s1, s2, s3, product

    # Handle cases where the synthon IDs are not as expected
    else:
        logger.warning(f"Unexpected synthon IDs format: {mel_ID}")
        return None  # Return None for unhandled synthon formats
    return result


def enumerate_after_generation_2_comp(synthons_formula_ID, rxn_mol, synthons_df, reaction_id ):

    # Parse synthon IDs from the formula
    _, synthons_formula_ID = synthons_formula_ID.split('-')
    s1_ID, s2_ID, s3_ID = synthons_formula_ID.split('_')

    # If it's a 3-component reaction, we don't handle it in this function
    if s3_ID != '0':
        logger.warning(f"Unexpected 3rd synthon ID ({s3_ID}) for {mel_ID}. Skipping.")
        return None  # Skip 3-component reactions

    # Case where s1 is the variable synthon
    if s1_ID == "s1":
        anchor_synthon_ID = int(s2_ID)
        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID]

        if anchor_synthon_row.empty:
            logger.warning(f"Anchor synthon ID {anchor_synthon_ID} not found. Skipping {mel_ID}.")
            return None  # Skip if anchor synthon not found

        variable_synthons = synthons_df[synthons_df['synton#'] == 1].sample(1)

        result = enumerator_engine_for_2_component(
            synthons_1=variable_synthons,
            synthons_2=anchor_synthon_row,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id,
        )
        s1, s2,  product = result
        s1.update_role('Variable')
        s2.update_role('Anchor')
        result = s1, s2, product
       
    # Case where s2 is the variable synthon
    elif s2_ID == "s2":
        anchor_synthon_ID = int(s1_ID)
        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID]

        if anchor_synthon_row.empty:
            logger.warning(f"Anchor synthon ID {anchor_synthon_ID} not found. Skipping {mel_ID}.")
            return None  # Skip if anchor synthon not found

        variable_synthons = synthons_df[synthons_df['synton#'] == 2].sample(1)

        result = enumerator_engine_for_2_component(
            synthons_1=anchor_synthon_row,
            synthons_2=variable_synthons,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id,
        )
        s1, s2,  product = result
        s1.update_role('Anchor')
        s2.update_role('Variable')
        result = s1, s2, product
      
    # Handle cases where the synthon IDs are not as expected
    else:
        logger.warning(f"Unexpected synthon IDs format: {mel_ID}")
        return None  # Return None for unhandled synthon formats

    return result  # Return the generated molecules as a dictionary


def enumerate_after_generation_3_comp_first_iter(synthons_formula_ID, rxn_mol, synthons_df, reaction_id, minimal_caps_df ):


    # Parse synthon IDs from the formula
    _, synthons_formula_ID = synthons_formula_ID.split('-')
    s1_ID, s2_ID, s3_ID = synthons_formula_ID.split('_')

    if ("s2" == s2_ID) and ('s3' == s3_ID):

        #Enumeration 1: synthon_1 (anchor)  + synthon_2 + synthon_3
        anchor_synthon_ID = int(s1_ID)
        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID] 

        #Combinatorial Case 1:  synthon_1 (anchor) + synthon_2 (variable_full) + synthon_3 (cap)
        variable_synthons = synthons_df[synthons_df['synton#'] == 2]
        cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 3]

        result1 = enumerator_engine_for_3_component(
            synthons_1= anchor_synthon_row,
            synthons_2= variable_synthons,
            synthons_3= cap_synthons,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id, 
            cap_position=3
        )
        s1, s2, s3, product = result1
        s1.update_role('Anchor')
        s2.update_role('Variable')
        s3.update_role('Cap')
        result1 = s1, s2, s3, product

        #Combinatorial Case 2:  synthon_1 (anchor) + synthon_2 (cap) + synthon_3 (variable_full) 
        variable_synthons = synthons_df[synthons_df['synton#'] == 3]
        cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 2]

        result2 = enumerator_engine_for_3_component(
            synthons_1= anchor_synthon_row,
            synthons_2= cap_synthons,
            synthons_3= variable_synthons,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id, 
            cap_position = 2
        )

        s1, s2, s3, product = result2
        s1.update_role('Anchor')
        s2.update_role('Cap')
        s3.update_role('Variable')
        result2 = s1, s2, s3, product

    elif ("s1" == s1_ID) and ('s2' == s2_ID):

        #Enumeration 2:  synthon_1 + synthon_2 + synthon_3 (anchor)
        anchor_synthon_ID = int(s3_ID)
        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID] 

        #Combinatorial Case 1:  synthon_1 (cap) + synthon_2 (variable_full) + synthon_3 (anchor) 
        variable_synthons = synthons_df[synthons_df['synton#'] == 2]
        cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 1]

        result1 = enumerator_engine_for_3_component(
            synthons_1= cap_synthons, 
            synthons_2= variable_synthons,
            synthons_3= anchor_synthon_row,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id, 
            cap_position= 1
        )
        s1, s2, s3, product = result1
        s1.update_role('Cap')
        s2.update_role('Variable')
        s3.update_role('Anchor')
        result1 = s1, s2, s3, product

        #Combinatorial Case 2:  synthon_1 (variable_full)  + synthon_2 (cap) + synthon_3 (anchor) 
        variable_synthons = synthons_df[synthons_df['synton#'] == 1]
        cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 2]

        result2 = enumerator_engine_for_3_component(
            synthons_1= variable_synthons,
            synthons_2= cap_synthons, 
            synthons_3= anchor_synthon_row,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id, 
            cap_position = 2
        )
        s1, s2, s3, product = result2
        s1.update_role('Variable')
        s2.update_role('Cap')
        s3.update_role('Anchor')
        result2 = s1, s2, s3, product
        
    

    elif ("s1" == s1_ID) and ('s3' == s3_ID):

        #Enumeration 3: synthon_1 + synthon_2 (anchor) + synthon_3 
        anchor_synthon_ID = int(s2_ID)
        anchor_synthon_row = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID] 

        #Combinatorial Case 1:  synthon_1 (cap) + synthon_2 (anchor) + synthon_3 (variable_full)
        variable_synthons = synthons_df[synthons_df['synton#'] == 3]
        cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 1]

        result1 = enumerator_engine_for_3_component(
            synthons_1= cap_synthons, 
            synthons_2= anchor_synthon_row, 
            synthons_3= variable_synthons,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id, 
            cap_position=1
        )
        s1, s2, s3, product = result1
        s1.update_role('Cap')
        s2.update_role('Anchor')
        s3.update_role('Variable')
        result1 = s1, s2, s3, product

        #Combinatorial Case 2:  synthon_1 (variable_full) + synthon_2 (anchor) + synthon_3 (cap)
        variable_synthons = synthons_df[synthons_df['synton#'] == 1]
        cap_synthons = minimal_caps_df[minimal_caps_df['synton#'] == 3]

        result2 = enumerator_engine_for_3_component(
            synthons_1= variable_synthons,
            synthons_2= anchor_synthon_row, 
            synthons_3= cap_synthons,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id, 
            cap_position=3
        )
        s1, s2, s3, product = result2
        s1.update_role('Variable')
        s2.update_role('Anchor')
        s3.update_role('Cap')
        result2 = s1, s2, s3, product

    else:
        logger.warning(f"Unexpected synthon IDs format: {mel_ID}")
        return {}
    
    return result1, result2 

def enumerate_after_generation_3_comp_second_iter(synthons_formula_ID, rxn_mol, synthons_df, reaction_id ):

    # Parse synthon IDs from the formula
    _, synthons_formula_ID = synthons_formula_ID.split('-')
    s1_ID, s2_ID, s3_ID = synthons_formula_ID.split('_')

    if s3_ID == 's3':

        #Enumeration 1: synthon_1 (anchor)  + synthon_2 (anchor) + synthon_3 (variable_full)

        anchor_synthon_ID_1 = int(s1_ID)
        anchor_synthon_row_1 = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID_1] 

        anchor_synthon_ID_2 = int(s2_ID)
        anchor_synthon_row_2 = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID_2] 

        variable_synthons = synthons_df[synthons_df['synton#'] == 3]

        result = enumerator_engine_for_3_component(
            synthons_1= anchor_synthon_row_1,
            synthons_2= anchor_synthon_row_2,
            synthons_3= variable_synthons,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id,
            cap_position=0
        )
        s1, s2, s3, product = result
        s1.update_role('Anchor')
        s2.update_role('Anchor')
        s3.update_role('Variable')
        result = s1, s2, s3, product

    elif s2_ID == 's2':

        #Enumeration 1: synthon_1 (anchor)  + synthon_2 (variable_full)  + synthon_3  (anchor)

        anchor_synthon_ID_1 = int(s1_ID)
        anchor_synthon_row_1 = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID_1] 

        anchor_synthon_ID_2 = int(s3_ID)
        anchor_synthon_row_2 = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID_2] 

        variable_synthons = synthons_df[synthons_df['synton#'] == 2]
        result = enumerator_engine_for_3_component(
            synthons_1= anchor_synthon_row_1,
            synthons_2= variable_synthons,
            synthons_3= anchor_synthon_row_2,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id,
            cap_position=0
        )
        s1, s2, s3, product = result
        s1.update_role('Anchor')
        s2.update_role('Variable')
        s3.update_role('Anchor')
        result = s1, s2, s3, product

    elif s1_ID == 's1':

        #Enumeration 1: synthon_1 (variable_full)  + synthon_2  (anchor)  + synthon_3  (anchor)

        anchor_synthon_ID_1 = int(s2_ID)
        anchor_synthon_row_1 = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID_1] 

        anchor_synthon_ID_2 = int(s3_ID)
        anchor_synthon_row_2 = synthons_df[synthons_df['synton_id'] == anchor_synthon_ID_2] 

        variable_synthons = synthons_df[synthons_df['synton#'] == 1]

        result = enumerator_engine_for_3_component(
            synthons_1= variable_synthons,
            synthons_2= anchor_synthon_row_1,
            synthons_3= anchor_synthon_row_2,
            rxn_mol=rxn_mol,
            reaction_id=reaction_id,
            cap_position=0
        )
        s1, s2, s3, product = result
        s1.update_role('Variable')
        s2.update_role('Anchor')
        s3.update_role('Anchor')
        result = s1, s2, s3, product

    # Handle cases where the synthon IDs are not as expected
    else:
        logger.warning(f"Unexpected synthon IDs format: {mel_ID}")
        return None  # Return None for unhandled synthon formats

    return result

def enumerator_engine_for_3_component(
    synthons_1: pd.DataFrame, 
    synthons_2: pd.DataFrame, 
    synthons_3: pd.DataFrame, 
    rxn_mol,
    reaction_id: str, 
    cap_position: int 
):
    # Sample and prepare molecules
    synthons_1 = add_mols(synthons_1.sample(1))
    synthons_2 = add_mols(synthons_2.sample(1))
    synthons_3 = add_mols(synthons_3.sample(1))

    for i, row1 in synthons_1.iterrows():
        for j, row2 in synthons_2.iterrows():
            for k, row3 in synthons_3.iterrows():
                try:
                    # Run the reaction
                    product_tuples = rxn_mol.RunReactants((row1["mol"], row2["mol"], row3["mol"]))
                    product = product_tuples[0][0]

                    product_copy = Chem.Mol(product)
                    Chem.SanitizeMol(product_copy)
                    charged_product = add_charge(product_copy)

                    result_smile = Chem.MolToSmiles(product)
                    charged_smile = Chem.MolToSmiles(charged_product)

                    # Synthon ID formatting
                    synthon_parts = [row1["synton_id"], row2["synton_id"], row3["synton_id"]]
                    if cap_position in {1, 2, 3}:
                        synthon_parts[cap_position - 1] = f's{cap_position}'
                    result_synthon_id = f'{reaction_id}-{"_".join(map(str, synthon_parts))}'

                    # Build MoleculeInfo objects
                    synthon_1_info = MoleculeInfo(
                        ID=row1["synton_id"],
                        smiles=row1["SMILES"],
                        role="synthon_1"
                    )

                    synthon_2_info = MoleculeInfo(
                        ID=row2["synton_id"],
                        smiles=row2["SMILES"],
                        role="synthon_2"
                    )

                    synthon_3_info = MoleculeInfo(
                        ID=row3["synton_id"],
                        smiles=row3["SMILES"],
                        role="synthon_3"
                    )

                    product_info = MoleculeInfo(
                        ID=result_synthon_id,
                        smiles=result_smile,
                        role="product"
                    )

                    return synthon_1_info, synthon_2_info, synthon_3_info, product_info

                except Exception as e:
                    try:
                        rxn_smarts = AllChem.ReactionToSmarts(rxn_mol)
                    except Exception:
                        rxn_smarts = "[could not convert rxn_mol to SMARTS]"

                    print(f"[Error] Reaction failed for: {row1['SMILES']} + {row2['SMILES']} + {row3['SMILES']} | Reaction_ID: {reaction_id} | Error: {e} | RXN_SMARTS: {rxn_smarts}")

    # If no product could be generated, return None or dummy MoleculeInfos
    return MoleculeInfo(role="synthon_1"), MoleculeInfo(role="synthon_2"), MoleculeInfo(role="synthon_3"), MoleculeInfo(role="product")


def enumerator_engine_for_2_component(
    synthons_1: pd.DataFrame, 
    synthons_2: pd.DataFrame, 
    rxn_mol,
    reaction_id
    ):
    # Add RDKit mol objects to synthons
    synthons_1 = add_mols(synthons_1)
    synthons_2 = add_mols(synthons_2)

    # Iterate over all pairs, return the first successful result
    for _, row1 in synthons_1.iterrows():
        for _, row2 in synthons_2.iterrows():
            try:
                # Apply reaction
                product_tuples = rxn_mol.RunReactants([row1["mol"], row2["mol"]])
                if not product_tuples or not product_tuples[0]:
                    continue

                product = product_tuples[0][0]

                # Sanitize and charge the product
                product_copy = Chem.Mol(product)
                Chem.SanitizeMol(product_copy)
                charged_product = add_charge(product_copy)

                # Convert to SMILES
                result_smile = Chem.MolToSmiles(product)
                charged_smile = Chem.MolToSmiles(charged_product)

                # Build MoleculeInfo objects
                synthon_1_info = MoleculeInfo(
                    ID=row1["synton_id"],
                    smiles=row1["SMILES"],
                    role="synthon_1"
                )

                synthon_2_info = MoleculeInfo(
                    ID=row2["synton_id"],
                    smiles=row2["SMILES"],
                    role="synthon_2"
                )

                product_info = MoleculeInfo(
                    ID=f"{reaction_id}-{row1['synton_id']}_{row2['synton_id']}",
                    smiles=result_smile,
                    role="product"
                )

                return synthon_1_info, synthon_2_info, product_info

            except Exception as e:
                print(f"Error processing pair: {row1['SMILES']} + {row2['SMILES']} for Reaction_ID {reaction_id}. Exception: {str(e)}")
                continue

    # Fallback if nothing worked
    return MoleculeInfo(role="synthon_1"), MoleculeInfo(role="synthon_2"), MoleculeInfo(role="product")

