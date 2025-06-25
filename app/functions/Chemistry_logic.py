import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

SYNTHON_LABELS = ['U','Np','Pu','Am']


def add_charge(mol, charge_tertiary=True, charge_imidazole=True):
    # Define SMARTS queries for different nitrogen-containing substructures
    Chem.SanitizeMol(mol)
    substructure_queries = {
        "amine": '[$([NX3]);!$(N~[!C]);!$(NC~[!#6]);!$(N[#6][#6]([F,Cl,Br])[F,Cl,Br]);!$(N(C=*)@C=*);!$(N[#6]=[#6][#6]=[*])]',

        "imidazole": "[$([nX2]);$([nr5]:[cr5]:[nr5]:[cr5]:[cr5])]",

        "guanidine": "[$([NX2]);$([#7]=[#6]([#7])[#7])]",

        "amidine": "[$([NX2]);$([#7]=[#6]([#7])[C])]",

        "tetrazole": "[$([nX3;H1]);$([nr5]:[nr5]:[cr5]:[nr5]:[nr5]),$([nr5]:[cr5]:[nr5]:[nr5]:[nr5])]"
    }

    # Loop through each subgroup and modify the charges for matching substructures
    for subgroup, smarts in substructure_queries.items():
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        
        # Apply charge to nitrogen atoms in the matching substructures
        for match in matches:
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0:  # Check if it's a nitrogen atom with no charge
                    if subgroup == "tetrazole":
                        atom.SetFormalCharge(-1)  # Set negative charge for tetrazole
                        atom.SetNumExplicitHs(0)  # Remove explicit hydrogen to maintain valence
                    else:
                        atom.SetFormalCharge(+1)  # Set positive charge for other subgroups

    return mol 

def add_mols(current_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a new column to the DataFrame with RDKit molecule objects generated from the SMILES column
    """
    synthon_df = current_df.copy()
    synthon_df['mol'] = synthon_df['SMILES'].apply(Chem.MolFromSmiles)
    return synthon_df
