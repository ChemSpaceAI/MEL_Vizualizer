# === Standard Library Imports ===
from rdkit import Chem
import pandas as pd
import os
from app.config import Config
from pathlib import Path

def get_reaction_ids(mel_type: str, scope: str = 'unique') -> list[str]:

    df = pd.read_csv(Config.IDS_FILE_PATH, sep='\t')
    df = df[df['mel_type'] == mel_type]
    if scope == 'unique':
        df = df.drop_duplicates(subset='reaction_id_type')
    return df['reaction_id'].dropna().astype(str).unique().tolist()

def safe_mol_from_smiles(smiles):
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception as e:
        return Chem.MolFromSmiles('*')


def extract_single_value(df, column):
    """Extracts a single scalar value from a one-row DataFrame."""
    return df.iloc[0][column]


def cleanup_generated_images(directory, max_files=20):
    path = Path(directory)
    png_files = sorted(
        [f for f in path.glob("*.png") if f.is_file()],
        key=lambda f: f.stat().st_mtime  # Sort by modified time
    )

    if len(png_files) > max_files:
        for f in png_files[:-max_files]:  # Keep newest N, delete the rest
            f.unlink()
