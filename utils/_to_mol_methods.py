
import os
import libchebipy
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from dGbyG.config import kegg_additions_csv_path, kegg_compound_data_path, metanetx_database_path, hmdb_database_path, recon3d_mol_dir_path, chebi_database_path, lipidmaps_database_path

cache = {}


def inchi_to_mol(inchi:str, sanitize=True):
    mol = Chem.MolFromInchi(inchi, removeHs=False, sanitize=sanitize)
    return mol

def smiles_to_mol(smiles:str, sanitize=True):
    mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    return mol

def file_to_mol(path:str, sanitize=True):
    mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=sanitize)
    return mol

def kegg_compound_to_mol(entry:str, sanitize=True):
    path = os.path.join(kegg_compound_data_path, entry+'.mol')
    if os.path.exists(path):
        mol = file_to_mol(path)
    #elif download_kegg_compound(entry):
    #    mol = file_to_mol(path)
    else:
        kegg_additions_df = pd.read_csv(kegg_additions_csv_path, index_col=0)
        if entry in kegg_additions_df.index:
            inchi = kegg_additions_df.loc[entry, 'inchi']
            mol = inchi_to_mol(inchi) if pd.notna(inchi) else None
    return mol

def metanetx_id_to_mol(id:str, sanitize=True):
    if 'metanetx' not in cache:
        cache['metanetx'] = pd.read_csv(metanetx_database_path, sep='\t', header=351, index_col=0)
    metanetx_df = cache['metanetx']
    smiles = (metanetx_df.loc[id, 'SMILES'])
    mol = smiles_to_mol(smiles)
    return mol

def hmdb_id_to_mol(id:str, sanitize=True):
    if 'hmdb' not in cache:
        cache['hmdb'] = pd.read_csv(hmdb_database_path, index_col=0, dtype={33: object})
    hmdb_df = cache['hmdb']
    if len(id.replace('HMDB', '')) < 7:
        id = 'HMDB' + '0'*(7-len(id.replace('HMDB', ''))) + id.replace('HMDB', '')
    smiles = hmdb_df.loc[id, 'SMILES']
    mol = smiles_to_mol(smiles)
    return mol

def chebi_id_to_mol(id:str, sanitize=True):
    libchebipy.set_download_cache_path(chebi_database_path)
    chebi_entity = libchebipy.ChebiEntity(str(id), )
    smiles = chebi_entity.get_smiles()
    mol = smiles_to_mol(smiles)
    return mol

def lipidmaps_id_to_mol(id:str, sanitize=True):
    if 'lipidmaps' not in cache:
        cache['lipidmaps'] = pd.read_csv(lipidmaps_database_path, sep='\t', index_col=0)
    lipidmaps_df = cache['lipidmaps']
    smiles = lipidmaps_df.loc[id, 'smiles']
    mol = smiles_to_mol(smiles)
    return mol

def recon3d_id_to_mol(id:str, sanitize=True):
    if id+'.mol' in os.listdir(recon3d_mol_dir_path):
        path = os.path.join(recon3d_mol_dir_path, id+'.mol')
        mol = file_to_mol(path)
        return mol
    else:
        return None
    
def inchi_key_to_mol(inchi_key:str, sanitize=True):
    if 'hmdb' not in cache:
        cache['hmdb'] = pd.read_csv(hmdb_database_path, index_col=0, dtype={33: object})
    if 'metanetx' not in cache:
        cache['metanetx'] = pd.read_csv(metanetx_database_path, sep='\t', header=351, index_col=0)
    
    if smiles := cache['hmdb'].loc[cache['hmdb'].INCHI_KEY == inchi_key, 'SMILES'].to_list():
        return smiles_to_mol(smiles[0])
    elif smiles := cache['metanetx'].loc[cache['metanetx'].InChIKey == 'InChIKey='+inchi_key, 'SMILES'].to_list():
        return smiles_to_mol(smiles[0])
    elif smiles := cache['hmdb'].loc[cache['hmdb'].INCHI_KEY.apply(lambda x:x[:-2]) == inchi_key[:-2], 'SMILES'].to_list():
        return smiles_to_mol(smiles[0])
    elif smiles := cache['metanetx'].loc[cache['metanetx'].InChIKey.apply(lambda x:x[:-2] if pd.notna(x) else x) == 'InChIKey='+inchi_key[:-2], 'SMILES'].to_list():
        return smiles_to_mol(smiles[0])
    else:
        return None

def name_to_mol(name:str, sanitize=True):
    if 'name' not in cache:
        cache['name'] = pd.read_csv(metanetx_database_path, sep='\t', header=351, index_col=1)
    metanetx_df = cache['name']
    smiles = (metanetx_df.loc[name, 'SMILES'])
    mol = smiles_to_mol(smiles)
    return mol