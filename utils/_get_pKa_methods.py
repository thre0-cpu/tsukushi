
import os
import json
import jpype
import requests
import multiprocessing
from copy import deepcopy
import pandas as pd
from typing import List
from tqdm import tqdm

from dGbyG.utils.constants import *
from dGbyG.utils.CustomError import NoLicenseError
from dGbyG.config import chemaxon_jar_dir, chemaxon_pka_json_path, chemaxon_pka_csv_path


pka_cache = {}


def get_pka_from_chemaxon_rest(smiles, temperature):
    # 
    chemaxon_pka_api = 'https://jchem-microservices.chemaxon.com/jws-calculations/rest-v1/calculator/calculate/pka'
    headers = {'accept': '*/*', 'Content-Type': 'application/json'}
    pka_req_body = json.dumps({
        "inputFormat": "smiles",
        "micro": False,
        "outputFormat": "smiles",
        "outputStructureIncluded": False,
        "pKaLowerLimit": -20,
        "pKaUpperLimit": 20,
        "prefix": "STATIC",
        "structure": smiles,
        "temperature": temperature,
        "types": "acidic, basic",
        })
    try:
        pka = requests.post(chemaxon_pka_api, data=pka_req_body, headers=headers).json()
        return pka if not pka.get('error') else None
    except:
        return None



def _batch_get_pKa_from_chemaxon(smiles_list:List[str], temperature:float=default_T) -> List[str|dict]:
    jar_dir = chemaxon_jar_dir 
    fileList = [os.path.join(jar_dir,i)  for i in os.listdir(jar_dir)]

    jpype.startJVM('-Dchemaxon.license.url=/home/fanwc/.chemaxon/license.cxl')
    for p in fileList:
        jpype.addClassPath(p)

    MolImporter = jpype.JClass('chemaxon.formats.MolImporter')
    pKaPlugin = jpype.JClass('chemaxon.marvin.calculations.pKaPlugin')
    pKa = pKaPlugin()
    if not pKa.isLicensed():
        return ["ChemAxon license not found"]
    
    output = []
    pKa.setTemperature(temperature)
    pKa.setpKaPrefixType(2) # 'acidic,basic'
    pKa.setAcidicpKaUpperLimit(20)
    pKa.setBasicpKaLowerLimit(-20)

    if len(smiles_list) > 1:
        smiles_list = tqdm(smiles_list)
        
    for smiles in smiles_list:
        mol=MolImporter.importMol(smiles)
        pKa.setMolecule(mol)
        #print(smiles)
        if not pKa.run():
            output.append("pKa calculation failed")
            continue

        apka = []
        bpka = []
        for i in range(mol.getAtomCount()):
            apka.append({'atomIndex': i, 'value': float(pKa.getpKa(i, pKa.ACIDIC))})
            bpka.append({'atomIndex': i, 'value': float(pKa.getpKa(i, pKa.BASIC))})
        res = {'acidicValuesByAtom':apka, 'basicValuesByAtom':bpka}
        output.append(res)
        
    jpype.shutdownJVM()
    print('shut down JVM!')
    return output



def get_pKa_from_chemaxon(smiles:str, temperature:float=default_T) -> dict|None:
    if not isinstance(smiles, str):
        raise InputValueError("get_pKa_from_chemaxon(smiles:str, temperature:float=default_T), smiles must be a string")
    queue = multiprocessing.Queue()
    func = lambda queue, smiles, temperature: queue.put(_batch_get_pKa_from_chemaxon([smiles], temperature))
    p = multiprocessing.Process(target=func, args=(queue, smiles, temperature, ))
    p.start()
    p.join()

    res = queue.get()[0]
    if res == "ChemAxon license not found":
        raise NoLicenseError("ChemAxon license not found")
    elif res == "pKa calculation failed":
        print(smiles, res)
        return None
    elif isinstance(res, dict):
        return res
    else:
        raise Exception(f"Unknown error, return value: {res}")
    


def batch_calculation_pKa_to_json(smiles_list:list, temperature:float=default_T) -> None:
    # 
    if not isinstance(smiles_list, list):
        raise InputValueError("batch_calculation_pKa_to_json(), smiles must be a list")
    queue = multiprocessing.Manager().Queue()
    func = lambda queue, smiles_list, temperature: queue.put(_batch_get_pKa_from_chemaxon(smiles_list, temperature))
    p = multiprocessing.Process(target=func, args=(queue, smiles_list, temperature, ))
    p.start()
    p.join()
    print('ChemAxon pKa calculation finished')

    pKa_list = queue.get()
    if pKa_list[0] == "ChemAxon license not found":
        raise NoLicenseError("ChemAxon license not found")
    
    print('Reading json file...')
    with open(chemaxon_pka_json_path, "r", encoding="utf-8") as f:
        file = f.read()
        if len(file) > 0:
            data = json.loads(file)
        else:
            data = {}
    
    assert len(smiles_list) == len(pKa_list)
    for smiles, pka in zip(smiles_list, pKa_list):
        if pka == "pKa calculation failed":
            pka = None
        smiles_data = data.get(smiles, {})
        smiles_data[str(temperature)] = pka
        data[smiles] = smiles_data
    print('Writing json file...')
    with open(chemaxon_pka_json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, sort_keys=True, )
    print('Done')
    return None



def get_pka_from_json(smiles:str, temperature:float=default_T) -> dict|None:
    # 
    if 'chemaxon_pKa_json' not in pka_cache.keys():
        with open(chemaxon_pka_json_path, "r", encoding="utf-8") as f:
            pka_cache['chemaxon_pKa_json'] = json.load(f)
    
    pKa_json = pka_cache['chemaxon_pKa_json']
    if str(temperature) not in pKa_json.get(smiles, {}).keys():
        print(f'{smiles}, {temperature} K pKa not in file')
        batch_calculation_pKa_to_json([smiles], temperature)
        with open(chemaxon_pka_json_path, "r", encoding="utf-8") as f:
            pka_cache['chemaxon_pKa_json'] = json.load(f)
            pKa_json = pka_cache['chemaxon_pKa_json']
        
    return deepcopy(pKa_json[smiles][str(temperature)])




def get_pka_from_file(smiles, temperature):
    # 
    if 'chemaxon_csv' not in pka_cache.keys():
        pka_cache['chemaxon_csv'] = pd.read_csv(chemaxon_pka_csv_path, index_col=0)
    pKa_df = pka_cache['chemaxon_csv']

    if smiles not in pKa_df.index:
        print(smiles, 'pka not in file')
        return None

    pka_copy = {'acidicValuesByAtom':[], 'basicValuesByAtom': []}
    pka = pKa_df.loc[smiles, :]
    atoms_idx = 0
    for k, v in pka.items():
        if pd.notna(v) and k.startswith('apKa'):
            atomIndex = pka['atoms'].split(',')[atoms_idx]
            pka_copy['acidicValuesByAtom'].append({'atomIndex':int(atomIndex), 'value':float(v)})
            atoms_idx += 1
        elif pd.notna(v) and k.startswith('bpKa'):
            atomIndex = pka['atoms'].split(',')[atoms_idx]
            pka_copy['basicValuesByAtom'].append({'atomIndex':int(atomIndex), 'value':float(v)})
            atoms_idx += 1

    return pka_copy