import os, re
import numpy as np
import pandas as pd
from copy import deepcopy
from typing import Dict, Tuple, List

import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from dGbyG.config import chemaxon_pka_csv_path, chemaxon_pka_json_path
from dGbyG.utils.constants import *
from dGbyG.utils.CustomError import *
from ._to_mol_methods import *
from ._get_pKa_methods import get_pKa_from_chemaxon, get_pka_from_json, get_pka_from_file

pka_cache = {}


def parse_equation(equation:str, eq_sign=None) -> dict:
    # 
    eq_Signs = [' = ', ' <=> ', ' -> ']
    if eq_sign:
        equation = equation.split(eq_sign)
    else:
        for eq_sign in eq_Signs:
            if eq_sign in equation:
                equation = equation.split(eq_sign)
                break
            
        
        
    if not type(equation)==list:
        return {equation:1}
    
    equation_dict = {}
    equation = [x.split(' + ') for x in equation]

    for side, coefficient in zip(equation, (-1,1)):
        for t in side:
            if (reactant := re.match(r'^(\d+) (.+)$', t)) or (reactant := re.match(r'^(\d+\.\d+) (.+)$', t)):
                value = float(reactant.group(1))
                entry = reactant.group(2)
                equation_dict[entry] = equation_dict.get(entry,0) + value * coefficient
            else:
                equation_dict[t] = equation_dict.get(t,0) + coefficient
                
    if '' in equation_dict:
        equation_dict.pop('')

    return equation_dict



def build_equation(equation_dict:dict, eq_sign:str='=') -> str:
    # 
    left, right = [], []
    for comp, coeff in equation_dict.items():
        if coeff < 0:
            x = comp if coeff==-1 else str(-coeff)+' '+comp
            left.append(x)
        elif coeff > 0:
            x = comp if coeff==1 else str(coeff)+' '+comp
            right.append(x)
        elif coeff == 0:
            left.append(comp)
            right.append(comp)

    equation = ' + '.join(left)+' '+eq_sign.strip()+' '+' + '.join(right)
    return equation



def to_mol_methods():
    methods = {'inchi': inchi_to_mol,
               'smiles': smiles_to_mol,
               'file': file_to_mol,
               'kegg': kegg_compound_to_mol,
               'kegg.compound': kegg_compound_to_mol,
               'metanetx': metanetx_id_to_mol,
               'metanetx.chemical': metanetx_id_to_mol,
               'hmdb': hmdb_id_to_mol,
               'chebi': chebi_id_to_mol,
               'lipidmaps': lipidmaps_id_to_mol,
               'recon3d': recon3d_id_to_mol,
               'inchi-key': inchi_key_to_mol,
               'name': name_to_mol,
               }
    return methods



def to_mol(cid:str, cid_type:str, Hs=True, sanitize=True) -> rdkit.Chem.rdchem.Mol:
    # 
    if not isinstance(cid, str):
        raise InputTypeError('cid must be String type, but got {0}'.format(type(cid)))
    elif not isinstance(cid_type, str):
        raise InputTypeError('cid_type must be String type, but got {0}'.format(type(cid_type)))

    # the main body
    methods = to_mol_methods()

    if cid_type.lower() not in methods.keys():
        raise InputValueError(f'cid_type must be one of {list(methods.keys())}, {cid_type} id cannot be recognized')
    
    cid_type = cid_type.lower()
    if cid_type=='auto':
        _to_mols = methods
    else:
        _to_mols = {cid_type: methods[cid_type]}

    output = {}
    for _cid_type, _to_mol in _to_mols.items():
        try:
            mol = _to_mol(cid)
        except:
            mol = None
        if mol:
            if Hs==True:
                mol = Chem.AddHs(mol)
            elif Hs==False:
                mol = Chem.RemoveHs(mol)
            output[_cid_type] = mol
    if len(output)>1:
        raise ValueError(f'Which {cid} is {tuple(output.keys())}?')
    return tuple(output.values())[0] if output else None



def equation_dict_to_mol_dict(equation_dict:dict, cid_type) -> Dict[rdkit.Chem.rdchem.Mol, float]|None:
    if type(cid_type)==str:
        cid_Types = [cid_type.lower()] * len(equation_dict)
    elif type(cid_type)==list:
        assert len(cid_type)==len(equation_dict)
        cid_Types = cid_type

    mol_dict = dict(map(lambda item: (to_mol(item[0][0], item[1]), item[0][1]),
                        zip(equation_dict.items(), cid_Types)))
    
    return mol_dict if not None in mol_dict else None



def equation_to_mol_dict(equation:str, cid_type:str):
    equation_dict = parse_equation(equation)
    
    return equation_dict_to_mol_dict(equation_dict, cid_type)



def equations_to_S(equations) -> pd.DataFrame:
    comps = set()
    for cs in [parse_equation(x).keys() for x in equations]:
        for c in cs:
            comps.add(c)
    comps = list(comps)
    comps.sort()

    S = pd.DataFrame(index=comps, columns=range(len(equations)), dtype=float, data=0)
    for i in range(len(equations)):
        for comp, coeff in parse_equation(equations[i]).items():
            S.loc[comp,i] = coeff

    return S



def S_to_equations(S, mets) -> list:
    if len(mets) not in S.shape:
        raise ValueError('S.shape not match mets length')
    elif S.shape[0]==S.shape[1]:
        print('S is square, make sure dim 0 match mets')
    elif S.shape[0] == len(mets):
        S = S.T
    elif S.shape[1] == len(mets):
        pass
    else:
        print('S.shape=', S.shape)

    equations = []
    for s in S:
        reactants = np.array(mets)[s!=0]
        coeff = s[s!=0]
        equations.append(build_equation(dict(zip(reactants, coeff))))

    return equations



def normalize_mol(mol):
    #
    mol = rdMolStandardize.Normalize(mol)
    #mol = rdMolStandardize.CanonicalTautomer(mol)
    #mol = rdMolStandardize.Reionize(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    #mol = rdMolStandardize.Cleanup(mol)
    return Chem.AddHs(mol)



def atom_bag(mol:rdkit.Chem.rdchem.Mol) -> Dict[str, int|float]:
    atom_bag = {}
    charge = 0
    for atom in mol.GetAtoms():
        atom_bag[atom.GetSymbol()] = 1 + atom_bag.get(atom.GetSymbol(), 0)
        charge += atom.GetFormalCharge()
    atom_bag['charge'] = charge
        
    return atom_bag



def is_balanced(reaction:Dict[rdkit.Chem.rdchem.Mol|None, float|int], ignore_H_ion=False, ignore_H2O=False) -> bool:
    # 
    if None in reaction.keys():
        return None
    diff_atom = {}
    for mol, coeff in reaction.items():
        for atom, num in atom_bag(mol).items():
            diff_atom[atom] = diff_atom.get(atom, 0) + coeff * num

    if (diff_atom.get('R', 0) + diff_atom.get('*', 0)) == 0:
        diff_atom.pop('R', None)
        diff_atom.pop('*', None)

    if ignore_H_ion:
        diff_atom['H'] = diff_atom.get('H', 0) - diff_atom.get('charge', 0)
        diff_atom['charge'] = 0
    if ignore_H2O:
        diff_atom['H'] = diff_atom.get('H', 0) - 2 * diff_atom.get('O', 0)
        diff_atom['O'] = 0
        
    unbalanced_atom = {}
    for atom, num in diff_atom.items():
        if num!=0:
            unbalanced_atom[atom] = num

    return False if unbalanced_atom else True



def get_pKa_methods():
    methods = {}
    if os.path.isfile(chemaxon_pka_json_path):
        methods = {'chemaxon_pKa_json':get_pka_from_json}
    elif os.path.isfile(chemaxon_pka_csv_path):
        methods.update({'chemaxon_pKa_csv':get_pka_from_file})
    methods.update({'chemaxon':get_pKa_from_chemaxon})
    return methods



def get_pKa(smiles, temperature:float=default_T, source='chemaxon_pKa_json') -> dict:
    # source: 
    methods = get_pKa_methods()
    # the main body of this function
    if source=='auto':
        for source in methods:
            if pKa := methods[source](smiles, temperature=temperature):
                break
    elif source in methods.keys():
        pKa = methods[source](smiles, temperature=temperature)
    else:
        raise InputValueError('source must be one of', methods.keys())
    pKa = deepcopy(pKa)
    if pKa is None:
        return None
    else:
        for xpKa in pKa.values():
            for atom_pKa in xpKa.copy():
                if np.isnan(atom_pKa['value']):
                    xpKa.remove(atom_pKa)
        return pKa



def debye_hueckel(sqrt_ionic_strength: float, T_in_K: float) -> float:
    """Compute the ionic-strength-dependent transformation coefficient.

    For the Legendre transform to convert between chemical and biochemical
    Gibbs energies, we use the extended Debye-Hueckel theory to calculate the
    dependence on ionic strength and temperature.

    Parameters
    ----------
    sqrt_ionic_strength : float
        The square root of the ionic-strength in M
    temperature : float
        The temperature in K


    Returns
    -------
    Quantity
        The DH factor associated with the ionic strength at this
        temperature in kJ/mol

    """
    _a1 = 9.20483e-3  # kJ / mol / M^0.5 / K
    _a2 = 1.284668e-5  # kJ / mol / M^0.5 / K^2
    _a3 = 4.95199e-8  # kJ / mol / M^0.5 / K^3
    B = 1.6  # 1 / M^0.5
    alpha = _a1 * T_in_K - _a2 * T_in_K ** 2 + _a3 * T_in_K ** 3  # kJ / mol
    return alpha / (1.0 / sqrt_ionic_strength + B)  # kJ / mol



def ddGf_to_aqueous(pH: float, pMg: float, I: float, T: float, net_charge: float, num_H: float, num_Mg: float) -> float:
    RT = R * T
    H_term = num_H * RT * np.log(10) * pH

    _dg_mg = (T / default_T) * standard_formation_dg_Mg + (1.0 - T / default_T) * standard_formation_dh_Mg
    Mg_term = num_Mg * (RT * np.log(10) * pMg - _dg_mg)

    if I > 0:
        dh_factor = debye_hueckel(I ** 0.5, T)
        is_term = dh_factor * (net_charge ** 2 - num_H - 4 * num_Mg)
    else:
        is_term = 0.0

    return (H_term + Mg_term - is_term)



def iter_pseudo(dG_dis, d_charge, pH, T, pKa, n):
    # 
    acidicV = dict([(x['atomIndex'],x['value']) for x in pKa['acidicValuesByAtom'] if not np.isnan(x['value'])])
    basicV = dict([(x['atomIndex'],x['value']) for x in pKa['basicValuesByAtom'] if not np.isnan(x['value'])])
    t = list(set(acidicV.keys())|set(basicV.keys()))
    t.sort()

    if n==len(t):
        assert len(dG_dis)==len(d_charge)
        return dG_dis, d_charge
    else:
        i = t[n]
        pka = acidicV.get(i)
        pkb = basicV.get(i)
        if pka is not None:
            dG_dis_a = dG_dis - R * T * np.log(10) * (pH - pka)
            d_charge_a = d_charge - 1
        else:
            dG_dis_a = np.array([])
            d_charge_a = np.array([])
        if pkb is not None:
            dG_dis_b = dG_dis - R * T * np.log(10) * (pkb - pH)
            d_charge_b = d_charge + 1
        else:
            dG_dis_b = np.array([])
            d_charge_b = np.array([])
        dG_dis = np.concatenate([dG_dis, dG_dis_a, dG_dis_b])
        d_charge = np.concatenate([d_charge, d_charge_a, d_charge_b])
        return iter_pseudo(dG_dis, d_charge, pH, T, pKa, n+1)



def pseudoisomers_ddGf(chemaxon_pKa, pH:float, T:float):
    # 
    dG_dis, d_charge = np.array([0]), np.array([0])

    return iter_pseudo(dG_dis, d_charge, pH, T, chemaxon_pKa, 0)[0]



def pseudoisomers_delta_charge(chemaxon_pKa, pH:float, T:float):
    # 
    dG_dis, d_charge = np.array([0]), np.array([0])

    return iter_pseudo(dG_dis, d_charge, pH, T, chemaxon_pKa, 0)[1]



def pseudoisomers_ratio(pKa, pH:float, T:float):
    # 
    RT = R * T
    ddGf_standard_prime_j_array = pseudoisomers_ddGf(pKa, pH, T)
    ratio = np.exp(-ddGf_standard_prime_j_array/RT)/np.sum(np.exp(-ddGf_standard_prime_j_array/RT))

    return ratio



def ddGf_to_dissociation(pH: float, T: float, pKa:dict):
    # 
    RT = R * T
    
    ddGf_standard_prime_j_array = pseudoisomers_ddGf(pKa, pH, T)
    ddGf_standard_prime = -RT * np.log(np.sum(np.exp(-ddGf_standard_prime_j_array/RT)))

    return ddGf_standard_prime



def ddGf_to_elec(charge, e_potential):
    # 
    dg = FARADAY * charge * e_potential
    return dg



def ddGf(compound, condition1, condition2):
    num_H = compound.atom_bag.get('H', 0)
    num_Mg = compound.atom_bag.get('Mg', 0)
    net_charge = compound.atom_bag.get('charge', 0)

    ddGf_1_2 = []
    for condition in [dict(condition1), dict(condition2)]:
        pH = condition.get('pH', default_pH)
        T = condition.get('T', default_T)
        I = condition.get('I', default_I)
        pMg = condition.get('pMg', default_pMg)
        e_potential = condition.get('e_potential', default_e_potential)

        if compound.Smiles == '[H+]':
            return 0

        pKa = compound.pKa(T)
        if pKa == None:
            return None

        ratio = pseudoisomers_ratio(pKa, pH, T)
        d_charge = pseudoisomers_delta_charge(pKa, pH, T)
        charge = net_charge# + np.sum(ratio * d_charge)
        
        num_H_condition = num_H# + np.sum(ratio * d_charge)

        ddGf_1_2.append(ddGf_to_dissociation(pH, T, pKa) + ddGf_to_aqueous(pH, pMg, I, T, charge, num_H_condition, num_Mg) + ddGf_to_elec(charge, e_potential))
    
    return ddGf_1_2[1] - ddGf_1_2[0]



def ddGf_H_ion(compound, condition1, condition2):
    assert compound.Smiles == '[H+]'
    ddGf_1_2 = []
    for condition in [dict(condition1), dict(condition2)]:
        pH = condition.get('pH', default_pH)
        T = condition.get('T', default_T)
        e_potential = condition.get('e_potential', default_e_potential)
        ddGf_1_2.append(R * T * np.log(10) * pH + ddGf_to_elec(1, e_potential))
    return ddGf_1_2[1] - ddGf_1_2[0]



def ddGr(reaction, condition1, condition2):
    ddGf_list = [coeff * ddGf(comp, condition1, condition2) for comp, coeff in reaction.items()]
    if None in ddGf_list:
        print(ddGf_list)
        return None
    return sum(ddGf_list)



def remove_map_num(mol):
    mol = deepcopy(mol)
    [atom.ClearProp('molAtomMapNumber') for atom in mol.GetAtoms()]
    return mol



def reacting_map_num_of_rxn(rxn):
    rxn.Initialize()
    map_num = []
    for mol, atoms_idx in zip(rxn.GetReactants(), rxn.GetReactingAtoms()):
        map_num = map_num + [mol.GetAtomWithIdx(idx).GetAtomMapNum() for idx in atoms_idx]
    assert len(set(map_num))==len(map_num)
    return tuple(map_num)



def map_num_to_idx_with_radius(mol:rdkit.Chem.rdchem.Mol, map_num:tuple, radius:int) -> tuple:
    atoms_idx_radius = [set([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in map_num]), ]
    for r in range(radius):
        atoms_idx_r = set()
        for atom_idx_i in atoms_idx_radius[r]:
            for atom_j in mol.GetAtomWithIdx(atom_idx_i).GetNeighbors():
                atoms_idx_r.add(atom_j.GetIdx())
        for exist_atoms in atoms_idx_radius:
            atoms_idx_r = atoms_idx_r - exist_atoms
        atoms_idx_radius.append(atoms_idx_r)

    idxs = set()
    for idxs_set in atoms_idx_radius:
        idxs |= idxs_set

    return tuple(idxs)
    
