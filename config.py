import os

package_path = os.path.split(os.path.abspath(__file__))[0]

# models
inference_model_state_dict_dir = os.path.join(package_path, 'network', 'model_params', 'MP_network')
os.path.join(package_path, 'network', 'best_model_params')

# training data
train_data_path = os.path.join(package_path, 'data', 'TrainingData.csv')

# compound databases
kegg_database_path = os.path.join(package_path, 'data', 'KEGG')
kegg_compound_data_path = os.path.join(package_path, 'data', 'kegg_compound')
kegg_additions_csv_path = os.path.join(package_path, 'data', 'kegg_compound', 'kegg_additions.csv')
metanetx_database_path = os.path.join(package_path, 'data', 'MetaNetX', 'chem_prop.tsv')
hmdb_database_path = os.path.join(package_path, 'data', 'HMDB', 'structures.csv')
chebi_database_path = os.path.join(package_path, 'data', 'libChEBI')
lipidmaps_database_path = os.path.join(package_path, 'data', 'LIPID_MAPS', 'lipidmaps_ids_cc0.tsv')
recon3d_mol_dir_path = os.path.join(package_path, 'data', 'Recon3D', 'mol')
chemaxon_pka_csv_path = os.path.join(package_path, 'data', 'chemaxon_pKa.csv')
chemaxon_pka_json_path = os.path.join(package_path, 'data', 'chemaxon_pKa.json')
chemaxon_jar_dir = "/opt/chemaxon/jchemsuite/lib" 