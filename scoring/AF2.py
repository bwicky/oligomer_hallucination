#!/software/conda/envs/SE3/bin/python
# Predict structures from a FASTA file with AlphaFold2.
# A bunch of scores are appended at the end of the generated PDBs.
# Basile Wicky -- April 21, 2022

import sys
from time import time
import argparse
import mock
import numpy as np
# import pickle
# from typing import Dict
from Bio import SeqIO
from jax.lib import xla_bridge
print('Using:', xla_bridge.get_backend().platform.upper())

# AlphaFold2 imports
sys.path.append('/projects/ml/alphafold/alphafold_git/')
from alphafold.common import protein
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model
from alphafold.relax import relax


# ARGUMENTS
models = ['model_1_ptm', 'model_2_ptm', 'model_3_ptm', 'model_4_ptm', 'model_5_ptm']

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=" * Predict structres from single sequences using AlphaFold2.\n"
                    " * A series of scores get appended at the end of the generated PDBs.\n"
                    " * Extra outptus can be optionally generated.\n"
        )
parser.add_argument(
        '-f', '--fasta',
        help='FASTA file containing the sequence(s) to be predicted. For oligomers, chains need to be separated by a forward slash (/)',
        action='store',
        type=str,
        required=True
        )
parser.add_argument(
        '-m', '--models',
        help='Model(s) used for prediction (Default: all).',
        action='store',
        type=str,
        choices=models,
        nargs='+',
        default=['model_1_ptm', 'model_2_ptm', 'model_3_ptm', 'model_4_ptm', 'model_5_ptm']
        )
parser.add_argument(
        '-r', '--recycle',
        help='Number of recyles through the network (Default: 3).',
        action='store',
        type=int,
        default=3
        )
parser.add_argument(
        '--relax',
        help='Peform Amber relax on the outputed structure. This is slow.',
        action='store_true',
        )
parser.add_argument(
        '--plddt',
        help='Output the pLDDT array as a .npz file.',
        action='store_true',
        )
parser.add_argument(
        '--write_plddt',
        help='Write plDDT array to PDB file.',
        action='store_true',
        )
parser.add_argument(
        '--pae',
        help='Output the PAE array as a .npz file.',
        action='store_true',
        )
parser.add_argument(
        '--write_pae',
        help='Write PAE array to PDB file.',
        action='store_true',
        )
parser.add_argument(
        '--distogram',
        help='Output the distogram array as a .npz file.',
        action='store_true',
        )

args = parser.parse_args()

# FUNCTIONS
def make_mock_template_features(
    query_sequence: str) -> dict:
    ''' Make mock template features.'''
    
    output_templates_sequence = []
    output_confidence_scores = []
    templates_all_atom_positions = []
    templates_all_atom_masks = []
    
    for _ in query_sequence:
        templates_all_atom_positions.append(np.zeros((templates.residue_constants.atom_type_num, 3)))
        templates_all_atom_masks.append(np.zeros(templates.residue_constants.atom_type_num))
        output_templates_sequence.append('-')
        output_confidence_scores.append(-1)
    
    output_templates_sequence = ''.join(output_templates_sequence)
    templates_aatype = templates.residue_constants.sequence_to_onehot(output_templates_sequence, templates.residue_constants.HHBLITS_AA_TO_ID)

    template_features = {
        'template_all_atom_positions': np.array(templates_all_atom_positions)[None],
        'template_all_atom_masks': np.array(templates_all_atom_masks)[None],
        'template_sequence': [f'none'.encode()],
        'template_aatype': np.array(templates_aatype)[None],
        'template_confidence_scores': np.array(output_confidence_scores)[None],
        'template_domain_names': [f'none'.encode()],
        'template_release_date': [f'none'.encode()]
    }

    return template_features


def gen_seq_features(
    query_sequence: str) -> dict:
    '''Generate features for given sequence.'''
    
    query_sequence = query_sequence.replace("/", " ")
    query_sequence = query_sequence.split()
    Ls = list()
    for seq in query_sequence:
        Ls.append(len(seq))
    
    query_sequence = "".join(query_sequence)

    # Mock pipeline for testing.
    data_pipeline_mock = mock.Mock()
    data_pipeline_mock.process.return_value = {
        **pipeline.make_sequence_features(sequence=query_sequence,
                                          description='none',
                                          num_res=len(query_sequence)),
        **pipeline.make_msa_features(msas=[[query_sequence]],
                                     deletion_matrices=[[[0]*len(query_sequence)]]),
        **make_mock_template_features(query_sequence)
    }
    
    # Get features.
    feature_dict = data_pipeline_mock.process()

    # Add large number to residue index to indicate chain breaks.
    idx_res = feature_dict['residue_index']
    L_prev = 0
    for L_i in Ls[:-1]:
        idx_res[L_prev+L_i:] += 200
        L_prev += L_i
    
    feature_dict['residue_index'] = idx_res
    
    return feature_dict, Ls


def predict_structure(
    name: str,
    seq_feature_dict: dict,
    model_runner: model.RunModel,
    random_seed: int):
    '''Predict structure with AlphaFold2.'''

    # Run prediction.   
    start_pred = time()
    
    # Process the sequence features.
    processed_feature_dict = model_runner.process_features(seq_feature_dict, random_seed=random_seed)
    
    # Generate the prediction.
    prediction_result = model_runner.predict(processed_feature_dict)
    
    # Generate the protein object.
    unrelaxed_protein = protein.from_prediction(processed_feature_dict, prediction_result)

    print( f"Predicting {name} w/ {model_name} took: {time() - start_pred:.1f} s" )
    
    return prediction_result, unrelaxed_protein


def gen_scores(
    prediction_result: dict,
    Ls: list) -> str:
    '''Generate score string to be appended to the PDB.'''

    # Extract prediction arrays.
    scores = {
        'plddt':prediction_result['plddt'],
        'ptm':prediction_result['ptm'],
        'pae':prediction_result['predicted_aligned_error']
    }

    # PAE off-diagonal matrices.
    sub_mats = []
    prev1, prev2 = 0, 0
    for L1 in Ls:
        Lcorr1 = prev1 + L1

        for L2 in Ls:
            Lcorr2 = prev2 + L2
            sub_mats.append(scores['pae'][prev1:Lcorr1, prev2:Lcorr2])
            prev2 = Lcorr2

        prev2 = 0
        prev1 = Lcorr1

    off_diag = [sub_m for sub_m in np.arange(len(Ls)**2) if sub_m%(len(Ls)+1)!=0] # exclude diagonal matrices
    if off_diag != []:
        scores['pae_off_diag'] = np.hstack([m.flatten() for m in np.array(sub_mats)[off_diag]]) # flatten all off-diagonal matrices
        
    # PAE inter-chain contacts.
    pae_contacts_names = []
    for dist in [8.0, 10.0, 12.0]:
        ij_contacts = get_interface_ij(prediction_result, Ls, dist)
        pae_contacts = []
        for ij in ij_contacts:
            pae_contacts.append(scores['pae'][ij])
        
        pae_contacts_name = f'pae_contacts_{dist:.0f}A'
        pae_contacts_names.append(pae_contacts_name)
        scores[pae_contacts_name] = pae_contacts

    # Compute statistics.
    not_ptm = [sc for sc in scores.keys() if sc!='ptm']
    for sc_name in not_ptm:
        
        if len(scores[sc_name]) == 0:
            pass
        
        else:
            
            if sc_name in pae_contacts_names:
                scores[sc_name + '_num'] = len(scores[sc_name])
                
            scores[sc_name + '_mean'] = np.mean(scores[sc_name])
            scores[sc_name + '_median'] = np.median(scores[sc_name])
            scores[sc_name + '_std'] = np.std(scores[sc_name])
            sc_min, sc_max = np.min(scores[sc_name]), np.max(scores[sc_name]) 
            scores[sc_name + '_min'] = sc_min
            scores[sc_name + '_max'] = sc_max
            scores[sc_name + '_range'] = sc_max - sc_min

    # Generate score string.
    score_str = ''
    for sc_name, sc in scores.items():
        
        if sc_name == 'plddt':
            if args.write_plddt:
                score_str += 'plddt_array ' + ','.join(np.array(scores['plddt'], dtype=str)) + '\n'
        
        elif sc_name == 'pae':
            if args.write_pae:
                score_str += 'pae_array ' + ','.join(np.array(scores['pae'], dtype=str).flatten()) + '\n'
        
        elif sc_name == 'pae_off_diag':
            pass
        
        elif sc_name in pae_contacts_names:
            pass
            
        else:
            score_str += sc_name + ' ' + str(sc) + '\n'

    return score_str
    

def get_interface_ij(
    prediction_result: dict, 
    Ls: list,
    d: float) -> list:
    '''Return i,j pairs between chains that are under the distance cutoff.'''

    L0 = 0
    chains_xyz = {}
    for ch_ID, L in enumerate(Ls):
        chains_xyz[ch_ID] = prediction_result['structure_module']['final_atom_positions'][L0:L0+L,1,:] # 3rd index of atom_types is CB
        L0 += L
    
    ij_pairs = []
    for ch1_id, xyz_chain1 in chains_xyz.items():
        
        for ch2_id, xyz_chain2 in chains_xyz.items(): 

            if ch1_id != ch2_id: # only consider inter-chain interactions
                
                dist = np.linalg.norm(xyz_chain1.reshape(-1, 3)[:, None] - xyz_chain2.reshape(-1, 3)[None, :], axis=-1)
                i, j = np.where(dist <= d) # NB only consider CB-CB distance
                ij_pairs += list(zip(i+np.sum(Ls[:ch1_id], dtype=int), j+np.sum(Ls[:ch2_id], dtype=int)))
    
    return ij_pairs



# Make model runners: Models [1,2] and [3,4,5] have the same number of parameters, so only need 2 runners.
model_runners = {}
for model_name in ["model_1_ptm", "model_3_ptm"]:
    model_config = config.model_config(model_name)
    model_config.data.eval.num_ensemble = 1 # Huge impact on compute time. From the paper they say it does not really improve prediction and that they use 1.      
    model_config.model.num_recycle = args.recycle # AF2 default is 3. Effect on computing time is linear.
    model_config.data.common.num_recycle = args.recycle # AF2 default is 3. Effect on computing time is linear.
    model_config.data.common.max_extra_msa = 1  # AF2 default is 5120. Not sure how much faster this makes prediction.
    model_params = data.get_model_haiku_params(model_name=model_name, data_dir="/projects/ml/alphafold")
    model_runner = model.RunModel(model_config, model_params)
    model_runners[model_name] = model_runner

# Model parameters.
model_params = {}
for model_name in models:
    model_param = data.get_model_haiku_params(model_name=model_name, data_dir="/projects/ml/alphafold")
    model_params[model_name] = model_param


# START HERE <---
fasta = list(SeqIO.parse(args.fasta, 'fasta'))
sequences = {}
for f in fasta:
    sequences[f.id] = str(f.seq)

# Order sequences by length to avoid un-necessary re-compliations.
sorted_sequences = {k:v for k, v in sorted(sequences.items(), key=lambda item:len(item[1]))}
print('Number of structures to predict:', len(sorted_sequences))
print('Number of models used for prediction:', len(args.models))
   

# Generate sequence features.
featurized_sequences = {}
for name, sequence in sorted_sequences.items():
    featurized_sequences[name] = gen_seq_features(sequence)
    


# Cycle over the models.
for model_name in args.models:
    
    print(f' >>> Predicting sequence(s) with {model_name} ...')
    start_model = time()

    # Swap params to avoid recompiling.
    if any(str(m) in model_name for m in [1,2]): model_runner = model_runners['model_1_ptm']
    if any(str(m) in model_name for m in [3,4,5]): model_runner = model_runners['model_3_ptm']
    model_runner.params = model_params[model_name]
    
    # Cycle over the sequences.
    for name, featurized_sequence in featurized_sequences.items():
    
        prediction_result, unrelaxed_protein = predict_structure(
                                    name=name,
                                    seq_feature_dict=featurized_sequence[0],
                                    model_runner=model_runner,
                                    random_seed=np.random.randint(42)
        )
        
        scores = gen_scores(prediction_result, featurized_sequence[1])
        
        if args.plddt:
            np.savez_compressed(f'{name}_plddt.npz', a=prediction_result['plddt'])
        
        if args.pae:
            np.savez_compressed(f'{name}_pae.npz', a=prediction_result['predicted_aligned_error'])
            
        if args.distogram:
            np.savez_compressed(f'{name}_dist.npz', a=prediction_result['distogram'])

        if args.relax: # Relax the prediction. 
            start_relax = time()
            amber_relaxer = relax.AmberRelaxation(max_iterations=0, tolerance=2.39,
                                                stiffness=10.0, exclude_residues=[],
                                                max_outer_iterations=20)      
            relaxed_pdb, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
            
            print(f'Amber relax took: {time() - start_relax} s')
            
            with open(f'{name}_relaxed_{model_name}.pdb', 'w') as f:
                f.write(relaxed_pdb)
                f.write(scores)
        
        else:
            test=protein.to_pdb(unrelaxed_protein)
            with open(f'{name}_unrelaxed_{model_name}.pdb', 'w') as f:
                f.write(protein.to_pdb(unrelaxed_protein))
                f.write(scores)

        
    print(f'Predictions with {model_name} took: {time() - start_model} s')
    
print('Done!')

