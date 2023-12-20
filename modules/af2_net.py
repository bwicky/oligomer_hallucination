# af2 network module
import sys
import mock
import numpy as np
from timeit import default_timer as timer
from jax.lib import xla_bridge
print(xla_bridge.get_backend().platform)

# AF2-specific libraries
try:
    import alphafold
except:
    # change this to your local clone
    sys.path.append('/projects/ml/alphafold/alphafold_git/')
from alphafold.common import protein
from alphafold.data import pipeline, templates
from alphafold.data.parsers import Msa
from alphafold.model import data, config, model
from alphafold.relax import relax


def setup_models(oligo_names, model_id=4, recycles=1, msa_clusters=1):
    '''Setup AlphaFold2 models.'''

    mod = f'model_{model_id}_ptm'  # _ptm series of models necessary for returning pTM and PAE.
    print ("Setting up model: " , mod)
    model_config = config.model_config(mod)
    model_config.model.num_recycle = recycles # AF2 default is 3. Effect on computing time is linear.
    model_config.data.common.num_recycle = recycles # AF2 default is 3. Effect on computing time is linear.
    # model 1,2 do not accept msa_clusters < 512, manually correct:
    min_msa_clusters_model1and2 = 5
    if any(str(m) in mod for m in [1,2]) and msa_clusters < min_msa_clusters_model1and2 :
        print( f"Setting MSA clusters to {min_msa_clusters_model1and2} as using: ", mod)
        msa_clusters = int(min_msa_clusters_model1and2)


    model_config.data.common.max_extra_msa = msa_clusters  # AF2 default is 5120. Turning off is about 8x faster.
    model_config.data.eval.max_msa_clusters = msa_clusters # AF2 default is 512. Turning off is about 8x faster.
    model_config.data.eval.num_ensemble = 1
    
    model_params = data.get_model_haiku_params(model_name=mod, data_dir='/projects/ml/alphafold')
    model_runner = model.RunModel(model_config, model_params)
    # model_runner.model_params = model_params

    # Make separate model runners for each oligomer to avoid re-compilation and save time.
    model_runners = {}
    for oligo in oligo_names:
        model_runners[oligo] = model_runner

    return model_runners


def mk_mock_template(query_sequence):
    '''Generate mock template features from the input sequence.'''

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
    templates_aatype = templates.residue_constants.sequence_to_onehot(output_templates_sequence,
                                                                    templates.residue_constants.HHBLITS_AA_TO_ID)

    template_features = {'template_all_atom_positions': np.array(templates_all_atom_positions)[None],
                         'template_all_atom_masks': np.array(templates_all_atom_masks)[None],
                         'template_sequence': [f'none'.encode()],'template_aatype': np.array(templates_aatype)[None],
                         'template_confidence_scores': np.array(output_confidence_scores)[None],
                         'template_domain_names': [f'none'.encode()],
                         'template_release_date': [f'none'.encode()]
                        }

    return template_features


def predict_structure(oligo_object,
                    single_chain,
                    model_runner: model.RunModel,
                    random_seed=0):
    '''Predicts structure for a given oligomer using AlphaFold2.'''

    query_sequence = oligo_object.try_seq
    Ls = oligo_object.chain_Ls

    # Mock pipeline.
    data_pipeline_mock = mock.Mock()
    data_pipeline_mock.process.return_value = {
        **pipeline.make_sequence_features(sequence=query_sequence,
                                          description="none",
                                          num_res=len(query_sequence)),
        **pipeline.make_msa_features([Msa(sequences=[query_sequence],
                                     deletion_matrix=[[0]*len(query_sequence)],
                                     descriptions=["none"])]),
        **mk_mock_template(query_sequence)
    }

    # Get features.
    feature_dict = data_pipeline_mock.process()

    if single_chain == False:
        # Add big enough number to residue index to indicate chain breaks.
        idx_res = feature_dict['residue_index']
        L_prev = 0
        # Ls: number of residues in each chain.
        for L_i in Ls[:-1]:
            idx_res[L_prev+L_i:] += 200
            L_prev += L_i
        feature_dict['residue_index'] = idx_res

    # Run AlphaFold2 prediction.
    start = timer()
    processed_feature_dict = model_runner.process_features(feature_dict, random_seed=random_seed)
    prediction_results = model_runner.predict(processed_feature_dict, random_seed=random_seed)
    unrelaxed_protein = protein.from_prediction(processed_feature_dict, prediction_results)
    end = timer()
    # scale pLDDT to be between 0 and 1
    prediction_results['plddt'] = prediction_results['plddt'] / 100.0

    print(f'({oligo_object.name} prediction took {(end - start):.1f} s)')

    return prediction_results, unrelaxed_protein


def amber_relax(unrelaxed_protein):
    '''Relax the prediction, return Amber-relaxed protein'''

    start = timer()
    amber_relaxer = relax.AmberRelaxation(max_iterations=0,tolerance=2.39,
                                            stiffness=10.0,exclude_residues=[],
                                            max_outer_iterations=20)
    relaxed_protein, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
    end = timer()

    print(f' AMBER relax took {(end - start):.1f} s')

    return relaxed_protein
