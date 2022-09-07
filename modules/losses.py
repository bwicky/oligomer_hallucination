# losses module

import numpy as np
import sys; sys.path.append('/projects/ml/alphafold/alphafold_git/')
from alphafold.common import protein
# dssp loss imports
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
# to run tmalign
import subprocess
from scipy import linalg


######################################################
# COMMON FUNCTIONS USED BY DIFFERENT LOSSES.
######################################################

def get_coord(atom_type, oligo_object):
    '''
    General function to get the coordinates of an atom type in a pdb. For geometric-based losses.
    Returns an array [[chain, resid, x, y, z]]
    '''
    coordinates = []
    pdb_lines = protein.to_pdb(oligo_object.try_unrelaxed_structure).split('\n')
    for l in pdb_lines: # parse PDB lines and extract atom coordinates
        if 'ATOM' in l and atom_type in l:
            s = l.split()
            if len(s[4]) > 1: # residue idx and chain id are no longer space-separated at high id values
                coordinates.append([s[4][0], int(s[4][1:]), np.array(s[5:8], dtype=float)])
            else:
                coordinates.append([s[4], int(s[5]), np.array(s[6:9], dtype=float)])

    coord = np.array(coordinates, dtype=object)

    # Find chain breaks.
    ch_breaks = np.where(np.diff(coord[:, 1]) > 1)[0]
    ch_ends = np.append(ch_breaks, len(coord) - 1)
    ch_starts = np.insert(ch_ends[:-1], 0, 0)

    chain_list = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for k, start_finish in enumerate(list(zip(ch_starts, ch_ends))):
        coord[start_finish[0] + 1 : start_finish[1]+1 , 0] = chain_list[k] # re-assign chains based on chain breaks

    return coord


def dssp_wrapper(pdbfile):
    '''Compute DSSP string on structure.'''

    dssp_tuple = dssp_dict_from_pdb_file(pdbfile, DSSP="/home/lmilles/lm_bin/dssp")[0]

    dssp_list = []
    for key in dssp_tuple.keys():
        dssp_list.append(dssp_tuple[key][2])

    return dssp_list

def calculate_dssp_fractions(dssp_list):
    '''Compute DSSP fraction based on a DSSP list.'''

    N_residues = len(dssp_list)
    fraction_beta  = float(dssp_list.count("E") ) / float(N_residues)
    fraction_helix = float(dssp_list.count("H") ) / float(N_residues)
    fraction_other = float(1.0 - fraction_beta-fraction_helix)
    # print(dssp_list, fraction_beta, fraction_helix, fraction_other)

    return fraction_beta, fraction_helix, fraction_other


def tmalign_wrapper(template, temp_pdbfile, force_alignment=None):
    if force_alignment == None:
        p = subprocess.Popen(f'/home/lmilles/lm_bin/TMalign {template} {temp_pdbfile} | grep -E "RMSD|TM-score=" ', stdout=subprocess.PIPE, shell=True)
    else:
        p = subprocess.Popen(f'/home/lmilles/lm_bin/TMalign {template} {temp_pdbfile} -I {force_alignment} | grep -E "RMSD|TM-score=" ', stdout=subprocess.PIPE, shell=True)
    output, __ = p.communicate()
    tm_rmsd  = float(str(output)[:-3].split("RMSD=")[-1].split(",")[0] )
    tm_score = float(str(output)[:-3].split("TM-score=")[-1].split("(if")[0] )

    return tm_rmsd, tm_score


############################
# LOSS COMPUTATION
############################

def compute_loss(losses, oligo, args, loss_weights):
    '''
    Compute the loss of a single oligomer.
    losses: list of list of losses and their associated arguments (if any).
    oligo: an Oligomer object.
    args: the whole argument namespace (some specific arguments are required for some specific losses).
    loss_weights: list of weights associated with each loss.
    '''
    # intialize scores
    scores = []
    # iterate over all losses
    for loss_idx, current_loss in enumerate(losses) :
        loss_type, loss_params  = current_loss # assign loss and its arguments if any.

        if loss_type == 'plddt':
            # NOTE:
            # Using this loss will optimise plddt (predicted lDDT) for the sequence(s).
            # Early benchmarks suggest that this is not appropriate for forcing the emergence of complexes.
            # Optimised sequences tend to be folded (or have good secondary structures) without forming inter-chain contacts.
            score = 1. - np.mean(oligo.try_prediction_results['plddt'])


        elif loss_type == 'ptm':
            # NOTE:
            # Using this loss will optimise ptm (predicted TM-score) for the sequence(s).
            # Early benchmarks suggest that while it does force the apparition of inter-chain contacts,
            # it might be at the expense of intra-chain contacts and therefore folded protomer structures.
            score = 1. - np.mean(oligo.try_prediction_results['ptm'])


        elif loss_type == 'pae':
            # NOTE:
            # Using this loss will optimise the mean of the pae matrix (predicted alignment error).
            # This loss has not been properly benchmarked, but some early results suggest that it might suffer from the same problem as ptm.
            # During optimisation, off-digonal contacts (inter-chain) may get optimsed at the expense of the diagonal elements (intra-chain).
            norm = np.mean(oligo.init_prediction_results['predicted_aligned_error'])
            score = np.mean(oligo.try_prediction_results['predicted_aligned_error']) / norm


        elif loss_type == 'entropy':
            # CAUTION:
            # This loss is unlikely to yield anything useful at present.
            # i,j pairs that are far away from each other, or for which AF2 is unsure, have max prob in the last bin of their respective distograms.
            # This will generate an artifically low entropy for these positions.
            # Need to find a work around this issue before using this loss.
            print('Entropy definition is most likley improper for loss calculation. Use at your own risk...')

            # Distogram from AlphaFold2 represents pairwise Cb-Cb distances, and is outputted as logits.
            probs = softmax(oligo.try_prediction_results['distogram']['logits'], -1. ) # convert logit to probs

            score = np.mean(-np.sum((np.array(probs) * np.log(np.array(probs))), axis=-1))


        elif loss_type == 'dual':
            # NOTE:
            # This loss jointly optimises ptm and plddt (equal weights).
            # It attemps to combine the best of both worlds -- getting folded structures that are in contact.
            # This loss is currently recommended unless cyclic geometries are desired (tends to generate linear oligomers).
            score = 1. - (np.mean(oligo.try_prediction_results['plddt']) / 2.) - (oligo.try_prediction_results['ptm'] / 2.)


        elif loss_type == 'pae_sub_mat':
            # NOTE:
            # This loss optimises the mean of the pae sub-matrices means.
            # The value of the loss will be different to pae in the case of hetero-oligomers that have chains of different lenghts, but identical otherwise.
            # The mean of the sub matrices' means is different from the overall mean if the sub matrices don't all have the same shape.

            sub_mat_init = []
            sub_mat = []
            prev1, prev2 = 0, 0
            for L1 in oligo.chain_Ls:
                Lcorr1 = prev1 + L1

                for L2 in oligo.chain_Ls:
                    Lcorr2 = prev2 + L2
                    sub_mat_init.append(oligo.init_prediction_results['predicted_aligned_error'][prev1:Lcorr1, prev2:Lcorr2]) # means of the initial sub-matrices
                    sub_mat.append(oligo.try_prediction_results['predicted_aligned_error'][prev1:Lcorr1, prev2:Lcorr2]) # means of the tried move sub-matrices
                    prev2 = Lcorr2

                prev2 = 0
                prev1 = Lcorr1

            norm =  np.mean([np.mean(sub_m) for sub_m in sub_mat_init])
            score = np.mean([np.mean(sub_m) for sub_m in sub_mat]) / norm


        elif loss_type == 'pae_asym':
            # NOTE:
            # This loss has different weights associated to the means of the different PAE sub matrices (asymmetric weighting).
            # The idea is enforcing loss optimisation for adjacent units to force cyclisation.
            # Off-diagonal elements (+/-1 from the diagaonl, and opposite corners) have higher weights.
            # The weight correction is scaled with the shape of the matrix of sub matrices.
            # By default is scales so that the re-weighted terms count as much as the rest (irrespective of the size of the matrix of sub matrices)

            contribution = 1 # if set to one, the re-weighting is done such that diagonal/corner elements count as much as the rest.

            sub_mat_means_init = []
            sub_mat_means = []
            prev1, prev2 = 0, 0
            for L1 in oligo.chain_Ls:
                Lcorr1 = prev1 + L1

                for L2 in oligo.chain_Ls:
                    Lcorr2 = prev2 + L2
                    sub_mat_means_init.append(np.mean(oligo.init_prediction_results['predicted_aligned_error'][prev1:Lcorr1, prev2:Lcorr2])) # means of the initial sub-matrices
                    sub_mat_means.append(np.mean(oligo.try_prediction_results['predicted_aligned_error'][prev1:Lcorr1, prev2:Lcorr2])) # means of the tried move sub-matrices
                    prev2 = Lcorr2

                prev2 = 0
                prev1 = Lcorr1

            w_corr = contribution * (oligo.oligo_L**2) / (2. * oligo.oligo_L) # correction scales with the size of the matrix of sub matrices and the desired contribution.

            # Weight matrix
            W = np.ones((len(oligo.subunits), len(oligo.subunits)))
            W[0, -1] = 1 * w_corr
            W[-1, 0] = 1 * w_corr
            W[np.where(np.eye(*W.shape, k=-1) == 1)] = 1 * w_corr
            W[np.where(np.eye(*W.shape, k=1) == 1)] = 1 * w_corr

            norm = np.mean( W * np.reshape(sub_mat_means_init, (len(oligo.subunits), len(oligo.subunits))))
            score = np.mean( W * np.reshape(sub_mat_means, (len(oligo.subunits), len(oligo.subunits)))) / norm


        elif loss_type == 'dual_cyclic':
            # NOTE:
            # This loss is based on dual, but adds a geometric term that forces a cyclic symmetry.
            # At each step the PDB is generated, and the distance between the center of mass of adjacent units computed.
            # The standard deviation of these neighbour distances is added to the loss.

            c = get_coord('CA', oligo) # get CA atoms

            # Compute center of mass (CA) of each chain.
            chains = set(c[:,0])
            center_of_mass = {ch:float for ch in chains}
            for ch in chains:
                center_of_mass[ch] = np.mean(c[c[:,0]==ch][:,2:], axis=0)[0]

            # Compare distances between adjacent chains, including first-last.
            chain_order = sorted(center_of_mass.keys())
            next_chain = np.roll(chain_order, -1)

            proto_dist = []
            for k, ch in enumerate(chain_order):
                proto_dist.append(np.linalg.norm(center_of_mass[next_chain[k]]-center_of_mass[ch])) # compute separation distances.

            separation_std = np.std(proto_dist) # the standard deviation of the distance separations.

            # Compute the score, which is an equal weighting between plddt, ptm and the geometric term.
            score = 1. - (np.mean(oligo.try_prediction_results['plddt']) / 2. )  \
                       - (oligo.try_prediction_results['ptm'] / 2.)  \
                       +  separation_std

        elif loss_type == 'cyclic':
            # NOTE:
            # At each step the PDB is generated, and the distance between the center of mass of adjacent units computed.
            # The standard deviation of these neighbour distances is added to the loss.

            c = get_coord('CA', oligo) # get CA atoms

            # Compute center of mass (CA) of each chain.
            chains = set(c[:,0])
            center_of_mass = {ch:float for ch in chains}
            for ch in chains:
                center_of_mass[ch] = np.mean(c[c[:,0]==ch][:,2:], axis=0)[0]

            # Compare distances between adjacent chains, including first-last.
            chain_order = sorted(center_of_mass.keys())
            next_chain = np.roll(chain_order, -1)

            proto_dist = []
            for k, ch in enumerate(chain_order):
                proto_dist.append(np.linalg.norm(center_of_mass[next_chain[k]]-center_of_mass[ch])) # compute separation distances.

            separation_std = np.std(proto_dist) # the standard deviation of the distance separations.

            # Compute the score, which is an equal weighting between plddt, ptm and the geometric term.
            score = 0.0 + separation_std

        elif loss_type == "tmalign":
            # NOTE:
            # a loss to enforce tmscore against a given template.

            # write temporary pdbfile to compute tmscore.
            temp_pdbfile = f'{args.out}_models/tmp.pdb'
            with open( temp_pdbfile , 'w') as f:
                f.write( protein.to_pdb(oligo.try_unrelaxed_structure) )

            force_alignment = None

            tm_rmsd, tm_score = tmalign_wrapper(args.template, temp_pdbfile, args.template_alignment)
            print("   tm_RMSD, tmscore " , tm_rmsd, tm_score)

            score = 1. - tm_score


        elif loss_type == "dual_tmalign":
            # NOTE:
            # This loss jointly optimises plddt, ptm, and tmscore against a template (equal weights).

            # write temporary pdbfile to compute tmscore.
            temp_pdbfile = f'{args.out}_models/tmp.pdb'
            with open( temp_pdbfile , 'w') as f:
                f.write( protein.to_pdb(oligo.try_unrelaxed_structure) )

            tm_rmsd, tm_score = tmalign_wrapper(args.template, temp_pdbfile, args.template_alignment)
            print("   tm_RMSD, tmscore" , tm_rmsd, tm_score)

            score = 1. - (np.mean(oligo.try_prediction_results['plddt']) / 3.) - (oligo.try_prediction_results['ptm'] / 3.) - (tm_score / 3.)


        elif loss_type == "frac_dssp":
            # NOTE:
            # a loss to enforce a fraction of secondary structure elements (e.g. 80% of fold must be beta sheet) is enforced

            # write temporary pdbfile for computing DSSP.
            temp_pdbfile = f'{args.out}_models/tmp.pdb'
            with open( temp_pdbfile , 'w') as f:
                f.write( protein.to_pdb(oligo.try_unrelaxed_structure) )

            dssp_list = dssp_wrapper(temp_pdbfile)
            fraction_beta, fraction_helix, fraction_other = calculate_dssp_fractions(dssp_list)
            print (f" fraction E|H|notEH: {fraction_beta:2.2f} | {fraction_helix:2.2f} | {fraction_other:2.2f}")
            
            measured_fractions = {"E" : fraction_beta, "H" : fraction_helix, "L": fraction_other }
            
            # DSSP assigns eight states: 310-helix (represented by G), alpha-helix (H), pi-helix (I), helix-turn (T), extended beta sheet (E), beta bridge (B), bend (S) and other/loop (L).

            frac_delta = []
            for specified_secstruc in args.dssp_fractions_specified.keys() :
                if args.dssp_fractions_specified[specified_secstruc] != None:
                    frac_delta.append( np.abs( 
                        args.dssp_fractions_specified[specified_secstruc] - measured_fractions[specified_secstruc] ) 
                     ) 

            score = 0. + np.mean(frac_delta)   

        elif loss_type == "min_frac_dssp":
            # NOTE:
            # a loss to enforce a minium fraction of secondary structure elements (e.g. loss is maximal if more than 80% of fold must be beta sheet) is enforced
            # IN DEVELOPMENT

            # write temporary pdbfile for computing DSSP.
            temp_pdbfile = f'{args.out}_models/tmp.pdb'
            with open( temp_pdbfile , 'w') as f:
                f.write( protein.to_pdb(oligo.try_unrelaxed_structure) )

            dssp_list = dssp_wrapper(temp_pdbfile)
            fraction_beta, fraction_helix, fraction_other = calculate_dssp_fractions(dssp_list)
            print (f" fraction E|H|notEH: {fraction_beta:2.2f} | {fraction_helix:2.2f} | {fraction_other:2.2f}")
            
            measured_fractions = {"E" : fraction_beta, "H" : fraction_helix, "L": fraction_other }
            
            # DSSP assigns eight states: 310-helix (represented by G), alpha-helix (H), pi-helix (I), helix-turn (T), extended beta sheet (E), beta bridge (B), bend (S) and other/loop (L).

            frac_delta = []
            for specified_secstruc in args.dssp_fractions_specified.keys() :
                if args.dssp_fractions_specified[specified_secstruc] != None:
                    
                    curr_delta = args.dssp_fractions_specified[specified_secstruc] - measured_fractions[specified_secstruc] 
                    if curr_delta < 0.0 :
                        curr_delta = 0.0
                    frac_delta.append( np.abs( curr_delta )) 

            score = 0. + np.mean(frac_delta)  

        elif loss_type == 'aspect_ratio':
            # NOTE:
            # This loss adds a geometric term that forces an aspect ratio of 1 (spherical protomers) to prevent extended structures.
            # At each step, the PDB is generated, and a singular value decomposition is performed on the coordinates of the CA atoms.
            # The ratio of the two largest values is taken as the aspect ratio of the protomer.
            # For oligomers, the aspect ratio is calculated for each protomer independently, and the average returned.

            c = get_coord('CA', oligo) # get CA atoms, returns array [[chain, resid, x, y, z]]
            aspect_ratios = []
            chains = set(c[:,0])
            for ch in chains:
                coords = np.array([a[2:][0] for a in c[c[:,0]==ch]])
                coords -= coords.mean(axis=0) # mean-center the protomer coordinates
                s = linalg.svdvals(coords) # singular values of the coordinates
                aspect_ratios.append(s[1] / s[0])

            score = 1. - np.mean(aspect_ratios) # average aspect ratio across all protomers of an oligomer

        else:
            sys.exit("specified loss:  {loss_type}  \n not found in available losses!")

        scores.append(score)

    # Normalize loss weights vector.
    loss_weights_normalized = np.array(loss_weights) / np.sum(loss_weights)

    # Total loss for this oligomer is the sum of its weighted scores.
    final_score = np.sum(np.array(scores) * loss_weights_normalized)

    # The loss counts positively or negatively to the overall loss depending on whether this oligomer is positively or negatively designed.
    if oligo.positive_design == True:
        loss = float(final_score)
    else:
        loss = float(final_score) + 1

    return loss
