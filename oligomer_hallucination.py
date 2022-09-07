#!/software/conda/envs/SE3/bin/python

## Script for performing multistate-state design using AlphaFold2 MCMC hallucination.
## <bwicky@uw.edu> and <lmilles@uw.edu>
## Started: 2021-08-11
## Re-factored: 2021-08-20

#######################################
# LIBRARIES
#######################################

import os, sys
import numpy as np
import copy

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir + '/modules/') # import modules

from arg_parser import *
from mutations import *
from af2_net import *
from losses import *


########################################
# PROTOMERS / OLIGOMER CLASSES
########################################

class Protomers:
    '''
    Object for keeping track of protomer sequences during hallucination.
    Contains the mutation method function.
    unique_protomer: set of unique letters, one for each protomer
    lengths: length of each protomer (must be in alphabetic order)
    aa_freq: dictonary containing the frequencies of each aa
    sequences: sequences corresponding to each protomer (in alphabetic order). Optional.
    '''

    def __init__(self, unique_protomers=[], lengths=[], aa_freq={}, sequences=None, position_weights=None):

        # Initialise sequences.
        self.init_sequences = {}
        self.position_weights = {}
        if sequences is None:
            for z in list(zip(unique_protomers, lengths)):
                self.init_sequences[z[0]] = ''.join(np.random.choice(list(AA_freq.keys()), size=z[1], p=list(AA_freq.values())))
                self.position_weights[z[0]] = np.ones(z[1]) / z[1]

        else:
            for p, proto in enumerate(unique_protomers):
                if sequences[p] == '': # empty sequence
                    self.init_sequences[proto] = ''.join(np.random.choice(list(AA_freq.keys()), size=lengths[p], p=list(AA_freq.values())))
                    self.position_weights[proto] = np.ones(lengths[p]) / lengths[p]

                else:
                    self.init_sequences[proto] = sequences[p]
                    if position_weights is None:
                        self.position_weights[proto] = np.ones(lengths[p]) / lengths[p]

                    else:
                        if position_weights[p] == '':
                            self.position_weights[proto] = np.ones(lengths[p]) / lengths[p]

                        else:
                            self.position_weights[proto] = np.array(position_weights[p])

        # Initialise lengths.
        self.lengths = {}
        for proto, seq in self.init_sequences.items():
            self.lengths[proto] = len(seq)

        self.current_sequences = {p:s for p, s in self.init_sequences.items()}
        self.try_sequences = {p:s for p, s in self.init_sequences.items()}

    # Method functions.
    def assign_mutable_positions(self, mutable_positions):
        '''Assign dictonary of protomers with arrays of mutable positions.'''
        self.mutable_positions = mutable_positions

    def assign_mutations(self, mutated_protomers):
        '''Assign mutated sequences to try_sequences.'''
        self.try_sequences = mutated_protomers

    def update_mutations(self):
        '''Update current sequences to try sequences.'''
        self.current_sequences = copy.deepcopy(self.try_sequences)


class Oligomer:
    '''
    Object for keeping track of oligomers during hallucination.
    Also keeps track of AlphaFold2 scores and structure.
    Sort of like a Pose object in Rosetta.
    '''

    def __init__(self, oligo_string:str, protomers:Protomers):

        self.name = oligo_string
        self.positive_design = (lambda x: True if x=='+' else False)(oligo_string[-1])
        self.subunits = oligo_string[:-1]

        # Initialise oligomer sequence (concatanation of protomers).
        self.init_seq = ''
        for unit in self.subunits:
            self.init_seq += protomers.init_sequences[unit]

        # Initialise overall length and protomer lengths.
        self.oligo_L = len(self.init_seq)

        self.chain_Ls = []
        for unit in self.subunits:
            self.chain_Ls.append(len(protomers.init_sequences[unit]))

        self.current_seq = str(self.init_seq)
        self.try_seq = str(self.init_seq)

    def init_prediction(self, af2_prediction):
        '''Initalise scores/structure'''
        self.init_prediction_results, self.init_unrelaxed_structure = af2_prediction
        self.current_prediction_results = copy.deepcopy(self.init_prediction_results)
        self.current_unrelaxed_structure = copy.deepcopy(self.init_unrelaxed_structure)
        self.try_prediction_results = copy.deepcopy(self.init_prediction_results)
        self.try_unrelaxed_structure = copy.deepcopy(self.init_unrelaxed_structure)

    def init_loss(self, loss):
        '''Initalise loss'''
        self.init_loss = loss
        self.current_loss = float(self.init_loss)
        self.try_loss = float(self.init_loss)

    def assign_oligo(self, protomers):
        '''Make try oligomer sequence from protomer sequences'''
        self.try_seq = ''
        for unit in self.subunits:
            self.try_seq += protomers.try_sequences[unit]

    def update_oligo(self):
        '''Update current oligomer sequence to try ones.'''
        self.current_seq = str(self.try_seq)

    def assign_prediction(self, af2_prediction):
        '''Assign try AlphaFold2 prediction (scores and structure).'''
        self.try_prediction_results, self.try_unrelaxed_structure = af2_prediction

    def update_prediction(self):
        '''Update current scores/structure to try scores/structure.'''
        self.current_unrelaxed_structure = copy.deepcopy(self.try_unrelaxed_structure)
        self.current_prediction_results = copy.deepcopy(self.try_prediction_results)

    def assign_loss(self, loss):
        '''Assign try loss.'''
        self.try_loss = float(loss)

    def update_loss(self):
        '''Update current loss to try loss.'''
        self.current_loss = float(self.try_loss)


##################################
# INITALISATION
##################################

args = get_args(); print('#', args)

os.makedirs(f'{args.out}_models', exist_ok=True) # where all the outputs will go.

# Notes.
print(f'> Git commit: {args.commit}')
if args.single_chain == True:
    print('> Design(s) will be generated as sequence-symmetric repeat proteins, NOT oligomers.')
    print('> The following repeat proteins will be designed:')
else:
    print(f'> The following oligomers will be designed:')
for i, oligo in enumerate(args.oligo.strip(',').split(',')):
    print(f' >> {oligo[:-1]} ({(lambda x: "positive" if x=="+" else "negative")(oligo[-1])} design), contributing {args.oligo_weights[i]} to the global loss')
print(f'> Simulated annealing will be performed over {args.steps} steps with a starting temperature of {args.T_init} and a half-life for the temperature decay of {args.half_life} steps.')
print(f'> The mutation rate at each step will go from {args.mutation_rate.split("-")[0]} to {args.mutation_rate.split("-")[1]} over {args.steps} steps (stepped linear decay).')
if args.tolerance is not None:
    print(f'> A tolerance setting of {args.tolerance} was set, which might terminate the MCMC trajectory early.')
print(f'> The choice of position to mutate at each step will be based on {args.select_positions}, with parameter(s): {args.select_position_params}.')
print(f'> At each step, selected positions will be mutated based on {args.mutation_method}.')
print(f'> Predictions will be performed with AlphaFold2 model_{args.model}_ptm, with recyling set to {args.recycles}, and {args.msa_clusters} MSA cluster(s).')
print(f'> The loss function used during optimisation was set to: {args.loss}, with respective weights: {args.loss_weights}.')

# Amino-acid frequencies taken from background frequencies of BLOSUM62.
# Data from https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/composition_adjustment/matrix_frequency_data.c
AA_freq = {'A': 0.07421620506799341,
 'R': 0.05161448614128464,
 'N': 0.044645808512757915,
 'D': 0.05362600083855441,
 'C': 0.02468745716794485,
 'Q': 0.03425965059141602,
 'E': 0.0543119256845875,
 'G': 0.074146941452645,
 'H': 0.026212984805266227,
 'I': 0.06791736761895376,
 'L': 0.09890786849715096,
 'K': 0.05815568230307968,
 'M': 0.02499019757964311,
 'F': 0.04741845974228475,
 'P': 0.038538003320306206,
 'S': 0.05722902947649442,
 'T': 0.05089136455028703,
 'W': 0.013029956129972148,
 'Y': 0.03228151231375858,
 'V': 0.07291909820561925}

for aa in args.exclude_AA:
    del AA_freq[aa]

# Re-compute frequencies to sum to 1.
sum_freq = np.sum(list(AA_freq.values()))
adj_freq = [f/sum_freq for f in list(AA_freq.values())]
AA_freq = dict(zip(AA_freq, adj_freq))

print(f'> Allowed amino acids: {len(AA_freq.keys())} [{" ".join([aa for aa in list(AA_freq.keys())])}]')
print(f'> Excluded amino acids: {len(args.exclude_AA)} [{" ".join([aa for aa in args.exclude_AA])}]')

# Initialise Protomer object (one for the whole simulation).
if args.proto_sequences is None:
    protomers = Protomers(unique_protomers=args.unique_protomers, lengths=args.proto_Ls, aa_freq=AA_freq)

else:
    protomers = Protomers(unique_protomers=args.unique_protomers, lengths=args.proto_Ls, aa_freq=AA_freq, sequences=args.proto_sequences, position_weights=args.position_weights)

for proto, seq in protomers.init_sequences.items():
    print(f' >> Protomer {proto} init sequence: {seq}')
    print(f' >> Protomer {proto} position-specific weights: {protomers.position_weights[proto]}')

# Initialise Oligomer objects (one for each specified oligomer).
oligomers = {}
for o in args.oligo.split(','):
    oligomers[o] = Oligomer(o, protomers)

# Setup AlphaFold2 models.
model_runners = setup_models(args.oligo.split(','), model_id=args.model, recycles=args.recycles, msa_clusters=args.msa_clusters)

# Start score file.
with open(f'{args.out}_models/{os.path.splitext(os.path.basename(args.out))[0]}.out', 'w') as f:
    print_str = f'# {args}\n'
    print_str += 'step accepted temperature mutations loss plddt ptm pae'
    for oligo in oligomers.keys():
        print_str += f' sequence_{oligo} loss_{oligo} plddt_{oligo} ptm_{oligo} pae_{oligo}'
    print_str += '\n'
    f.write(print_str)


####################################
# MCMC WITH SIMULATED ANNEALING
####################################

Mi, Mf = args.mutation_rate.split('-')
M = np.linspace(int(Mi), int(Mf), args.steps) # stepped linear decay of the mutation rate

current_loss = np.inf
rolling_window = []
rolling_window_width = 100
for i in range(args.steps):

    if args.tolerance is not None and i > rolling_window_width: # check if change in loss falls under the tolerance threshold for terminating the simulation.

        if np.std(rolling_window[-rolling_window_width:]) < args.tolerance:
            print(f'The change in loss over the last 100 steps has fallen under the tolerance threshold ({args.tolerance}). Terminating the simulation...')
            sys.exit()
    else:

        # Update a few things.
        T = args.T_init * (np.exp(np.log(0.5) / args.half_life) ** i) # update temperature
        n_mutations = round(M[i]) # update mutation rate
        accepted = False # reset
        try_losses = []

        if i == 0: # do a first pass through the network before mutating anything -- baseline
            print('-' * 100)
            print('Starting...')
            for name, oligo in oligomers.items():
                af2_prediction = predict_structure(oligo,
                                                    args.single_chain,
                                                    model_runners[name],
                                                    random_seed=np.random.randint(42)) # run AlphaFold2 prediction
                oligo.init_prediction(af2_prediction) # assign
                loss = compute_loss(args.loss,
                                    oligo,
                                    args,
                                    args.loss_weights) # calculate the loss
                oligo.init_loss(loss) # assign
                try_losses.append(loss) # increment global loss

        else:

            # Mutate protomer sequences and generate updated oligomer sequences
            protomers.assign_mutable_positions(select_positions(n_mutations,
                                                                protomers,
                                                                oligomers,
                                                                args.select_positions,
                                                                args.select_position_params)) # define mutable positions for each protomer
            protomers.assign_mutations(mutate(args.mutation_method,
                                                protomers,
                                                AA_freq)) # mutate those positions

            for name, oligo in oligomers.items():
                oligo.assign_oligo(protomers) # make new oligomers from mutated protomer sequences
                oligo.assign_prediction(predict_structure(oligo,
                                                            args.single_chain,
                                                            model_runners[name],
                                                            random_seed=np.random.randint(42))) # run AlphaFold2 prediction
                loss = compute_loss(args.loss, oligo, args, args.loss_weights) # calculate the loss for that oligomer
                oligo.assign_loss(loss) # assign the loss to the object (for tracking)
                try_losses.append(loss) # increment the global loss

        # Normalize oligo weights vector.
        oligo_weights_normalized = np.array(args.oligo_weights) / np.sum(args.oligo_weights)

        # Global loss is the weighted average of the individual oligomer losses.
        try_loss = np.mean( np.array(try_losses) * oligo_weights_normalized )

        delta = try_loss - current_loss # all losses must be defined such that optimising equates to minimising.

        # If the new solution is better, accept it.
        if delta < 0:
            accepted = True

            print(f'Step {i:05d}: change accepted >> LOSS {current_loss:2.3f} --> {try_loss:2.3f}')

            current_loss = float(try_loss) # accept loss change
            protomers.update_mutations() # accept sequence changes

            for name, oligo in oligomers.items():
                print(f' > {name} loss  {oligo.current_loss:2.3f} --> {oligo.try_loss:2.3f}')
                print(f' > {name} plddt {np.mean(oligo.current_prediction_results["plddt"]):2.3f} --> {np.mean(oligo.try_prediction_results["plddt"]):2.3f}')
                print(f' > {name} ptm   {oligo.current_prediction_results["ptm"]:2.3f} --> {oligo.try_prediction_results["ptm"]:2.3f}')
                print(f' > {name} pae   {np.mean(oligo.current_prediction_results["predicted_aligned_error"]):2.3f} --> {np.mean(oligo.try_prediction_results["predicted_aligned_error"]):2.3f}')
                oligo.update_oligo() # accept sequence changes
                oligo.update_prediction() # accept score/structure changes
                oligo.update_loss() # accept loss change

            print('=' * 70)

        # If the new solution is not better, accept it with a probability of e^(-cost/temp).
        else:

            if np.random.uniform(0, 1) < np.exp( -delta / T):
                accepted = True

                print(f'Step {i:05d}: change accepted despite not improving the loss >> LOSS {current_loss:2.3f} --> {try_loss:2.3f}')

                current_loss = float(try_loss)
                protomers.update_mutations() # accept sequence changes

                for name, oligo in oligomers.items():
                    print(f' > {name} loss  {oligo.current_loss:2.3f} --> {oligo.try_loss:2.3f}')
                    print(f' > {name} plddt {np.mean(oligo.current_prediction_results["plddt"]):2.3f} --> {np.mean(oligo.try_prediction_results["plddt"]):2.3f}')
                    print(f' > {name} ptm   {oligo.current_prediction_results["ptm"]:2.3f} --> {oligo.try_prediction_results["ptm"]:2.3f}')
                    print(f' > {name} pae   {np.mean(oligo.current_prediction_results["predicted_aligned_error"]):2.3f} --> {np.mean(oligo.try_prediction_results["predicted_aligned_error"]):2.3f}')
                    oligo.update_oligo() # accept sequence changes
                    oligo.update_prediction() # accept score/structure changes
                    oligo.update_loss() # accept loss change

                print('=' * 70)

            else:
                accepted = False
                print(f'Step {i:05d}: change rejected >> LOSS {current_loss:2.3f} !-> {try_loss:2.3f}')
                print('-' * 70)

        sys.stdout.flush()


        # Save PDB if move was accepted.
        if accepted == True:

            for name, oligo in oligomers.items():

                with open(f'{args.out}_models/{os.path.splitext(os.path.basename(args.out))[0]}_{oligo.name}_step_{str(i).zfill(5)}.pdb', 'w') as f:
                    # write pdb
                    if args.amber_relax == 0 :
                        pdb_lines = protein.to_pdb(oligo.current_unrelaxed_structure).split('\n')
                    elif args.amber_relax == 1 :
                        pdb_lines = amber_relax(oligo.current_unrelaxed_structure).split('\n')

                    # Identify chain breaks and re-assign chains correctly before generating PDB file.
                    split_lines = [l.split() for l in pdb_lines if 'ATOM' in l]
                    split_lines = np.array([l[:4] + [l[4][0]] + [l[4][1:]] + l[5:] if len(l)<12 else l for l in split_lines]) # chain and resid no longer space-separated at high resid.
                    splits = np.argwhere(np.diff(split_lines.T[5].astype(int))>1).flatten() + 1 # identify idx of chain breaks based on resid jump.
                    splits = np.append(splits, len(split_lines))
                    chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 
                            'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN', 'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY', 'ZZ', 
                            'AAA', 'BBB', 'CCC', 'DDD', 'EEE', 'FFF', 'GGG', 'HHH', 'III', 'JJJ', 'KKK', 'LLL', 'MMM', 'NNN', 'OOO', 'PPP', 'QQQ', 'RRR', 'SSS', 'TTT', 'UUU', 'VVV', 'WWW', 'XXX', 'YYY', 'ZZZ']
                    chain_str = ''
                    prev = 0
                    for ch, resid in enumerate(splits): # make new chain string
                        length = resid - prev
                        chain_str += chains[ch] * length
                        prev = resid
                    atom_lines = [l for l in pdb_lines if 'ATOM' in l]
                    new_lines = [l[:21]+chain_str[k]+l[22:] for k, l in enumerate(atom_lines) if 'ATOM' in l] # generate chain-corrected PDB lines.

                    # write PDB file and append scores at the end of it.
                    f.write('MODEL     1\n')
                    f.write('\n'.join(new_lines))
                    f.write('\nENDMDL\nEND\n')
                    f.write(f'plddt_array {",".join(oligo.current_prediction_results["plddt"].astype(str))}\n')
                    f.write(f'plddt {np.mean(oligo.current_prediction_results["plddt"])}\n')
                    f.write(f'ptm {oligo.current_prediction_results["ptm"]}\n')
                    f.write(f'pae {np.mean(oligo.current_prediction_results["predicted_aligned_error"])}\n')
                    f.write(f'loss {oligo.current_loss}\n')
                    f.write(f'# {str(args)}\n')

                # Optionally save the PAE matrix
                if args.output_pae == True:
                    np.save(f'{args.out}_models/{os.path.splitext(os.path.basename(args.out))[0]}_{oligo.name}_step_{str(i).zfill(5)}.npy', oligo.current_prediction_results['predicted_aligned_error'])


        # Save scores for the step (even if rejected).
        # step accepted temperature mutations loss plddt ptm pae '
        score_string = f'{i:05d} '
        score_string += f'{accepted} '
        score_string += f'{T} '
        score_string += f'{n_mutations} '
        score_string += f'{try_loss} '
        score_string += f'{np.mean([np.mean(r.try_prediction_results["plddt"]) for r in oligomers.values()])} '
        score_string += f'{np.mean([r.try_prediction_results["ptm"] for r in oligomers.values()])} '
        score_string += f'{np.mean([np.mean(r.try_prediction_results["predicted_aligned_error"]) for r in oligomers.values()])} '

        for name, oligo in oligomers.items():
            breaked_seq = ''
            Lprev = 0
            for L in oligo.chain_Ls:
                Lcorr = Lprev + L
                breaked_seq += oligo.try_seq[Lprev:Lcorr] + '/'
                Lprev = Lcorr

            score_string += f'{breaked_seq[:-1]} '
            score_string += f'{oligo.try_loss} '
            score_string += f'{np.mean(oligo.try_prediction_results["plddt"])} '
            score_string += f'{oligo.try_prediction_results["ptm"]} '
            score_string += f'{np.mean(oligo.try_prediction_results["predicted_aligned_error"])} '

        with open(f'{args.out}_models/{os.path.splitext(os.path.basename(args.out))[0]}.out', 'a') as f:
            f.write(score_string + '\n')

        rolling_window.append(current_loss)

print('Done')
