# arg parser module

import argparse
import os, sys, datetime
import numpy as np
import subprocess
script_dir = os.path.dirname(os.path.realpath(__file__))

def get_args():
    ''' Parse input arguments'''

    parser = argparse.ArgumentParser(
            description='AlphaFold2 multi-state hallucination using a Markov chain Monte Carlo protocol. Monomers, homo-oligomers and hetero-oligomers can be hallucinated. Both positive and negative multistate designs are possible.'
            )

    # General arguments.
    parser.add_argument(
            '--oligo',
            help='oligomer(s) definitions (comma-separated string, no space). Numbers and types of subunits (protomers), and design type (positive or negative) specifying each oligomer.\
             Protomers are defined by unique letters, and strings indicate oligomeric compositions. The last character of each oligomer has to be [+] or [-] to indicate positive or negative design of that oligomer (e.g. AAAA+,AB+).\
             The number of unique protomers must match --L / --seq. Must be specified.',
            action='store',
            type=str,
            required=True
            )

    parser.add_argument(
            '--L',
            help='lengths of each protomer (comma-separated, no space). Must be specified if not using --seq or a .af2h config file.',
            action='store',
            type=str
            )

    parser.add_argument(
            '--seq',
            help='seed sequence for each protomer (comma-separated, no space). Optional.',
            action='store',
            type=str
            )

    parser.add_argument(
            '--out',
            help='the prefix appended to the output files. Must be specified.',
            action='store',
            type=str,
            required=True
            )

    parser.add_argument(
            '--single_chain',
            default=False,
            action='store_true',
            help='this option will generate sequence-symmetric repeat proteins instead of oligomers by removing chain breaks (default: False).',
            )

    # Hallucination arguments.
    parser.add_argument(
            '--exclude_AA',
            default='C',
            action='store',
            type=str,
            help='amino acids to exclude during hallucination. Must be a continous string (no spaces) in one-letter code format (default: %(default)s).'
            )

    parser.add_argument(
            '--mutation_rate',
            default='3-1',
            action='store',
            help='number of mutations at each MCMC step (start-finish, stepped linear decay). Should probably be scaled with protomer length (default: %(default)s).'
            )

    parser.add_argument(
            '--select_positions',
            default='random',
            action='store',
            type=str,
            help='how to select positions for mutation at each step. Choose from [random, plddt::quantile, FILE.af2h::quantile].\
            FILE.af2h needs to be a file specifying the probability of mutation at each site.\
            Optional arguments can be given with :: e.g. plddt::0.25 will only mutate the 25%% lowest plddt positions (default: %(default)s).'
            )

    parser.add_argument(
            '--mutation_method',
            default='frequency_adjusted',
            action='store',
            type=str,
            help='how to mutate selected positions. Choose from [uniform, frequency_adjusted, blosum62, pssm] (default: %(default)s).'
            )

    parser.add_argument(
            '--loss',
            default='dual',
            type=str,
            help='the loss function used during optimization. Choose from \
            [plddt, ptm, pae, dual, cyclic, dual_cyclic, pae_sub_mat, pae_asym, tmalign (requires --template), dual_tmalign (requires --template), \
            aspect_ratio, frac_dssp, min_frac_dssp (requires --dssp_fractions_specified), pae_asym_tmalign (in development), entropy (in development)].\
            Multiple losses can be combined as a comma-separarted string of loss_name:args units (and weighed with --loss_weights).\
            loss_0_name::loss0_param0;loss0_param1,loss_1_name::[loss_1_configfile.conf] ... \
             (default: %(default)s).'
            )

    parser.add_argument(
            '--loss_weights',
            default=None,
            type=str,
            action='store',
            help='if a combination of losses is passed, specify relative weights of each loss to the globabl loss by providing a comma-separated list of relative weights. \
            E.g. 2,1 will make the first loss count double relative to the second one (default: equal weights).'
            )

    parser.add_argument(
            '--oligo_weights',
            default=None,
            type=str,
            help='contribution of the loss of each oligomer to the global loss, provided as a comma-separted list of relative weights (default: equal weights).'
            )

    # MCMC arguments.
    parser.add_argument(
            '--T_init',
            default=0.01,
            action='store',
            type=float,
            help='starting temperature for simulated annealing. Temperature is decayed exponentially (default: %(default)s).'
            )

    parser.add_argument(
            '--half_life',
            default=1000,
            action='store',
            type=float,
            help='half-life for the temperature decay during simulated annealing (default: %(default)s).'
            )

    parser.add_argument(
            '--steps',
            default=5000,
            action='store',
            type=int,
            help='number for steps for the MCMC trajectory (default: %(default)s).'
            )

    parser.add_argument(
            '--tolerance',
            default=None,
            action='store',
            type=float,
            help='the tolerance on the loss sliding window for terminating the MCMC trajectory early (default: %(default)s).'
            )

    # AlphaFold2 arguments.
    parser.add_argument(
            '--model',
            default=4,
            action='store',
            type=int,
            help='AF2 model (_ptm) used during prediction. Choose from [1, 2, 3, 4, 5] (default: %(default)s).'
            )

    parser.add_argument(
            '--amber_relax',
            default=0,
            action='store',
            type=int,
            help='amber relax pdbs written to disk, 0=do not relax, 1=relax every prediction (default: %(default)s).'
            )

    parser.add_argument(
            '--recycles',
            default=1,
            action='store',
            type=int,
            help='the number of recycles through the network used during structure prediction. Larger numbers increase accuracy but linearly affect runtime (default: %(default)s).'
            )

    parser.add_argument(
            '--msa_clusters',
            default=1,
            action='store',
            type=int,
            help='the number of MSA clusters used during feature generation (?). Larger numbers increase accuracy but significantly affect runtime (default: %(default)s).'
            )

    parser.add_argument(
            '--output_pae',
            default=False,
            action='store_true',
            help='output the pAE (predicted alignment error) matrix for each accepted step of the MCMC trajectory (default: %(default)s).'
            )

    parser.add_argument(
            '--timestamp',
            default=False,
            action='store_true',
            help='timestamp output and every PDB written to disk with: %%Y%%m%%d_%%H%%M%%S_%%f (default: %(default)s).'
            )

    parser.add_argument(
            '--template',
            default=None,
            type=str,
            action='store',
            help='template PDB for use with tmalign-based losses (default: %(default)s).'
            )

    parser.add_argument(
            '--dssp_fractions_specified',
            default=None,
            type=str,
            action='store',
            help='dssp fractions specfied for frac_dssp loss as E(beta sheet), H(alpha helix), notEH(other)\
            e.g. 0.8,None,None will enforce 80%% beta sheet; or 0.5,0,None will enforce 50%% beta sheet, no helices (default: %(default)s).'
            )

    parser.add_argument(
            '--template_alignment',
            default=None,
            type=str,
            action='store',
            help='enforce tmalign alignment with fasta file (default: %(default)s).'
            )

    args = parser.parse_args()


    ########################################
    # SANITY CHECKS
    ########################################

    # Errors.
    if args.oligo is None:
        print('ERROR: the definiton for each oligomer must be specified. System exiting...')
        sys.exit()

    if args.L is None and args.seq is None and '.af2h' not in args.select_positions:
        print('ERROR: either seed sequence(s) or length(s) must specified. System exiting...')
        sys.exit()

    if np.any([(lambda x: True if x not in ['+', '-'] else False)(d[-1]) for d in args.oligo.split(',')]):
        print('ERROR: the type of design (positive [+] or negative [-]) must be specified for each oligomer. System existing...')
        sys.exit()

    if "tmalign" in args.loss and args.template is None:
        print('ERROR: tmalign loss require a --template [PDBFILE]. System exiting...')
        sys.exit()
    if "frac_dssp" in args.loss and args.dssp_fractions_specified is None:
        print('ERROR: frac_dssp loss requires a --dssp_fractions_specified [E,H,notEH]. System exiting...')
        sys.exit()


    # Warnings.
    if (args.L is not None) and (args.seq is not None or '.af2h' in args.select_positions):
        print('WARNING: Both user-defined sequence(s) and length(s) were provided. Are you sure of what you are doing? The simulation will continue assuming you wanted to use the provided sequence(s) as seed(s).')

    ##########################################
    # ARGUMENT MODIFICATIONS BASED ON INPUTS
    ##########################################

    # GENERAL
    if args.timestamp == True:
        date_time_str = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        args.out = args.out + "_" + date_time_str

    # Add some arguments.
    args.commit = "" #subprocess.check_output(f'git --git-dir {script_dir}/../.git rev-parse HEAD', shell=True).decode().strip() # add git hash of current commit.

    args.unique_protomers = sorted(set(args.oligo.replace(',','').replace('+','-').replace('-','')))


    # LOSSSES
    # all losses stored in list of lists [ [ loss_name, [loss_param0, loss_param1] ], ...]
    # losses processed in order
    # if no parameters necessary list entry is empty list
    losses = []
    for curr_loss_str in args.loss.strip(',').split(','):
        loss_parameters = []

        if '::' in curr_loss_str:
            loss_name, loss_arguments = curr_loss_str.split('::')

            for curr_loss_param in loss_arguments.strip(';').split(';'):
                if "[" in curr_loss_param:
                    print("loss configfile NOT IMPLEMENTED YET. System exiting...")
                    sys.exit()
                else:
                    #assuming all parameters are number
                    loss_parameters.append( float(curr_loss_param) )
        else:
            loss_name = str(curr_loss_str)

        losses.append([loss_name, loss_parameters])

    # replace args.loss string with new dictionary
    args.loss = losses

    # loss weights only relevant if more than one loss declared
    if args.loss_weights != None:
        # split relative weight for losses from input string,
        # update as list
        loss_weights = []
        for curr_loss_weight in args.loss_weights.strip(',').split(','):
            if "-" in curr_loss_weight :
                # loss is ramped over course of run
                print("loss weight ramping NOT IMPLMENTED YET. System exiting...")
                sys.exit()
                loss_weights.append( [ float(i) for i in curr_loss_weight.split("-") ] )
                pass
            else:
                #loss is constant over course of run (ramped to idential values)
                # loss_weights.append( [ float(curr_loss), float(curr_loss) ] )
                loss_weights.append( float(curr_loss_weight) )

        args.loss_weights = loss_weights

        assert len(args.loss_weights) == len(args.loss)
    else:
        # no loss weights declared, all losses weighed equally
        args.loss_weights = [1.0] * len(args.loss) # processed list of losses

    if args.oligo_weights != None:
        oligo_weights = []
        for curr_oligo_weight in args.oligo_weights.strip(',').split(','):
            oligo_weights.append(float(curr_oligo_weight))

        args.oligo_weights = oligo_weights

    else:
        # no oligo weights declared, all oligos weighed equally.
        args.oligo_weights = [1.0] * len(args.oligo.strip(',').split(','))


    # SEQUENCES
    # process before reading in the update args of .af2h files.
    # intialise empty sequences in case none are given.
    # initalise position_weights for downward compatibility.
    args.proto_sequences = None
    args.position_weights = None
    if args.seq != None:
        args.proto_sequences = args.seq.split(',')

    elif '.af2h' in args.select_positions:
        # Parse .af2h file -- should be fasta, with an extra line after the sequence.
        # The line after the sequence should be a comma-separated list of values (of the same length as the sequence) that represents the probability of mutating each position.
        with open(args.select_positions.split('::')[0] , 'r') as f:
            lines = list(line for line in (l.strip() for l in f) if line) # strip empty lines.
            seq_prob = {}
            for entry in np.reshape(lines, (-1, 3)):
                freq = np.array(entry[2].split(','), dtype=float)
                if freq.sum() == 0: # if no mutations are desired
                    ajd_freq = freq
                else:
                    ajd_freq =  freq / freq.sum() # re-adjust frequencies to sum to 1 across the length of each protomer.
                #HAAAAAAACKY
                if "X" in entry[1]: # replate X with random aa:
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
                    sum_freq = np.sum(list(AA_freq.values()))
                    adj_freq = [f/sum_freq for f in list(AA_freq.values())]
                    AA_freq = dict(zip(AA_freq, adj_freq))
                    
                    new_sequence = ""
                    print("replace X with rand aa: ", entry[1])
                    for this_aa in entry[1]:
                        if this_aa == "X":
                            new_sequence += np.random.choice(list(AA_freq.keys()), 
                                                p=list(AA_freq.values() ) )
                        else:
                            new_sequence += this_aa
                    entry[1] = new_sequence
                    print(entry[1])
                # HAAAAACKy end

                seq_prob[entry[0][1:]] = {'seq':entry[1], 'prob':ajd_freq}

            for proto in args.unique_protomers: # complete with empty entries in case the user did not specify all protomers -- these will be replace by randomly generated sequences at initalisation.
                if proto not in list(seq_prob.keys()):
                    seq_prob[proto] = {'seq':'', 'prob':''}

        args.proto_sequences  = [seq_prob[proto]['seq'] for proto in args.unique_protomers]
        args.position_weights = [list(seq_prob[proto]['prob']) for proto in args.unique_protomers]

    if args.L is None:
        args.proto_Ls = [len(seq) for seq in args.proto_sequences]

    else:
        args.proto_Ls = [int(length) if length!='' else 0 for length in args.L.split(',')]


    # UPDATE / MUTATIONS
    # reading in additional arguments for updates, currently just option for
    # quantile to mutate when plddt or .af2h specified
    args.select_position_params = None
    if '::' in args.select_positions:
        # additional arguments in args.update
        args.select_positions, select_position_params = args.select_positions.split('::')
        args.select_position_params = float(select_position_params)
        print (" Select position params set to ", args.select_position_params)

    # Additional Errors.
    if args.seq is None:
        if len(args.unique_protomers) != len(args.proto_Ls):
            print('ERROR: the number of unique protomers and the number of specified lengths must match. System exiting...')
            sys.exit()

    if args.L is None:
        if len(args.unique_protomers) != len(args.proto_Ls):
             print('ERROR: the number of unique protomers and the number of specified sequences must match. System exiting...')
             sys.exit()


    if args.dssp_fractions_specified is not None:
        dssp_assignments = ["E", "H", "L"]
        dssp_dict = {}
        for sec_struc, entry in zip( dssp_assignments , args.dssp_fractions_specified.strip(',').split(',') ):
            print(sec_struc, entry)
            if entry == "None" or entry == None or entry == "none" :
                dssp_dict[sec_struc] = None
            else: 
                dssp_dict[sec_struc] = float(entry)

        args.dssp_fractions_specified = dssp_dict
        print(dssp_dict)

    return args
