# sequence mutation module

import numpy as np
import sys


##########################
# GENERAL FUNCTIONS
##########################

def extract_stacked_plddts(proto_object, oligo_dict):
    '''
    Returns dictonary of protomers with associated stacked plddt array (one array from each oligomer context).
    Used for plddt-based update methods.
    '''

    unique_protomers = list(proto_object.init_sequences.keys())

    # For each protomer, get the plddt arrays of its different states/oligomers.
    stacked_plddts = {proto:[] for proto in unique_protomers}
    for oligo in oligo_dict.values():

        prev = 0
        ch_ranges = []
        for L in oligo.chain_Ls:
            Lcorr = prev + L
            ch_ranges.append([prev, Lcorr])
            prev = Lcorr

        for n, proto in enumerate(oligo.subunits):
            stacked_plddts[proto].append(oligo.current_prediction_results['plddt'][ch_ranges[n][0]:ch_ranges[n][1]])

    return stacked_plddts


##############################
# SELECT POSITIONS
##############################

def select_positions(n_mutations, proto_object, oligo_dict, select_positions, select_position_params):
    '''
    Select mutable positions in each protomer based on a specific method.
    Returns a dictionary of protomers with associated arrays indicating mutable positions.
    '''

    mutable_positions = {}

    if select_positions == 'random':
        # Choose positions randomly.
        for proto, seq in proto_object.current_sequences.items():
            mutable_positions[proto] = np.random.choice(range(len(seq)), size=n_mutations, replace=False)


    elif select_positions == 'plddt':
        # Choose positions based on lowest plddt (taken across all states/oligomer for each protomer).
        # First/last three positions of each protomers are choice frequency adjusted to avoid picking N/C term every time (they tend to score much lower).

        if select_position_params != None:
            mutate_plddt_quantile = float(select_position_params) # manually specified quantile.
        else:
            mutate_plddt_quantile = 0.25 # default worst pLDDT quantile to mutate.

        stacked_plddts = extract_stacked_plddts(proto_object, oligo_dict)

        for proto, plddts in stacked_plddts.items():

            proto_L = proto_object.lengths[proto]

            # Weights associated with each position in the protomer.
            # to account for termini systematically scoring worse in pLDDT.
            weights = np.array([0.25, 0.5, 0.75] + [1] * (proto_L - 6) + [0.75, 0.5, 0.25])

            # Sub-select lowest % quantile of plddt positions.
            n_potential = round(proto_L * mutate_plddt_quantile)
            consensus_min = np.min(plddts, axis=0)
            potential_sites = np.argsort(consensus_min)[:n_potential]

            # Select mutable sites
            sub_w = weights[potential_sites]
            sub_w = [w/np.sum(sub_w) for w in sub_w]
            sites = np.random.choice(potential_sites, size=n_mutations, replace=False, p=sub_w)

            mutable_positions[proto] = sites


    elif '.af2h' in select_positions:
        # Choice of protomer positions to mutate is based on weights provided in .af2h file.
        # If plddt parameter is given (as .af2h::0.xx), select worst pLDDT scoring quantile among user-specified positions.
        # Otherwise choose positions among user-defined positions with provided frequency-adjustment.

        # if selecting worst residues based on plddt, extract stacked pLDDT array across all protomers in the various oligomers.
        if select_position_params != None :
            quantile = float(select_position_params) # user-specified quantile.
            stacked_plddts = extract_stacked_plddts(proto_object, oligo_dict)

        for proto, seq in proto_object.current_sequences.items():

            position_weights = proto_object.position_weights[proto]
            if np.sum(position_weights) == np.nan or np.sum(position_weights) == 0.0 :
                mutable_positions[proto] = []
            else:
                if select_position_params != None :
                    # position weights is a np array
                    # make position weights 1.0 if non-zero for this protomer
                    # np.nonzero returns tuple of arrays along each axis, here choose for axis=0
                    allowed_positions = np.nonzero( np.ceil( proto_object.position_weights[proto] ) )[0]
                    # extract stacked plddts for this protomer from stacked dict
                    plddts = stacked_plddts[proto]
                    # Sub-select lowest % quantile in pLDDT of allowed positions
                    n_potential = round(len(allowed_positions) * quantile)
                    # find absolute worst plddt containing protomer for given positions
                    # !!! needs review
                    consensus_min = np.min([a[allowed_positions] for a in plddts], axis=0)
                    # sort extract plddts for positions, sort, find worst quantile, extract indices of these
                    # from array of allowed positions
                    potential_sites = allowed_positions[np.argsort(consensus_min)[:n_potential]]
                    #randomly select indices to mutate within the worst quantile
                    mutable_positions[proto] = potential_sites[ np.random.choice(range(len(potential_sites)), size=n_mutations, replace=False) ]

                else:
                    position_weights = proto_object.position_weights[proto]
                    # Choice of positions is biased by the user-defined position-specific weights.
                    mutable_positions[proto] = np.random.choice(range(len(seq)), size=n_mutations, replace=False, p=position_weights)



    return mutable_positions

################################
# MUTATE POSITIONS
################################

def mutate(method, proto_object, aa_freq):
    '''Mutate mutable positions based on a specific method.'''

    mutated_proto_sequences = {}

    if method == 'uniform':
        # Mutate positions with uniform probabilities.
        for proto, seq in proto_object.current_sequences.items():

            for p in proto_object.mutable_positions[proto]:
                        seq = seq[:p] + np.random.choice(list(aa_freq.keys())) + seq[p+1:]

            mutated_proto_sequences[proto] = seq


    elif method == 'frequency_adjusted':
        # Mutate positions based on frequency-ajusted AA substitutions.
        for proto, seq in proto_object.current_sequences.items():

            for p in proto_object.mutable_positions[proto]:
                        seq = seq[:p] + np.random.choice(list(aa_freq.keys()), p=list(aa_freq.values())) + seq[p+1:]

            mutated_proto_sequences[proto] = seq


    elif method == 'blosum62':
        # Mutate positions based on BLOSUM62 substituation matrix.
        # BLOSUM62 probabilies taken from:
        # https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/composition_adjustment/matrix_frequency_data.c

        # Joint probabilities for BLOSUM62.
        joint_p = [
          [2.1497573378347484e-02, 2.3470224274721213e-03, 1.9493235258876179e-03,
           2.1674844853066858e-03, 1.5903351423026848e-03, 1.9242657898716525e-03,
           2.9879059292799641e-03, 5.8158526388051033e-03, 1.1076584657559144e-03,
           3.1880644746334580e-03, 4.4186245468471547e-03, 3.3466571942021082e-03,
           1.3412107617355408e-03, 1.6360627863999076e-03, 2.1568959784943114e-03,
           6.2524987419815400e-03, 3.7180506975672363e-03, 4.0281679108936688e-04,
           1.2999956675626666e-03, 5.0679056444508912e-03],
          [2.3470224274721213e-03, 1.7757465118386322e-02, 1.9786027128591904e-03,
           1.5865480081162602e-03, 3.9365984789376245e-04, 2.4858611089731411e-03,
           2.6933867548771758e-03, 1.7221140903704937e-03, 1.2407382229440791e-03,
           1.2435878276496955e-03, 2.4193952633248727e-03, 6.2339060289407083e-03,
           8.0309461712520876e-04, 9.3181986323789834e-04, 9.5783034332718700e-04,
           2.2660898636037261e-03, 1.7802796534180537e-03, 2.6571979312581875e-04,
           9.2634607111251918e-04, 1.5810185245264004e-03],
          [1.9493235258876179e-03, 1.9786027128591904e-03, 1.4140291972553610e-02,
           3.7201973506001745e-03, 4.3845466068066216e-04, 1.5304436972610567e-03,
           2.2097156829738759e-03, 2.8591871815612977e-03, 1.4301072616183181e-03,
           9.9437221166923172e-04, 1.3690958423974782e-03, 2.4402105140841090e-03,
           5.2943633069226512e-04, 7.5004227978192801e-04, 8.6016459857770028e-04,
           3.1466019144814608e-03, 2.2360795375444384e-03, 1.6159545671597605e-04,
           7.0048422794024819e-04, 1.2014015528772706e-03],
          [2.1674844853066858e-03, 1.5865480081162602e-03, 3.7201973506001745e-03,
           2.1274574617480089e-02, 3.9909227141697264e-04, 1.6481246723433428e-03,
           4.9158017471929655e-03, 2.5221102126636373e-03, 9.5384849402143984e-04,
           1.2347404942429857e-03, 1.5202051791453383e-03, 2.4453087721980561e-03,
           4.6429229320514104e-04, 7.6023722413111566e-04, 1.2373315413524663e-03,
           2.8035127901697272e-03, 1.8961512776990257e-03, 1.6218020183662784e-04,
           5.9842263937853702e-04, 1.3158365660538270e-03],
          [1.5903351423026848e-03, 3.9365984789376245e-04, 4.3845466068066216e-04,
           3.9909227141697264e-04, 1.1931352277704348e-02, 3.0937204045913537e-04,
           3.8338775043186374e-04, 7.6951976030099293e-04, 2.2976387481074697e-04,
           1.0956590131781735e-03, 1.5682982157153873e-03, 5.0124929379033781e-04,
           3.7717165634097634e-04, 5.1389991547056834e-04, 3.6111795849154795e-04,
           1.0432626586831986e-03, 9.3041313726939057e-04, 1.4474923964368156e-04,
           3.4603772624580643e-04, 1.3606607271146112e-03],
          [1.9242657898716525e-03, 2.4858611089731411e-03, 1.5304436972610567e-03,
           1.6481246723433428e-03, 3.0937204045913537e-04, 7.3292255467189687e-03,
           3.5385780499965817e-03, 1.3683038039160171e-03, 1.0489026828741754e-03,
           8.9102936026571569e-04, 1.6174411456311808e-03, 3.0968229715707327e-03,
           7.3993258722701268e-04, 5.4255147972143906e-04, 8.4668181752066874e-04,
           1.8931125300036275e-03, 1.3796838284921874e-03, 2.2737931366728891e-04,
           6.7584155312457842e-04, 1.1660966117775285e-03],
          [2.9879059292799641e-03, 2.6933867548771758e-03, 2.2097156829738759e-03,
           4.9158017471929655e-03, 3.8338775043186374e-04, 3.5385780499965817e-03,
           1.6133927472163669e-02, 1.9380952488713059e-03, 1.3667885452189439e-03,
           1.2192061706431622e-03, 2.0030316026648431e-03, 4.1322603720305197e-03,
           6.7909745467514783e-04, 8.5179405867513139e-04, 1.4216207127018586e-03,
           2.9539180653600089e-03, 2.0493063257644955e-03, 2.6488552587183780e-04,
           8.7044186256788659e-04, 1.6987763526262680e-03],
          [5.8158526388051033e-03, 1.7221140903704937e-03, 2.8591871815612977e-03,
           2.5221102126636373e-03, 7.6951976030099293e-04, 1.3683038039160171e-03,
           1.9380952488713059e-03, 3.7804346453413303e-02, 9.5813607255887238e-04,
           1.3849118546156933e-03, 2.0864716056392773e-03, 2.5392537741810947e-03,
           7.3281559749652399e-04, 1.1976708695723554e-03, 1.3641171883713547e-03,
           3.8342830901664762e-03, 2.1858459940987062e-03, 4.0740829083805248e-04,
           8.3467413018106177e-04, 1.8218235950233687e-03],
          [1.1076584657559144e-03, 1.2407382229440791e-03, 1.4301072616183181e-03,
           9.5384849402143984e-04, 2.2976387481074697e-04, 1.0489026828741754e-03,
           1.3667885452189439e-03, 9.5813607255887238e-04, 9.2802502369336622e-03,
           5.8089627083019206e-04, 9.8696608463236094e-04, 1.1873625842258938e-03,
           3.8264639620910225e-04, 8.1041076335565583e-04, 4.7770135861914477e-04,
           1.1052034635193162e-03, 7.4371746073077327e-04, 1.5168037757411286e-04,
           1.5213771111755425e-03, 6.4882907765797669e-04],
          [3.1880644746334580e-03, 1.2435878276496955e-03, 9.9437221166923172e-04,
           1.2347404942429857e-03, 1.0956590131781735e-03, 8.9102936026571569e-04,
           1.2192061706431622e-03, 1.3849118546156933e-03, 5.8089627083019206e-04,
           1.8441526588740136e-02, 1.1382470627796603e-02, 1.5655862274689192e-03,
           2.5081290988482057e-03, 3.0458868657559346e-03, 1.0068164685944146e-03,
           1.7225081689171561e-03, 2.6953622613315018e-03, 3.6183761166072852e-04,
           1.3821121844492116e-03, 1.1972663837662637e-02],
          [4.4186245468471547e-03, 2.4193952633248727e-03, 1.3690958423974782e-03,
           1.5202051791453383e-03, 1.5682982157153873e-03, 1.6174411456311808e-03,
           2.0030316026648431e-03, 2.0864716056392773e-03, 9.8696608463236094e-04,
           1.1382470627796603e-02, 3.7141460156350926e-02, 2.4634345023228079e-03,
           4.9293545515183088e-03, 5.4151301166464015e-03, 1.4146090399381900e-03,
           2.4277107072013821e-03, 3.3238031308707055e-03, 7.3206640617832933e-04,
           2.2096734692836624e-03, 9.4786263030457313e-03],
          [3.3466571942021082e-03, 6.2339060289407083e-03, 2.4402105140841090e-03,
           2.4453087721980561e-03, 5.0124929379033781e-04, 3.0968229715707327e-03,
           4.1322603720305197e-03, 2.5392537741810947e-03, 1.1873625842258938e-03,
           1.5655862274689192e-03, 2.4634345023228079e-03, 1.6113385590544604e-02,
           9.0876633395557617e-04, 9.4875149773685364e-04, 1.5773020912564391e-03,
           3.1016069999481111e-03, 2.3467014804084987e-03, 2.7198500003555514e-04,
           9.9908866586876396e-04, 1.9360424083099779e-03],
          [1.3412107617355408e-03, 8.0309461712520876e-04, 5.2943633069226512e-04,
           4.6429229320514104e-04, 3.7717165634097634e-04, 7.3993258722701268e-04,
           6.7909745467514783e-04, 7.3281559749652399e-04, 3.8264639620910225e-04,
           2.5081290988482057e-03, 4.9293545515183088e-03, 9.0876633395557617e-04,
           4.0477309321969848e-03, 1.1901770463553603e-03, 4.0824445213456919e-04,
           8.5603787638552766e-04, 1.0095451907679563e-03, 1.9872537223131380e-04,
           5.7145288352831449e-04, 2.3123361470140736e-03],
          [1.6360627863999076e-03, 9.3181986323789834e-04, 7.5004227978192801e-04,
           7.6023722413111566e-04, 5.1389991547056834e-04, 5.4255147972143906e-04,
           8.5179405867513139e-04, 1.1976708695723554e-03, 8.1041076335565583e-04,
           3.0458868657559346e-03, 5.4151301166464015e-03, 9.4875149773685364e-04,
           1.1901770463553603e-03, 1.8277684015431908e-02, 5.2528021756783813e-04,
           1.1939618185901600e-03, 1.1624184369750680e-03, 8.4917468952377874e-04,
           4.2392005745634370e-03, 2.5763052227920180e-03],
          [2.1568959784943114e-03, 9.5783034332718700e-04, 8.6016459857770028e-04,
           1.2373315413524663e-03, 3.6111795849154795e-04, 8.4668181752066874e-04,
           1.4216207127018586e-03, 1.3641171883713547e-03, 4.7770135861914477e-04,
           1.0068164685944146e-03, 1.4146090399381900e-03, 1.5773020912564391e-03,
           4.0824445213456919e-04, 5.2528021756783813e-04, 1.9066033679132538e-02,
           1.6662567934883051e-03, 1.3511005665728870e-03, 1.4152209821874487e-04,
           4.5224391125285910e-04, 1.2451325046931832e-03],
          [6.2524987419815400e-03, 2.2660898636037261e-03, 3.1466019144814608e-03,
           2.8035127901697272e-03, 1.0432626586831986e-03, 1.8931125300036275e-03,
           2.9539180653600089e-03, 3.8342830901664762e-03, 1.1052034635193162e-03,
           1.7225081689171561e-03, 2.4277107072013821e-03, 3.1016069999481111e-03,
           8.5603787638552766e-04, 1.1939618185901600e-03, 1.6662567934883051e-03,
           1.2585947097159817e-02, 4.7004857686835334e-03, 2.8731729176487776e-04,
           1.0299846310599138e-03, 2.3587292053265561e-03],
          [3.7180506975672363e-03, 1.7802796534180537e-03, 2.2360795375444384e-03,
           1.8961512776990257e-03, 9.3041313726939057e-04, 1.3796838284921874e-03,
           2.0493063257644955e-03, 2.1858459940987062e-03, 7.4371746073077327e-04,
           2.6953622613315018e-03, 3.3238031308707055e-03, 2.3467014804084987e-03,
           1.0095451907679563e-03, 1.1624184369750680e-03, 1.3511005665728870e-03,
           4.7004857686835334e-03, 1.2514818886617953e-02, 2.8575770858467209e-04,
           9.4161039895612720e-04, 3.6402328079338207e-03],
          [4.0281679108936688e-04, 2.6571979312581875e-04, 1.6159545671597605e-04,
           1.6218020183662784e-04, 1.4474923964368156e-04, 2.2737931366728891e-04,
           2.6488552587183780e-04, 4.0740829083805248e-04, 1.5168037757411286e-04,
           3.6183761166072852e-04, 7.3206640617832933e-04, 2.7198500003555514e-04,
           1.9872537223131380e-04, 8.4917468952377874e-04, 1.4152209821874487e-04,
           2.8731729176487776e-04, 2.8575770858467209e-04, 6.4699301717154852e-03,
           8.8744160259272527e-04, 3.5578318710317554e-04],
          [1.2999956675626666e-03, 9.2634607111251918e-04, 7.0048422794024819e-04,
           5.9842263937853702e-04, 3.4603772624580643e-04, 6.7584155312457842e-04,
           8.7044186256788659e-04, 8.3467413018106177e-04, 1.5213771111755425e-03,
           1.3821121844492116e-03, 2.2096734692836624e-03, 9.9908866586876396e-04,
           5.7145288352831449e-04, 4.2392005745634370e-03, 4.5224391125285910e-04,
           1.0299846310599138e-03, 9.4161039895612720e-04, 8.8744160259272527e-04,
           1.0246100213822419e-02, 1.5489827890922993e-03],
          [5.0679056444508912e-03, 1.5810185245264004e-03, 1.2014015528772706e-03,
           1.3158365660538270e-03, 1.3606607271146112e-03, 1.1660966117775285e-03,
           1.6987763526262680e-03, 1.8218235950233687e-03, 6.4882907765797669e-04,
           1.1972663837662637e-02, 9.4786263030457313e-03, 1.9360424083099779e-03,
           2.3123361470140736e-03, 2.5763052227920180e-03, 1.2451325046931832e-03,
           2.3587292053265561e-03, 3.6402328079338207e-03, 3.5578318710317554e-04,
           1.5489827890922993e-03, 1.9631915140537640e-02]]

        # Background frequencies for BLOSUM62.
        bg_freqs = [7.4216205067993410e-02, 5.1614486141284638e-02, 4.4645808512757915e-02,
         5.3626000838554413e-02, 2.4687457167944848e-02, 3.4259650591416023e-02,
         5.4311925684587502e-02, 7.4146941452644999e-02, 2.6212984805266227e-02,
         6.7917367618953756e-02, 9.8907868497150955e-02, 5.8155682303079680e-02,
         2.4990197579643110e-02, 4.7418459742284751e-02, 3.8538003320306206e-02,
         5.7229029476494421e-02, 5.0891364550287033e-02, 1.3029956129972148e-02,
         3.2281512313758580e-02, 7.2919098205619245e-02];

        # AA order.
        aas = np.array(list('ARNDCQEGHILKMFPSTWYV'))

        # Generate substituation probability matrix -- code from Chris Norn.
        M = np.zeros((20,20))
        S = np.zeros((20,20))
        for i in range(0, 20):
            for j in range(0, 20):
                M[i,j] = joint_p[i][j]
                S[i,j] = 2 * np.log(joint_p[i][j] / (bg_freqs[i] * bg_freqs[j])) / np.log(2)

        P_conditional = M / bg_freqs
        np.fill_diagonal(P_conditional, 0.0) # we don't want to make self substitutions.

        # Make BLOSUM62-informed substituations.
        for proto, seq in proto_object.current_sequences.items():

            for p in proto_object.mutable_positions[proto]:
                current_aa = seq[p]
                idx = np.argwhere(aas==current_aa)[0][0]
                sub_prob_renorm = P_conditional[:,idx] / P_conditional[:,idx].sum() # get subsitutation vector for that aa, and renormalise the vector.
                sub_prob = {a:f for a, f in list(zip(aas, sub_prob_renorm))}

                # Make mutations.
                seq = seq[:p] + np.random.choice(list(sub_prob.keys()), p=list(sub_prob.values())) + seq[p+1:]

            mutated_proto_sequences[proto] = seq

    elif method == 'pssm':
        # Mutate positions based on a user-defined PSSM.
        # PSSM needs to be a .csv of shape L,20
        # Each vector needs to sum to 1.
        # AA order.

        aas = np.array(list('ARNDCQEGHILKMFPSTWYV'))
        with open('pssm.csv', 'r') as f:
            lines = f.readlines().strip()
        pssm = np.array([l.split(',') for l in lines], dtype=float)

        for proto, seq in proto_object.current_sequences.items():

            for p in proto_object.mutable_positions[proto]:
                frequencies = pssm[p]
                seq = seq[:p] + np.random.choice(aas, p=frequencies) + seq[p+1:]

            mutated_proto_sequences[proto] = seq


    else:
        print(f'ERROR: sequence mutation method [{method}] does not exist. System exiting...')
        sys.exit()

    return mutated_proto_sequences
