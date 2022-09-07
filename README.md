Oligomer and Multi-State design using AlphaFold2 MCMC hallucination.
====================================================================================
<bwicky@uw.edu> 
<lmilles@uw.edu> 

Summary
-------
- Designs (hallucinations) are performed by MCMC searches in sequence space and optimizing (user-defined) losses composed of AlphaFold2 metrics, and/or geometric constraints, and/or secondary-structure definitions.
- Oligomers with arbitrary number of subunits can be designed.
- Multistate design (either positive or negative) can be specified.
- MCMC trajectories can either be seeded with input sequence(s), or started randomly (using a frequency-adjusted AA distribution).
- At each step, position(s) are chosen for mutation based on different options (see `modules/seq_mutation.py` for details).
- A 'resfile' (`.af2h` extension) can be employed to specify designable positions and associated probabilities of mutation.
- The oligomeric state (number of subunits) for each oligomer (state) can be specified.
- Repeat proteins (sequence-symmetric monomers) can be designed instead of oligomers by passing the `--single_chains` flag. 
- Specific amino acids can be exluded.
- MCMC paramters (initial temperature, annealing half-life, steps, tolerance) can be specified.
- Currently implemented loss functions are (see `modules/losses.py` for details):
  - `plddt`: plDDT seem to have trouble converging to complex formation.
  - `ptm`: pTM tends to 'melt' input structures.
  - `pae`: similar to result as ptm?
  - `dual`: combination of plddt and ptm losses with equal weights.
  - `entropy`: current implementation unlikely to work.
  - `pae_sub_mat`: initially made to enforce symmetry, but probably not working.
  - `pae_asym`: this loss has different weights associated with the means of the different PAE sub-matrices (asymmetric weighting of the different inter-chain contacts). Off-diagonal elements (+/-1 from the diagaonl, and opposite corners) have higher weights.
  - `cyclic`: geometric loss term to enforce cyclic symmetry, minimizes the standard deviation of protomer center of mass distances
  - `dual_cyclic`: dual with an added geometric loss term to enforce symmetry. 
  - `frac_dssp`: Enforcing an exact secondary structure percentage (E,H,notEH) as computed by DSSP on the structure.
  - `min_frac_dssp`: Enforcing a minimum fraction secondary structure content (E,H,notEH), loss is best if fraction of that sec.strucure is larger than specified as computed by DSSP on the structure.
  - `tmalign`: loss defined as TM-score to template PDB, given with `--template` , alignment of template (`tmalign -I`) can be forced with `--template_alignment [alignment].aln` if template has multiple chains remove the `TER` in the pdbfiles.
  - `dual_tmalign`: jointly optimises ptm, plddt and tmalign (see above) TM-score.
  - `pae_asym_tmalign`: in development.
  - `aspect_ratio`: geometric term that enforces protomers with aspect ratios close to 1 (i.e. spherical).  

Minimal inputs
--------------
- The number and type of subunits for each oligomer, also indicating whether it is a positive or negative design task.
- The length of each protomer or one seed sequence per protomer.

Examples
--------
- `./AF2_multistate_hallucination.py --oligo AAAAAA+ --L 30 --loss dual_cyclic --out C6` will perform design of an oligomeric assembly six protomers in C6 symmetry, each 30 amino-acids in length.
- `./AF2_multistate_hallucination.py --oligo AAAAAA+ --L 30 --single_chains` will perform single-state design of a monomeric repeat proteins containing six repeats, each 30 amino-acids in length.

Example `.af2h` file
--------------------
The following config file enables design at all positions set to 1 (equal probability of picking those sites for mutation), and disallow design at all positions that are set to 0.
```
>A
DEEQEKAEEWLKEAEEMLEQAKRAKDEEELLKLLVRLLELSVELAKIIQKTKDEEKKKELLEINKRLIEVIKELLRRLK
1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,0,0,1,1,0,1,1,0,0,1
>B
QEELAELIELILEVNEWLQRWEEEGLKDSEELVKEYEKIVEKIKELVKMAEEGHDEEEAEEEAKKLKKKAEEILREAEKG
1,1,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0
```

Example template alignment for tmalign loss
-------------------------------------------
Remove the `TER` in the template pdbfiles.
model1 (do not change the names) ist the template given in `--template` 
model2 should be the length of the protomer to be designed (sequence given here is irrelevant)
e.g. for a desing of length 130 with motifs placed at N- and C-termini
Do not change this order!
```
>model1
RSMSWDNEVAFN-----------------------------------------------------
----------------------------------------------------QHHLGGAKQAGAV

>model2
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

Outputs
-------
- PDB structures for each accepted move of the MCMC trajectory.
- A file (.out) containing the scores at each step of the MCMC trajectory (accepted and rejected).

