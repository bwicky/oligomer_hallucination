![header](imgs/oligo_hallucination.png)

# Oligomer hallucination with AlphaFold2

Design (hallucinate) cyclic-symmetric protein assemblies starting from only the specification of a homo-oligomeric valency.

Accompanying [oligomer hallucination paper](https://www.biorxiv.org/content/10.1101/2022.06.09.493773v1).

`HAL_design_models` contains models for each design that was experimentally tested.

`HAL_experimental` contains data about each design that was experimentally tested.

The [ProteinMPNN paper](https://www.biorxiv.org/content/10.1101/2022.06.03.494563v1) and [code](https://github.com/dauparas/ProteinMPNN).

## Getting started

1. Clone repo (NB you will also need to have [AlphaFold2](https://github.com/deepmind/alphafold) installed).
```
git clone https://github.com/bwicky/oligomer_hallucination
```

2. Create conda environment using the `SE3.yml` file.
```
cd oligomer_hallucination
conda env create -f SE3.yml
```

3. Change the shebang in `./oligomer_hallucination.py` to the location of your conda install.
```
#!/software/conda/envs/SE3/bin/python --> #!/location/of/your/conda/env/bin/python
```

4. Change paths to your AlphaFold2 install in `./modules/af2_net.py` and `./modules/losses.py`
```
sys.path.append('/projects/ml/alphafold/alphafold_git/') --> sys.path.append('/location/of/your/alphafold/installation/')
```

5. If using `tmalign` and `dssp` based losses, you will also need to install these packages ([TM-align](https://zhanggroup.org/TM-align/), [DSSP](https://anaconda.org/salilab/dssp)), and update the paths to their executables in `./modules/losses.py`


## Examples

- `./oligomer_hallucination.py --oligo AAA+ --L 100 --out example`

will perform design of a homo-trimer, with each protomer being composed of 100 amino acids.

- `./oligomer_hallucination.py --oligo AAAAAA+ --L 50 --loss dual_cyclic --out example` 

will perform design of a homo-hexamer, with each protomer being composed of 50 amino acids, and optimising for the `dual_cylic` loss.

- `./oligomer_hallucination.py --oligo AAAAAAAA+ --L 30 --single_chains` 

will perform design of a monomeric proteins containing eight repeats, each 30 amino acids in length.

- `./oligomer_hallucination.py --oligo AA+ --seq GDIQVQVNIDDNGKNFDYTYTVTTESELQKVLNELMDYIKKQGAKRVRISITARTKKEAEKFAAILIKVFAELGYNDINVTFDGDTVTVEGQLE`

will perform design of a homo-dimer, starting from the specified sequence.


## Outputs

- Structures (`.pdb`) for each accepted move of the MCMC trajectory.
- A log file (`.out`) containing the scores at each step of the MCMC trajectory (accepted and rejected).


## Features

Design (hallucination) is performed by MCMC search in sequence space, while optimizing (user-defined) losses composed of AlphaFold2 metrics, and/or geometric constraints, and/or secondary-structure definitions in order to match the design objective:

- Oligomers with arbitrary numbers of subunits can be designed.
- Multi-state design (either positive or negative) can be specified.
- MCMC trajectories can either be seeded with input sequence(s), or started randomly (using a frequency-adjusted AA distribution).
- At each step, position(s) are chosen for mutation (see `modules/mutations.py` for details about the different options).
- A configuration file (`.af2h` extension) can be used to specify designable positions and associated probabilities of mutation.
- Repeat proteins (sequence-symmetric monomers) can be designed instead of oligomers by passing the `--single_chains` flag.
- Specific amino acids can be exluded.
- MCMC paramters (initial temperature, annealing half-life, number of steps, tolerance for termination) can be specified.
- Currently implemented loss functions are (see `modules/losses.py` for details):
  - `plddt`: optimises pLDDT; seems to have trouble converging to assembled structures.
  - `ptm`: optmimises pTM; tends to 'melt' input structures.
  - `pae`: optimises pAE; similar to pTM.
  - `dual`: combination of `plddt` and `ptm` losses with equal weights.
  - `cyclic`: geometric loss term that enforce cyclic symmetry; minimises the standard deviation of the distances between center of masses of the protomers.
  - `dual_cyclic`: combination of dual and cyclic losses.
  - `aspect_ratio`: geometric term that enforces protomers with aspect ratios close to 1 (i.e. spherical and globular).
  - `pae_sub_mat`: initially implemented to enforce symmetry, under-tested.
  - `pae_asym`: this loss has different weights associated with the means of the different pAE sub-matrices (asymmetric weighting of the different inter-chain contacts). Off-diagonal elements (+/-1 from the diagaonl, and opposite corners) have higher weights.
  - `frac_dssp`: enforces an exact secondary structure content (`E, H, notEH`), computed by `DSSP` on the structure.
  - `min_frac_dssp`: enforces a minimum secondary structure content (`E, H, notEH`). The loss is minimal if the fraction of that secondary strucure (computed by `DSSP` on the structure) is larger than specified.
  - `tmalign`: loss defined as TM-score to a template PDB, given with `--template`. Alignment to the template (`tmalign -I`) can be forced with `--template_alignment [alignment].aln`. If the template has multiple chains, remove the `TER` in the PDB file.
  - `dual_tmalign`: jointly optimises ptm, plddt and tmalign losses.
  - `pae_asym_tmalign`: in development.
  - `entropy`: in development.


## Options
```
optional arguments:
  -h, --help            show this help message and exit
  --oligo OLIGO         oligomer(s) definitions (comma-separated string, no space). Numbers and types of subunits (protomers), and design type
                        (positive or negative) specifying each oligomer. Protomers are defined by unique letters, and strings indicate oligomeric
                        compositions. The last character of each oligomer has to be [+] or [-] to indicate positive or negative design of that
                        oligomer (e.g. AAAA+,AB+). The number of unique protomers must match --L / --seq. Must be specified.
  --L L                 lengths of each protomer (comma-separated, no space). Must be specified if not using --seq or a .af2h config file.
  --seq SEQ             seed sequence for each protomer (comma-separated, no space). Optional.
  --out OUT             the prefix appended to the output files. Must be specified.
  --single_chain        this option will generate sequence-symmetric repeat proteins instead of oligomers by removing chain breaks (default: False).
  --exclude_AA EXCLUDE_AA
                        amino acids to exclude during hallucination. Must be a continous string (no spaces) in one-letter code format (default: C).
  --mutation_rate MUTATION_RATE
                        number of mutations at each MCMC step (start-finish, stepped linear decay). Should probably be scaled with protomer length
                        (default: 3-1).
  --select_positions SELECT_POSITIONS
                        how to select positions for mutation at each step. Choose from [random, plddt::quantile, FILE.af2h::quantile]. FILE.af2h needs
                        to be a file specifying the probability of mutation at each site. Optional arguments can be given with :: e.g. plddt::0.25
                        will only mutate the 25% lowest plddt positions (default: random).
  --mutation_method MUTATION_METHOD
                        how to mutate selected positions. Choose from [uniform, frequency_adjusted, blosum62, pssm] (default: frequency_adjusted).
  --loss LOSS           the loss function used during optimization. Choose from [plddt, ptm, pae, dual, pae_sub_mat, pae_asym, entropy [not working
                        yet], dual_cyclic, tmalign (requires a --template), dual_tmalign (requires a --template), pae_asym_tmalign [not working yet],
                        aspect_ratio, frac_dssp or min_frac_dssp (requires a --dssp_fractions_specified), ]. Multiple losses can be combined as a
                        comma-separarted string of loss_name:args units (and weighed with --loss_weights).
                        loss_0_name::loss0_param0;loss0_param1,loss_1_name::[loss_1_configfile.conf] ... (default: dual).
  --loss_weights LOSS_WEIGHTS
                        if a combination of losses is passed, specify relative weights of each loss to the globabl loss by providing a comma-separated
                        list of relative weights. E.g. 2,1 will make the first loss count double relative to the second one (default: equal weights).
  --oligo_weights OLIGO_WEIGHTS
                        contribution of the loss of each oligomer to the global loss, provided as a comma-separted list of relative weights (default:
                        equal weights).
  --T_init T_INIT       starting temperature for simulated annealing. Temperature is decayed exponentially (default: 0.01).
  --half_life HALF_LIFE
                        half-life for the temperature decay during simulated annealing (default: 1000).
  --steps STEPS         number for steps for the MCMC trajectory (default: 5000).
  --tolerance TOLERANCE
                        the tolerance on the loss sliding window for terminating the MCMC trajectory early (default: None).
  --model MODEL         AF2 model (_ptm) used during prediction. Choose from [1, 2, 3, 4, 5] (default: 4).
  --amber_relax AMBER_RELAX
                        amber relax pdbs written to disk, 0=do not relax, 1=relax every prediction (default: 0).
  --recycles RECYCLES   the number of recycles through the network used during structure prediction. Larger numbers increase accuracy but linearly
                        affect runtime (default: 1).
  --msa_clusters MSA_CLUSTERS
                        the number of MSA clusters used during feature generation (?). Larger numbers increase accuracy but significantly affect
                        runtime (default: 1).
  --output_pae          output the PAE (predicted alignment error) matrix for each accepted step of the MCMC trajectory (default: False).
  --timestamp           timestamp output directly and every PDB written to disk with: %Y%m%d_%H%M%S_%f (default: False).
  --template TEMPLATE   template PDB for use with special loss functions (default: None).
  --dssp_fractions_specified DSSP_FRACTIONS_SPECIFIED
                        dssp fractions specfied for frac_dssp loss as E(beta sheet),H(alpha helix),notEH (other) e.g. 0.8,None,None will be enforce
                        80% beta sheet; or 0.5,0,None will enforce 50% beta sheet, no helices (default: None).
  --template_alignment TEMPLATE_ALIGNMENT
                        enforce tmalign alignment with fasta file (default: None).

```


## Example `.af2h` file

The following config file enables design at all positions set to 1 (equal probability of picking those sites for mutation), and disallow design at all positions that are set to 0.
```
>A
DEEQEKAEEWLKEAEEMLEQAKRAKDEEELLKLLVRLLELSVELAKIIQKTKDEEKKKELLEINKRLIEVIKELLRRLK
1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,0,0,1,1,0,1,1,0,0,1
>B
QEELAELIELILEVNEWLQRWEEEGLKDSEELVKEYEKIVEKIKELVKMAEEGHDEEEAEEEAKKLKKKAEEILREAEKG
1,1,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0
```

## Example template alignment for tmalign loss

Remove the `TER` in the template PDB file. `model1` (do not change the names) is the template given in `--template`, and `model2` should have the length of the protomer to be designed. The example below will design a 130 amino acid protein with motifs placed at the N- and C-termini (the sequence given here is arbitrary). Do not change this order!

```
>model1
RSMSWDNEVAFN-----------------------------------------------------
----------------------------------------------------QHHLGGAKQAGAV

>model2
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```


## Citing this work

If you use the code, please cite:

```bibtex
@article {HAL2022,
	author = {Wicky, B. I. M. and Milles, L. F. and Courbet, A. and Ragotte, R. J. and Dauparas, J. and Kinfu, E. and Tipps, S. and Kibler, R. D. and Baek, M. and DiMaio, F. and Li, X. and Carter, L. and Kang, A. and Nguyen, H. and Bera, A. K. and Baker, D.},
	title = {Hallucinating protein assemblies},
	year = {2022},
	doi = {10.1101/2022.06.09.493773},
	URL = {https://www.biorxiv.org/content/early/2022/06/09/2022.06.09.493773},
	eprint = {https://www.biorxiv.org/content/early/2022/06/09/2022.06.09.493773.full.pdf},
	journal = {bioRxiv}
```

## Acknowledgements

This work was made possible by the following separate libraries and packages:

*   [AlphaFold2](https://github.com/deepmind/alphafold)
*	[RosetTTAfold](https://github.com/RosettaCommons/RoseTTAFold)
*   [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
*	[Rosetta](https://www.rosettacommons.org/software)
*   [Biopython](https://biopython.org)
*   [Matplotlib](https://matplotlib.org/)
*   [Seaborn](https://seaborn.pydata.org/)
*   [NumPy](https://numpy.org)
*   [PyCORN](https://github.com/pyahmed/PyCORN)
*   [Pandas](https://pandas.pydata.org/)
*   [SciPy](https://scipy.org)
*   [Scikit-learn](https://scikit-learn.org/stable/)
*   [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
*	[PSIPRED](http://bioinf.cs.ucl.ac.uk/psipred/)
*   [TM-align](https://zhanggroup.org/TM-align/)
*   [MM-align](https://zhanggroup.org/MM-align/)

We thank all their contributors and maintainers!

## Get in touch

Questions and comments are welcome:

* Basile Wicky: [bwicky@uw.edu](mailto:bwicky@uw.edu)
* Lukas Milles: [lmilles@uw.edu](mailto:lmilles@uw.edu)
* Alexis Courbet [acourbet@uw.edu](mailto:acourbet@uw.edu)
