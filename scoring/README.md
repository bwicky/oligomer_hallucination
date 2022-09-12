# Scoring scripts

Fold and score query sequences with AlphaFold2, followed by minimization and scoring with Rosetta.

`AF2.py` predicts the structure with AlphaFold2 and appends confidence metrics at the end of the generated PDBs.

`ROSETTA_min_score.xml` minimizes and scores oligomers with Rosetta.

*NB* you will need to have [AlphaFold2](https://github.com/deepmind/alphafold), [Rosetta](https://www.rosettacommons.org/software), [PyRosetta](https://www.pyrosetta.org/),  and [PSIPRED](http://bioinf.cs.ucl.ac.uk/psipred/)
 installed.

1. Update paths to AlphaFold2, Rosetta and PSIPRED
```
sed -i 's/\/projects\/ml\/alphafold\/alphafold_git\//<path_to_your_alphafol2_install>/g' AF2.py
sed -i 's/\/software\//<path_to_your_rosetta_install>/g' ROSETTA_min_and_score.xml generate_rosetta_cmds.sh
sed -i 's/\/home\/bwicky\/tools\//<path_to_your_psipred4_install>/g' ROSETTA_min_and_score.xml
```

2. Predict your sequences (combined in a `FASTA` file) with `AF2.py` and extract confidence metrics.
```
./AF2.py -f your_sequences.fasta
bash extract_af2_scores.sh *.pdb
```

3. Run Rosetta on these models
```
bash generate_rosetta_cmds.sh
bash rosetta_cmds
```

4. Compute SAP (Spatial Aggregation Propensity) on the Rosetta-minimized models
```
./PyROSETTA_sap.py rosetta_min_and_scored_models/*.pdb
```

