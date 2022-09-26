# Experimental data

`HAL_exp_data.csv` contains all the experimental data collected on designs reported in the [oligomer hallucination paper](https://doi.org/10.1126/science.add1964) (*N* = 352).

The complete dataframe that also contains the raw SEC and CD data is available on Zenodo as a `.h5` file [![DOI](https://zenodo.org/badge/533942520.svg)](https://zenodo.org/badge/latestdoi/533942520).

*NB* for some designs not all entries are populated because they were dropped early in the characterization process due to unfavourable behaviours (e.g. low soluble expression). 

The dataframe contains the following features, grouped by categories:
- Design specifications: 
```
protomer_seq, protomer_length, seq_origin, symmetry,
pseudo_symmetry, oligomer_length, repeats_per_protomer,
repeat_length, repeat_seq, repeat_occurs, ID, large_HAL,
xtal, nsEM, cryoEM, high_res_struct, 
```
- Cloning results: 
```
eblock, into_plasmid, plasmid_seq, coding_seq, exp_prod,
exp_prod_length, protomer_mw, e280, OD280, charge@7.4, pI,
```
- Wetlab features:
```
culture_vol, soluble, sec_column, main_peak, integral,
nPeaks, monodisperse, CD_conc_mgmL, Kav, obs_mw_sec,
SECMALS_conc_mgmL, obs_mw_secmals, correct_secmals,
error_Da_secmals, correct_Vel, Kav_postmelt,
correct_Vel_postmelt, same_Vel_postmelt, main_peak_postmelt,
integral_postmelt, nPeaks_postmelt, monodisperse_postmelt,
tot_sol_yield, tot_sol_yield_per_Leq, oligo_mw,
log10_oligo_mw
```
- AF2 models and confidence metrics:
```
pdb_name, best_ptm_model, ptm_model_1, ptm_model_2,
ptm_model_3, ptm_model_4, ptm_model_5, pae_range_model_1,
pae_range_model_2, pae_range_model_3, pae_range_model_4,
pae_range_model_5, pae_min_model_1, pae_min_model_2,
pae_min_model_3, pae_min_model_4, pae_min_model_5,
pae_off_diag_range_model_1, pae_off_diag_range_model_2,
pae_off_diag_range_model_3, pae_off_diag_range_model_4,
pae_off_diag_range_model_5, pae_off_diag_mean_model_1,
pae_off_diag_mean_model_2, pae_off_diag_mean_model_3,
pae_off_diag_mean_model_4, pae_off_diag_mean_model_5,
plddt_min_model_1, plddt_min_model_2, plddt_min_model_3,
plddt_min_model_4, plddt_min_model_5,
pae_off_diag_min_model_1, pae_off_diag_min_model_2,
pae_off_diag_min_model_3, pae_off_diag_min_model_4,
pae_off_diag_min_model_5, pae_median_model_1,
pae_median_model_2, pae_median_model_3, pae_median_model_4,
pae_median_model_5, plddt_median_model_1,
plddt_median_model_2, plddt_median_model_3,
plddt_median_model_4, plddt_median_model_5,
pae_off_diag_max_model_1, pae_off_diag_max_model_2,
pae_off_diag_max_model_3, pae_off_diag_max_model_4,
pae_off_diag_max_model_5, plddt_range_model_1,
plddt_range_model_2, plddt_range_model_3,
plddt_range_model_4, plddt_range_model_5,
pae_off_diag_std_model_1, pae_off_diag_std_model_2,
pae_off_diag_std_model_3, pae_off_diag_std_model_4,
pae_off_diag_std_model_5, pae_mean_model_1, pae_mean_model_2,
pae_mean_model_3, pae_mean_model_4, pae_mean_model_5,
pae_std_model_1, pae_std_model_2, pae_std_model_3,
pae_std_model_4, pae_std_model_5, pae_max_model_1,
pae_max_model_2, pae_max_model_3, pae_max_model_4,
pae_max_model_5, plddt_mean_model_1, plddt_mean_model_2,
plddt_mean_model_3, plddt_mean_model_4, plddt_mean_model_5,
pae_off_diag_median_model_1, pae_off_diag_median_model_2,
pae_off_diag_median_model_3, pae_off_diag_median_model_4,
pae_off_diag_median_model_5, plddt_max_model_1,
plddt_max_model_2, plddt_max_model_3, plddt_max_model_4,
plddt_max_model_5, plddt_std_model_1, plddt_std_model_2,
plddt_std_model_3, plddt_std_model_4, plddt_std_model_5,
```
- Rosetta scores:
```
9mer, helix_sc, holes, loop_sc, nomega_off,
packstat, psipred_match, psipred_mismatch_prob,
psipred_probability, ss_sc, frac_DSSP_extended,
frac_DSSP_helix, pyrosetta_sap, frac_R, frac_H, frac_K,
frac_D, frac_E, frac_S, frac_T, frac_N, frac_Q,
frac_C, frac_G, frac_P, frac_A, frac_V, frac_I,
frac_L, frac_M, frac_F, frac_Y, frac_W,
total_score_per_res, buried_npsa_per_res, cav_vol_per_res,
dslf_fa13_per_res, exposed_hydrophobics_per_res,
exposed_npsa_per_res, fa_atr_per_res, fa_dun_dev_per_res,
fa_dun_rot_per_res, fa_dun_semi_per_res, fa_elec_per_res,
fa_intra_atr_xover4_per_res, fa_intra_elec_per_res,
fa_intra_rep_xover4_per_res, fa_intra_sol_xover4_per_res,
fa_rep_per_res, fa_sol_per_res, hbond_bb_sc_per_res,
hbond_lr_bb_per_res, hbond_sc_per_res, hbond_sr_bb_per_res,
hxl_tors_per_res, hydrophobic_sasa_per_res, lk_ball_per_res,
lk_ball_bridge_per_res, lk_ball_iso_per_res,
nomega_off_per_res, omega_per_res, p_aa_pp_per_res,
polar_sasa_per_res, pro_close_per_res, rama_prepro_per_res,
ref_per_res, sasa_per_res, total_charge_per_res,
total_sasa_per_res, pyrosetta_sap_per_res,
```
- BLAST results:
```
sseqid_blast_repeat, qlen_blast_repeat, slen_blast_repeat,
length_blast_repeat, pident_blast_repeat,
evalue_blast_repeat, bitscore_blast_repeat,
score_blast_repeat, qseq_blast_repeat, sseq_blast_repeat,
sseqid_blast_protomer, qlen_blast_protomer,
slen_blast_protomer, length_blast_protomer,
pident_blast_protomer, evalue_blast_protomer,
bitscore_blast_protomer, score_blast_protomer,
qseq_blast_protomer, sseq_blast_protomer, 
```
- TM-align against the PDB results:
```
tmalign_to, tmscore_design_norm, tmscore_pdb_norm,
length_pdb, tmaligned_length, tmalign_rmsd, tmalign_seqid,
```
- MM-align against the biounits from the PDB results:
```
mmalign_to, mmscore_design_norm, mmscore_biounit_norm, 
length_biounit, mmaligned_length, mmalign_rmsd, mmalign_seqid, 
```

- Only available in the `.h5` file (arrays):
```
vol, A280, vol_norm, A280_norm, peaks, peak_heights,
temp, mdeg@222nm, WL,
mdeg@25C, mdeg@35C, mdeg@45C, mdeg@55C, mdeg@65C,
mdeg@75C, mdeg@85C, mdeg@95C, mdeg@25C_postmelt,
vol_postmelt, A280_postmelt,
MRE@222nm, MRE@25C, MRE@35C, MRE@45C,
MRE@55C, MRE@65C, MRE@75C, MRE@85C, MRE@95C,
MRE@25C_postmelt, vol_norm_postmelt, A280_norm_postmelt,
peaks_postmelt, peak_heights_postmelt,
```
