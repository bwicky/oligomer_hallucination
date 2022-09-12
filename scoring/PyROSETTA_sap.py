#!/home/lmilles/.conda/envs/jupyter3/bin/python3
import sys
import os
import pyrosetta
import numpy as np
import pandas as pd

df = pd.DataFrame(columns=["pdb", "total_sap_score",  "sap_per_res", "no_res"])

pyrosetta.init("-ignore_unrecognized_res")

pdbfiles = sys.argv[1:]
print(pdbfiles)

with open("sap_scores.csv", "w") as wf:
    wf.write("#file,total_sap_score,sap_per_res,no_res\n")
    for pdb_path in sorted(pdbfiles):

        pose = pyrosetta.pose_from_file(pdb_path)

        true_sel = (
            pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
        )

        total_sap_score = pyrosetta.rosetta.core.pack.guidance_scoreterms.sap.calculate_sap(
            pose,
            true_sel ,
            true_sel ,
            true_sel ,
        )
        per_res_sap = pyrosetta.rosetta.core.pack.guidance_scoreterms.sap.calculate_per_res_sap(
            pose,
            true_sel ,
            true_sel ,
            true_sel ,
        )

        no_res = len(pose.sequence())
        sap_per_res = total_sap_score / float(no_res) 

        print(pdb_path, f":\n\t sap_score: {total_sap_score} | sap_per_res: {sap_per_res} | no_res = {no_res} ")
        wf.write(f"{pdb_path},{total_sap_score},{sap_per_res},{no_res}\n")
        
        df.loc[len(df.index)] = [ pdb_path, total_sap_score, sap_per_res, no_res ]


df = df.sort_values(by=['total_sap_score'], ascending=True)
df.to_csv("sap_scores_sorted.csv")
