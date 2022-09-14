#!/global/cfs/cdirs/m3962/conda/pyrosetta 
#!/global/cfs/cdirs/m3962/conda/SE3/bin/python 

import os, sys
import pyrosetta
import numpy as np
import glob
import shutil
import argparse
import string
import random
from sklearn.cluster import KMeans

####################################################################################################################

#generation_number = sys.argv[1]
#selection_pressure = sys.argv[1]

parser = argparse.ArgumentParser(description='This script sets up jobs for in silico evolution of hallucinated designs. Selection for new generation based on ptm and plddt, takes only the fittest designs that are 1sigma above the rest of the population.')

parser.add_argument('--generation', 
                    default=None,
                    required=True,
                    help='The evolution step number for the next generation')

parser.add_argument('--selection_pressure',
                    type=float,
                    required=False, 
                    help='Selection pressure, only the fittest survive... (in number of sigmas above the population mean) (default: %(default)s).'  )                                                                                                                                                                                            

parser.add_argument('--number_of_families',
                    type=float,
                    required=False,
                    help='How many of the fittest families will survive...? ')

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
generation_number = args.generation
selection_pressure = args.selection_pressure
number_of_families = args.number_of_families

####################################################################################################################

#Compute mean + sigma of plddt, ptm and fitness over design population
def find_the_fittest(pdbs, sel_pressure):
    dict_of_scores ={'plddt':[],'ptm':[], 'fitness':[]}
    for pdb in pdbs:
        with open(pdb,'r') as pdb_to_evaluate:
            lines = pdb_to_evaluate.readlines()
            for ln in lines:
                if ln.startswith("plddt "):
                    dict_of_scores['plddt'].append(float(ln[5:11]))
                if ln.startswith("ptm "):
                    dict_of_scores['ptm'].append(float(ln[5:9]))
                if ln.startswith("loss "):
                    dict_of_scores['fitness'].append(1 - float(ln[5:9]))
    avg_plus_std_Dict = {}
    for k,v in dict_of_scores.items(): 
        avg_plus_std_Dict[k] = np.mean(v) + sel_pressure*np.std(v)
    print(" All designs fitter than "+ str(avg_plus_std_Dict) + " will survive...")
    return avg_plus_std_Dict

#Clusters design population by fitness level, and select fitest individuals per clusters
def find_the_fittest_clusters(pdbs, n_clusters):
    data = []
    for pdb in pdbs:
        with open(pdb,'r') as pdb_to_evaluate:
            lines = pdb_to_evaluate.readlines()
            for ln in lines:
                if ln.startswith("loss "):
                    data.append(round(1 - float(ln[5:11]),3))
    kmeans = KMeans(int(n_clusters))
    kmeans.fit(np.array(data).reshape(-1,1))
    cluster_partitions = {}
    for centr in kmeans.cluster_centers_:
        centroid_label = kmeans.predict([centr])
        partition = []
        for k, v in zip(data, kmeans.labels_):
            if v == centroid_label:
                partition.append(k)
        cluster_partitions[centroid_label[0]] = partition

    max_cluster_partitions = {}
    for k,v in cluster_partitions.items():
        max_cluster_partitions[k] = np.max(v)
    print(" Clustered fittest families as " + str(max_cluster_partitions))
    return max_cluster_partitions

#ID generator
def id_generator(size=3, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

#Loop through pdbs and yield symmetry number and sequence of subunit of designs matching the selection pressure
def evaluate_pdbs(pdb, fitness_dict):
    with open(pdb,'r') as pdb_to_evaluate:
        lines = pdb_to_evaluate.readlines()
        scores = {}
        for ln in lines:
            if ln.startswith("plddt "):
                scores['plddt'] = (float(ln[5:11]))
            if ln.startswith("ptm "):
                scores['ptm'] = (float(ln[5:9]))
            if ln.startswith("loss "):
                scores['fitness'] = round(1 - float(ln[5:9]),3)
        if scores.get("fitness") is not None:
            if scores.get("fitness") >= fitness_dict.get("fitness"):
                fit_value = scores.get("fitness")
                print('Selected trajectory for evolution: \n' + pdb + ' \nWith fitness = '  + str(fit_value) )
                number_of_As = test = pdb.split("+",1)[0].rsplit('_')[-1]
                symmetry = number_of_As.count('A')
                pose = pyrosetta.pose_from_pdb(pdb) 
                sequence = pose.sequence()
                number_of_residues = pose.total_residue()
                number_of_residue_per_subunit = number_of_residues / symmetry
                sequence_of_subunit = sequence[:int(number_of_residue_per_subunit)] 
                print("number of subunits: " + number_of_As) 
                print("sequence: " + sequence_of_subunit)
                print("Symmetry: " + str(symmetry))
                return number_of_As, sequence_of_subunit, symmetry, fit_value
            else:
                print("Passing trajectory " + pdb + " with scores " + str(scores) )
        else:
            print("Passing trajectory " + pdb )

#Loop through pdbs and yield symmetry number and sequence of subunit of designs matching the selection pressure
def evaluate_clusters(pdb, fitness_dict):
    with open(pdb,'r') as pdb_to_evaluate:
        lines = pdb_to_evaluate.readlines()
        scores = {}
        for ln in lines:
            if ln.startswith("loss "):
                scores['fitness'] = round(1 - float(ln[5:11]),3)
        if scores.get("fitness") is not None:
            if scores.get("fitness") in fitness_dict.values():
                fit_value = scores.get("fitness")
                print('Selected trajectory for evolution: \n' + pdb + ' \nWith fitness = '  + str(fit_value) )
                number_of_As = test = pdb.split("+",1)[0].rsplit('_')[-1]
                symmetry = number_of_As.count('A')
                pose = pyrosetta.pose_from_pdb(pdb) 
                sequence = pose.sequence()
                number_of_residues = pose.total_residue()
                number_of_residue_per_subunit = number_of_residues / symmetry
                sequence_of_subunit = sequence[:int(number_of_residue_per_subunit)] 
                print("number of subunits: " + number_of_As) 
                print("sequence: " + sequence_of_subunit)
                print("Symmetry: " + str(symmetry))
                return number_of_As, sequence_of_subunit, symmetry, fit_value
            else:
                print("Passing trajectory " + pdb + " with scores " + str(scores) )
        else:
            print("Passing trajectory " + pdb )

#Write each command into task file into new evolution dir
def write_to_task_file(task_file, number_ofAs, sequence_of_sub, symm, va, gen_numb, ID):
    with open(task_file, 'a') as the_file:
        the_file.write('srun /global/cfs/cdirs/m3962/users/acourbet/AF2_multistate_hallucination.py --oligo ' 
                        + str(number_ofAs) + '+' 
                        + ' --loss dual_cyclic --seq ' + str(sequence_of_sub)
                        + ' --single_chain --mutation_rate 3-1 --T_init 0.01 --half_life 800 --steps 5000 --out HALC' + str(symm) 
                        + '_' + ID + '_generation' + str(gen_numb) + '_fitness' + str(va)
                        + ' --timestamp --recycles 1 --model 4 --select_positions random ; \n' )

#Write SLURM submission script and move everything into new dir
def make_new_dir(gen_number, task_file):
    jobs_file_name = 'HAL_jobs_gen' + str(gen_number) + '.sh'
    if not os.path.exists(jobs_file_name):
        with open(jobs_file_name, 'a') as the_otherfile:
            the_otherfile.write( 
                            '#!/bin/bash \n'
                            '#SBATCH --job-name=HAL_gen' + str(gen_number) + '\n'
                            '#SBATCH --constraint=gpu \n'
                            '#SBATCH --gpus=4   \n'
                            '#SBATCH --nodes=1  \n'
                            '#SBATCH --ntasks=4 #4 gpus per node and 1 task per gpu \n'
                            '#SBATCH --account=m3962_g  \n'
                            '#SBATCH --qos=regular \n'
                            '#SBATCH --time=2:00:00 \n'
                            '#SBATCH --gpu-bind=map_gpu:0,1,2,3 \n'
                            'module load cudatoolkit/21.9_11.4 \n'
                            'source ~/.bashrc \n'
                            'source activate SE3 \n'
                            'conda activate mlfold \n'
                            'CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ' + str(task_file) + ') \n'
                            'echo "${CMD}" | bash \n'
                            '#####sbatch -a 1-$(cat ' + str(task_file) + '|wc -l) *.sh'
                            )

    dirName = 'generation' + str(gen_number)

    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory " , dirName ,  "for new evolution step created ")
    else:    
        print("Directory " , dirName ,  " already exists")

    shutil.move(task_file, dirName)
    shutil.move(jobs_file_name, dirName)

####################################################################################################################

##################
pyrosetta.init("-mute all")
##################

def main():                      
    pdb_files = glob.glob( "**/*.pdb", recursive = True)
    tasks = open('tasks_gen' + str(generation_number), "x")
    name_of_task_file = tasks.name

    if args.number_of_families is None :
        fittest_designs = find_the_fittest(pdb_files, selection_pressure)
        for pdb in pdb_files:
            if evaluate_pdbs(pdb, fittest_designs) is not None:
                As, seq, sym, val = evaluate_pdbs(pdb, fittest_designs)
                unique_ID = id_generator()
                write_to_task_file(name_of_task_file, As, seq, sym, val, generation_number, unique_ID)
            else:
                continue
        make_new_dir(generation_number, name_of_task_file)

    elif args.selection_pressure is None :
        fittest_clusters = find_the_fittest_clusters(pdb_files, number_of_families)
        for pdb in pdb_files:
            if evaluate_clusters(pdb, fittest_clusters) is not None:
                As, seq, sym, val = evaluate_clusters(pdb, fittest_clusters)
                unique_ID = id_generator()
                write_to_task_file(name_of_task_file, As, seq, sym, val, generation_number, unique_ID)
            else:
                continue
        make_new_dir(generation_number, name_of_task_file)

main()

