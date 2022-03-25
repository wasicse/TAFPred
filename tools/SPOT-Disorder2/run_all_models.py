#!/usr/bin/env python
'''
12/2017
run_all_models.py
version 1.0

Usage:
python run_all_models.py --gpu 0 --input_list /filepath/to/testlist

If there is no gpu on your machine, or tensorflow is not built to use GPU, then
change run_spotcontact.sh to run without a GPU.

If any packages are missing, please install them on your machine.

if there is no filelist specified, or it can't be opened, then this will run 
over all proteins in the input directory.
'''
import argparse
import os, glob, re, tqdm
import pandas as pd
import numpy as np
import spotdis2

#Get input parameters
parser = argparse.ArgumentParser()
parser.add_argument('--gpu', default=0, type=int, help='Which GPU to use on the machine')
parser.add_argument('--input_list',default='', type=str, help="Protein ids")
parser.add_argument('--input_dir',default='', type=str, help="Full filepath to the list of proteins (blank if in current directory)")
parser.add_argument('--batch_size',default=30, type=int, help="Batch size for processing")
args = parser.parse_args()
os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)

NUM_MODELS = len([i for i in  os.listdir('dat') if i[-5:] =='index'])

#Read protein list
try:
    with open(args.input_list,'r') as f:
        protein_ids = f.read().splitlines()
except:
    all_files = []
    necessary_files = ['.fasta','.pssm','.hhm','.spotcon','spot1d']
    for i in necessary_files:
        protein_fnames = glob.glob('inputs/?*'+i)
        all_files.append([re.search('/(\w+)'+i,a).group(1) for a in protein_fnames])
    protein_ids = list(set(all_files[0]).intersection(*all_files))

#Obtain all predictions
for i in tqdm.tqdm(range(NUM_MODELS)):
    spotdis2.model_pass(protein_ids,i,args)

#now average the files
for I,i in enumerate(protein_ids):
    model_protein_outputs = []
    for j in range(NUM_MODELS):
        with open('outputs/'+i+'.spotdis2.'+str(j),'r') as f:
            model_protein_outputs.append(pd.read_csv(f,delim_whitespace=True,comment='#',header=None).dropna().values[:,2])
    model_outputs = np.mean(model_protein_outputs,0).astype(float)
    with open(args.input_dir+'inputs/'+i+'.fasta','r') as f:
        seq = ''.join(f.readlines()[1:]).rstrip()
    with open('outputs/'+i+'.spotd2','w') as f:
        f.write('#SPOT-Disorder2 output for file %s\n'%(i))
        f.write('#\tAA\tP(D)\tO/D\n')
        for J,j in enumerate(seq):
            label='O' if model_outputs[J] < 0.370 else 'D'
            f.write('%i\t%s\t%1.3f\t%s\n'%(J,j,model_outputs[J],label))

