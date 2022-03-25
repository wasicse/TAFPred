import numpy as np
import pandas as pd
import tqdm, os, glob, re
import cPickle as pickle
import tensorflow as tf
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
#------------------------FEATURE EXTRACTION FUNCTIONS--------------------------#
def spot1d_feature_sincos(x,seq,norm_ASA=True):
    ASA = x[:,0]
    if norm_ASA==True:
        rnam1_std = "ACDEFGHIKLMNPQRSTVWY-X"
        ASA_std = (115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
                            185, 160, 145, 180, 225, 115, 140, 155, 255, 230,1,1)
        dict_rnam1_ASA = dict(zip(rnam1_std, ASA_std))
        ASA_div =  np.array([dict_rnam1_ASA[i] for i in seq])
        ASA = (ASA/ASA_div)[:,None]
    else:
        ASA = ASA[:,None]
    HSE = x[:,1:4]
    angles = x[:,4:8]
    SSprob3_8 = x[:,8:] if np.max(x[:,8:])<1.01 else x[:,8:]/100. #in case probs are in 0-1 or 0-100
    angles = np.deg2rad(angles)
    angles = np.concatenate([np.sin(angles),np.cos(angles)],1)
    return np.concatenate([ASA,HSE,angles,SSprob3_8],1)

def read_pssm(fname,seq):
    try:
        num_pssm_cols = 44
        pssm_col_names = [str(j) for j in range(num_pssm_cols)]
        with open(fname,'r') as f:
            tmp_pssm = pd.read_csv(f,delim_whitespace=True,names=pssm_col_names).dropna().values[:,2:22].astype(float)
        if tmp_pssm.shape[0] != len(seq):
            raise ValueError('PSSM file is in wrong format or incorrect!')
    except: #Rare case that BLOSUM62 is used instead
        num_pssm_cols = 22
        pssm_col_names = [str(j) for j in range(num_pssm_cols)]
        with open(fname,'r') as f:
            tmp_pssm = pd.read_csv(f,delim_whitespace=True,names=pssm_col_names).dropna().values[:,2:].astype(float)
        if tmp_pssm.shape[0] != len(seq):
            raise ValueError('PSSM file is in wrong format or incorrect!')    
    return tmp_pssm

def read_hhm(fname,seq):
    num_hhm_cols = 22
    hhm_col_names = [str(j) for j in range(num_hhm_cols)]
    with open(fname,'r') as f:
        hhm = pd.read_csv(f,delim_whitespace=True,names=hhm_col_names)
    pos1 = (hhm['0']=='HMM').idxmax()+3
    num_cols = len(hhm.columns)
    hhm = hhm[pos1:-1].values[:,:num_hhm_cols].reshape([-1,44])
    hhm[hhm=='*']='9999'
    if hhm.shape[0] != len(seq):
        raise ValueError('HHM file is in wrong format or incorrect!')
    return hhm[:,2:-12].astype(float)

#USE THIS TO READ SPOT1D FILES
def read_spot1d_output(fname,seq):
    with open(fname,'r') as f:
        spot1d_features = pd.read_csv(f,delim_whitespace=True).values[:,4:].astype(float)
    tmp_spot1d = spot1d_feature_sincos(spot1d_features,seq)
    if tmp_spot1d.shape[0] != len(seq):
        raise ValueError('Spider3 file is in wrong format or incorrect!')
    return tmp_spot1d


def obtain_feats(id,mu,std,input_dir=''):
    features =      [['.spot1d',1,read_spot1d_output],
                     ['.hhm',1,read_hhm],
                     ['.pssm',1,read_pssm]]
    with open(input_dir+'inputs/'+id+'.fasta','r') as f:
        seq = ''.join(f.readlines()[1:]).rstrip()
    inputs_twoD = []
    inputs_oneD = []
    for i in features:
        tmp_feats = i[2](input_dir+'inputs/'+id+i[0],seq)
        if i[1] == 1:
            inputs_oneD.append(tmp_feats)
        if i[1] == 2:
            inputs_twoD.append(tmp_feats)
    inputs_oneD = np.concatenate(inputs_oneD,1)
    inputs_oneD = (inputs_oneD-mu)/std
    return inputs_oneD

#-------------------------------MODEL FUNCTIONS---------------------------------#
def softmax(x,axis=1):
    return np.exp(x)/(np.sum(np.exp(x),axis)[:,None])

def sigmoid(x):
    return 1/(1+np.exp(-x))
def model_pass(ids,model_id,args):
    os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)
    #---------------------------------MODEL------------------------------------#
    config = tf.ConfigProto()
    #config.gpu_options.allow_growth = True
    config.allow_soft_placement=True
    config.log_device_placement=False
    saver = tf.train.import_meta_graph('dat/model_'+str(model_id)+'.meta',clear_devices=True)
    feat_depth = tf.get_default_graph().get_tensor_by_name('oneD_feats:0').get_shape().as_list()[-1]
    with tf.Session(config=config) as sess:
        disout = []
        all_seq_lens = [0]
        all_seqs = []

        saver.restore(sess,'dat/model_'+str(model_id))
        model = tf.get_collection("output")
        placeholders = ['oneD_feats:0','twoD_feats:0','seq_lens:0','ph_dropout:0','train_bool:0','mask_bool:0','ln_mask_bool:0']
        batch_size = args.batch_size
        #---------------------------NORMALISATION DATA-----------------------------#
        with open('dat/norm_params.p','r') as f:
            normdic = pickle.load(f)
            normmu = normdic['mu1D']
            normstd = normdic['std1D']
        for i in range(0,len(ids),batch_size):   
            seqs = []    
            inputs = []
            for j in range(i,min(i+batch_size,len(ids))):
                with open(args.input_dir+'inputs/'+ids[j]+'.fasta','r') as f:
                    seqs.append(''.join(f.readlines()[1:]).rstrip())
                inputs.append(obtain_feats(ids[j],normmu,normstd,args.input_dir))
            all_seqs = all_seqs + seqs
            seq_lens = [j.shape[0] for j in inputs]
            maxseqlen = np.max(seq_lens)
            inputs = [np.concatenate([j,np.zeros([maxseqlen-seq_lens[J],j.shape[1]])],0) for J,j in enumerate(inputs)]
            mask = np.array([np.concatenate([np.ones(j),np.zeros(maxseqlen-j)]) for j in seq_lens])
            feed_dict = {placeholders[0]:inputs,placeholders[1]:np.zeros([0,0,0,1]),placeholders[2]:seq_lens,placeholders[3]:1,placeholders[4]:False,placeholders[5]:mask,placeholders[6]:mask}
            [tmp_out] = sess.run([model],feed_dict=feed_dict)
            disout.append(sigmoid(tmp_out[-1]))
            all_seq_lens += seq_lens
        all_out = np.concatenate(disout)
        all_seq_lens = np.cumsum(all_seq_lens)
        for p,prot in enumerate(ids):
            out = all_out[all_seq_lens[p]:all_seq_lens[p+1],:]            
            with open('outputs/'+prot+'.spotdis2.'+str(model_id),'w') as f:
                f.write("#SPOT-Disorder2 model %i's predictions for protein %s\n"%(model_id,prot))
                f.write('#\tAA\tP(D)\tO/D\n')
                for J,j in enumerate(all_seqs[p]):
                    label='O' if out[J] > 0.5 else 'D' #NOT TRUE PROBABILITY, JUST A PLACEHOLDER
                    f.write('%i\t%s\t%1.3f\t%s\n'%(J+1,j,out[J],label))



    tf.reset_default_graph()   
            #[tmp_out] = sess.run([model],feed_dict=feed_dict)
            #sigmoided_out = sigmoid(tmp_out[-1])
            #[tmp_out] = sess.run([model],feed_dict=feed_dict)
            #sigmoided_out = sigmoid(tmp_out[-1])







