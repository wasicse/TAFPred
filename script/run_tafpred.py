# %%
# Author: Md Wasi Ul Kabir
# Date: December 28 2020
import os
import pandas as pd
import subprocess
import pathlib 
import joblib
import pandas as pd
import numpy as np
import docker
from pathlib import Path
from optparse import OptionParser

# check if any files are missing
def checkmissing(source_dir):
    if not os.path.exists(source_dir):
        print("Missing", source_dir )
        missing_list.append(pid)
        
# check if datafarme has NAN values    
def checkforNAN(df,seqlength):
    if df.isnull().values.any():
        print("NAN")
        print(df[df.isnull().any(axis=1)])     
   
    if df.shape[0]!=seqlength:
        print("Length Mismatch", df.shape[0], seqlength)

# merge all features into one file     
def mergedata(PrintP,label):
    id_list =open("../Features/id_list.txt" ,"r")
    output_dir ="../output/merge_features/"
    output_file_end=".csv"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    fasta_dir="../Features/FASTA/"
        
    missing_list=[]
    flag=0
    for pid in id_list:
        # print(pid)
        pid=pid.strip()
        output_file_name=output_dir+pid+output_file_end
        
        read_fasta= open(fasta_dir+pid+".fasta", "r")
        fasta=read_fasta.readline()
        fasta=read_fasta.readline()
        # print(type(fasta))
        seqlength=fasta[:-1].__len__()
        print(pid, seqlength)   
        #Dispredict Features      

        source_dir="../Features/Features/DisPredict/Features/"+pid+"/"+pid+".57pfeatures"                     
        dispredict_column=['O/D(1)', 'AA(1)', 'PP(1)', 'PP(2)', 'PP(3)', 'PP(4)', 'PP(5)', 'PP(6)', 'PP(7)', 'PSSM(1)', 'PSSM(2)', 'PSSM(3)', 'PSSM(4)', 'PSSM(5)', 'PSSM(6)', 'PSSM(7)', 'PSSM(8)', 'PSSM(9)', 'PSSM(10)', 'PSSM(11)', 'PSSM(12)', 'PSSM(13)', 'PSSM(14)', 'PSSM(15)', 'PSSM(16)', 'PSSM(17)', 'PSSM(18)', 'PSSM(19)', 'PSSM(20)', 'SS(1)', 'SS(2)', 'SS(3)', 'ASA(1)', 'dphi(1)', 'dpsi(1)', 'MG(1)', 'BG(1)', 'BG(2)', 'BG(3)', 'BG(4)', 'BG(5)', 'BG(6)', 'BG(7)', 'BG(8)', 'BG(9)', 'BG(10)', 'BG(11)', 'BG(12)', 'BG(13)', 'BG(14)', 'BG(15)', 'BG(16)', 'BG(17)', 'BG(18)', 'BG(19)', 'BG(20)', 'sPSEE(1)', 't(1)']
        df1=pd.read_csv(source_dir,delim_whitespace=True,header=None,skiprows=1 )
        df1.columns=dispredict_column
        drop_col= ['O/D(1)']
        df1=df1.drop(drop_col, axis=1)
        df1.columns=["Dis_"+ s for s in  df1.columns]
        checkforNAN(df1,seqlength)
        if PrintP: print("Dispredict", df1.shape)  
                
        # #PSSM Features    
        source_dir="../Features/Features/PSSM_Parse/"+pid+".csv"  

        df3=pd.read_csv(source_dir )
        df3.columns=["PSI_"+ s for s in  df3.columns]
        checkforNAN(df3,seqlength)
        if PrintP: print("PSSM",df3.shape)
        
        #Spider Features
        source_dir="../Features/Features/Spider/"+pid.strip()+ ".i1"

        df=pd.read_csv(source_dir,delim_whitespace=True)  
        df4=df.drop([ "#","AA","SS","SS8"], axis=1)
        df4.columns=["Spi_"+ s for s in  df4.columns]
        checkforNAN(df4,seqlength)
        if PrintP: print("Spider",df4.shape)

        
        # #OPAL Features
        source_dir="../Features/Features/OPAL/"+pid+".txt"
        df=pd.read_csv(source_dir,delim_whitespace=True )   
        df5=df.drop([ "No:","residues"], axis=1)
        df5.columns=["Opa_"+ s for s in  df5.columns]
        checkforNAN(df5,seqlength)
        if PrintP: print("OPAL",df5.shape)

        # #CNCC Features
        source_dir="../Features/Features/CNCC/"+pid+".csv"
        df6=pd.read_csv(source_dir)  
        checkforNAN(df6,seqlength)
        if PrintP: print("CNCC",df6.shape)
        
        
        #SpotDisorder
        source_dir="../Features/Features/SpotDisorder/"+pid+".csv"   
        col= (   
            r"ASA,HSEa-u,HSEa-d,CN13,theta,tau,phi,psi,theta_c,tau_c,phi_c,psi_c,P(3-C),P(3-E),P(3-H),P(8-C),P(8-S),P(8-T),P(8-H),"
            r"P(8-G),P(8-I),P(8-E),P(8-B),HH_1,HH_2,HH_3,HH_4,HH_5,HH_6,HH_7,HH_8,HH_9,HH_10,HH_11,HH_12,HH_13,HH_14,HH_15,HH_16,"
            r"HH_17,HH_18,HH_19,HH_20,HH_21,HH_22,HH_23,HH_24,HH_25,HH_26,HH_27,HH_28,HH_29,HH_30,PSSMS_1,PSSMS_2,PSSMS_3,PSSMS_4,"
            r"PSSMS_5,PSSMS_6,PSSMS_7,PSSMS_8,PSSMS_9,PSSMS_10,PSSMS_11,PSSMS_12,PSSMS_13,PSSMS_14,PSSMS_15,PSSMS_16,PSSMS_17,PSSMS_18,PSSMS_19,PSSMS_20"
            )
        SpotDisorder_column=col.split(",")  
        SpotDisorder_column=[x.strip() for x in SpotDisorder_column ]
        df9=pd.read_csv(source_dir,index_col=0 )
        df9.columns=SpotDisorder_column
        checkforNAN(df9,seqlength)        
        if PrintP: print("SpotDisorder",df9.shape)  
        
        #SpotDisorder Probability
        source_dir="../Features/Features/SpotDisorderProba/"+pid+".spotd2"    
        df=pd.read_csv(source_dir,delim_whitespace=True,skiprows=1)  
        df10=df.drop([ "#","AA","O/D"], axis=1)
        df10.columns=["Spot_"+ s for s in  df10.columns]
        checkforNAN(df10,seqlength)
        if PrintP: print("SpotDisorderProba",df10.shape)

    # Target (Set Dummy Target for windowing code)    
        df0=pd.DataFrame()
        df0["Target"]=  range(df10.shape[0])
        print("Target",df0.shape)

    # code to add target
        # target_dir="../Feature_Extraction/Dataset/example/"+Angle+"target/"
        # target_dir_file_end=".csv"
        # source_dir=target_dir+pid+target_dir_file_end
        # df0=pd.read_csv(source_dir,index_col=0)   #,header=None
    
        # # df0.columns=["Target"]
        # if PrintP: print("Target",df0.shape)

        listdf=[df0,df10,df4,df3,df5,df6,df9,df1,df4[["Spi_Phi","Spi_Psi"]]]
        merged = pd.concat(listdf, axis=1)
        merged = merged.loc[:, ~merged.columns.str.contains('^Unnamed')] 
        if PrintP: print("Merge",merged.shape)  


        for dff in listdf:      
            if(dff.isnull().values.any()):
                prinf(dff)
                print("NAN")
                break;
            
        if(merged.isnull().values.any()):
            print("Merge files have NAN")
            
        merged.to_csv(output_dir+"/"+pid+'.csv',index=False,header=label) 
  
    return merged.shape

# run windowing code
def windowing(len,startwindow_size,endwindow_size):
    
    windowintput=open("./window_param.txt","w")    
    windowintput.write(len+ '\n')
    windowintput.write(str(startwindow_size)+ '\n')
    windowintput.write(str(endwindow_size)+ '\n')
    windowintput.close()
    output_dir ="../output/Windowed_file/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bashCommand='javac run_windowing.java '
    output = subprocess.check_output(['bash','-c', bashCommand])
    bashCommand='java run_windowing'
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 


def getprediction():
    np.set_printoptions(precision=3)
    workspace="../output"
    pathlib.Path(workspace).mkdir(parents=True, exist_ok=True) 

    X_test=pd.read_csv("../output/Windowed_file/feat_179_w_3.csv",header=None).to_numpy()
    X_test=X_test[:,1:]

    feature=open("../models/OPTGeneticFeat.txt","r").readline()
    mask=pd.DataFrame(list(feature)).astype(int)[0].to_numpy().astype(bool)
    fmask=[]   
    for i in range(3):
        fmask.extend(mask)

    X_test = X_test[:, fmask]


    scaler= joblib.load("../models/phi_scaler.pkl")
    X_test = scaler.transform(X_test)


    phisaved_model = joblib.load("../models/phi_model.pkl")
    phiproba = phisaved_model.predict(X_test)

    psisaved_model = joblib.load("../models/psi_model.pkl")
    psiproba = psisaved_model.predict(X_test)
    print(phiproba.shape)

    id_list =open("../Dataset/example/id_list.txt" ,"r")
    start=0
    for pid in id_list:
        pid=pid.strip()
        fastafile=open("../Dataset/example/FASTA/"+pid+".fasta","r")
        fasta=fastafile.readline().strip()
        fasta=fastafile.readline().strip()    
        end=start+fasta.__len__()
   
        print("S",start,"E",end)
        fphiproba=phiproba[start:end]
        print(fphiproba.shape)
        fpsiproba=psiproba[start:end]    
        result=np.hstack((  np.round(np.arange(1,fasta.__len__()+1).reshape(-1,1),0) ,np.round(fphiproba, 3).reshape(-1,1) ,np.round(fpsiproba, 3).reshape(-1,1) )) 
        with open("../output/"+pid+"_result.txt", "w") as f:
            f.write("Residue No\tΔPhi\tΔPsi\n")
        with open("../output/"+pid+"_result.txt", "ab") as f:
            np.savetxt(f,  result, delimiter='\t\t',fmt="%s ")        
        start=fasta.__len__()

# check if a container is running or not
def checkcontainerstatus(containername):
    cli = docker.APIClient()
    try:
        inspect_dict = cli.inspect_container(containername)
        state = inspect_dict['State']
        print(state)
        is_running = state['Status'] == 'running'

        if is_running:
            print("My container is running!")
            return False
    except:
        print("My container is not running!")
        return True

# collect features from docker container
def collectfeatures(containername):


    print("Pulling docker image")
    bashCommand="docker pull wasicse/featureextract:1.0"
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 

    print("Checking container status")
    bashCommand="docker ps -q -f name={"+containername+"}"
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 
    if checkcontainerstatus(containername):
        print("Creating container")
        bashCommand="docker run -itd -v "+os.getcwd()+"/Databases:/opt/common --name "+containername+"  wasicse/featureextract:1.0"
        print(bashCommand)
        output = subprocess.check_output(['bash','-c', bashCommand])
        print(output.decode('utf-8')) 
    else:
        print("Starting container")
        bashCommand="docker start "+containername
        print(bashCommand)
        output = subprocess.check_output(['bash','-c', bashCommand])
        print(output.decode('utf-8'))    

    print("Copying files to container")
    bashCommand="docker cp ../Dataset/example/ "+containername+":/opt/FeatureExtractionDocker/Dataset"
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 

    print("Removing old features")
    bashCommand="rm -rf ../Features/"
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 

    print("Running feature extraction. It might take a while to complete. and depends on the number of sequences in the input file.")
    bashCommand= "docker exec "+containername+" bash -c \"cd /opt/FeatureExtractionDocker/FeatureExtractTool ; /opt/.pyenv/versions/miniconda3-3.9-4.10.3/envs/feat/bin/python runScript.py\""
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 
    
    print("copying features from container")
    bashCommand="docker cp "+containername+":/opt/FeatureExtractionDocker/Dataset/example/ ../Features/"
    output = subprocess.check_output(['bash','-c', bashCommand])
    print(output.decode('utf-8')) 

    # print("Stop container")
    # bashCommand="docker stop "+containername
    # print(bashCommand)
    # output = subprocess.check_output(['bash','-c', bashCommand])
    # print(output.decode('utf-8'))   


if __name__ == "__main__":  
    parent_path = str(Path(__file__).resolve().parents[1])
    print("Parent Dir",parent_path)

    parser = OptionParser()
    parser.add_option("-f", "--containerName", dest="containerName", help="Container Name.", default="taffeaturesv1")
    parser.add_option("-o", "--output_path", dest="output_path", help="Path to output.", default=parent_path+'/output/')

    (options, args) = parser.parse_args()

    print("Container Name:",options.containerName)
    print("Output Path:",options.output_path)

    #Print features options
    PrintP=True
    #Add label to features
    label=True

    # Collect features by running the docker container
    collectfeatures(options.containerName)

    # Merge features
    merge_shape=mergedata(PrintP,label)

    # Windowing 
    startwindow_size=3
    endwindow_size=3
    windowing(str(merge_shape[1]-1),startwindow_size,endwindow_size)

    # Get prediction from the windowed features
    getprediction()
    


