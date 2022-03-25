import math
import numpy as np
import pandas as pd
import os
def cncc(df,windows_size,lebel):

	middle=int(window_size / 2)
	size=df.shape
	# print(size)
	cn=[]
	# print(size[0])
	for i in range (middle,size[0]-middle):
		# i=2

		# print(df[:3])
		k_index = range(i - middle, i + middle + 1)
		# print(k_index)
		cs=[]
		l=0
		for k in k_index:
			if k!=i :

				# print(i,k)
				
				pijpkj=df.loc[i,:]*df.loc[k,:]
				pij=df.loc[i,:]*df.loc[i,:]
				pkj=df.loc[k,:]*df.loc[k,:]
				# print(pijpkj)
				# print(pij)
				# print(pkj)
				spijpkj=pd.DataFrame.sum(pijpkj)
				spij=pd.DataFrame.sum(pij)
				spkj=pd.DataFrame.sum(pkj)

				

				if(spij== 0 or spkj==0):
					cs.append(0)	
				else:
					s=spijpkj/math.sqrt(spij*spkj)
					cs.append(spijpkj/math.sqrt(spij*spkj))	

				l=l+1
				# print(spijpkj)
				# print(spkj)
				# print(spij)
				# print()
				# cs.append()
		# print(lebel)
		cn.append(cs)
		
		# print(cn)	

	# print(cs)
	# print(df.loc[1,:])
	# print(pd.DataFrame(cn))
	file_name=output_dir+file.strip()+".csv"
	
	data_df = pd.DataFrame(cn)
	# list_of_names=['Leader', 'Time']
	data_df.columns = lebel

	# print([lebel])
	# print(file_name)
	data_df.to_csv(file_name,index=False)


	# return cn





# /home/mkabir3/Research_min/Research/all_features/Dispredict_feature/Features/DP00012

Dataset_Name_F =open("./Dataset_Name.txt" ,"r")

Dataset_Name=Dataset_Name_F.readline()
Dataset_Name=Dataset_Name.strip()


input_dir ="../../../../Dataset/"+Dataset_Name+"/PSSM_Parse/"
output_dir ="../../../../Dataset/"+Dataset_Name+"/CNCC/"

if not os.path.exists(output_dir):
	os.makedirs(output_dir)


id_list =open("../../../../Dataset/"+Dataset_Name+"/id_list.txt" ,"r")


window_size=5

lebel=[]
for i in range(window_size-1):
	# if(i!=int(window_size / 2)):
		lebel.append("cncc("+str(i)+")")
# lebel=lebel[:-2]
# print(lebel)

# Target=[ 'CNCC(1)', 'CNCC(2)', 'CNCC(3)', 'CNCC(4)', 'CNCC(5)']

for file in id_list:
	# z= id_list.readline()[:-1]

	file_loc = input_dir + file.strip() + ".csv"
	# pssm_file = open(file_loc, "r")
	print(file_loc)
	df = pd.read_csv(file_loc) 

	# print(df)
	pad=pd.DataFrame([[0]*20])


	
	pad.columns = df.columns


	# print(pad)

	for i in range(int(window_size/2)):
		df=pd.concat([pad,df,pad],axis=0)
	
	df = df.reset_index(drop=True)
	# print(df)

	cncc(df,window_size,lebel)


	