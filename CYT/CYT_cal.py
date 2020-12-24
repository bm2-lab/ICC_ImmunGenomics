import pandas as pd 

data=pd.read_csv("/home/zhouchi/HCC/Expression/RelapseLiver.gene.TPM.matrix.txt",header=0,sep='\t')


sample_name_list=data.columns.values.tolist()[1:]



GZMA_exp_list=[]
PRF1_exp_list=[]
for sample in sample_name_list:
	GZMA_exp=float(data[[sample]].iloc[53986])
	PRF1_exp=float(data[[sample]].iloc[5780])
	GZMA_exp_list.append(GZMA_exp)
	PRF1_exp_list.append(PRF1_exp)






data_exp=pd.DataFrame()
data_exp["sample_name"]=sample_name_list
data_exp["GZMA_exp"]=GZMA_exp_list
data_exp["PRF1_exp"]=PRF1_exp_list

data_exp.to_csv("Relapse_CYT.txt",header=1,sep='\t',index=0)