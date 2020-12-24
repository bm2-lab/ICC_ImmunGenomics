import pandas as pd 

import numpy as np

data=pd.read_csv("/home/zhouchi/HCC/Immunoscore/RF5Y/RF5Y_input.txt",header=0,sep='\t')


medain_CD8=float(np.median(data[["T cells CD8"]]))

median_CD3=float(np.median(data[["CD3 T cell"]]))


immunoscore_list=[]

for i in range(len(data[["T cells CD8"]])):
	cd8=float(data[["T cells CD8"]].iloc[i])
	cd3=float(data[["CD3 T cell"]].iloc[i])
	if cd8<medain_CD8 and cd3<median_CD3:
		immunoscore=0
	elif cd8>=medain_CD8 and cd3<median_CD3:
		immunoscore=1
	elif cd8<medain_CD8 and cd3>=median_CD3:
		immunoscore=2
	else:
		immunoscore=3
	immunoscore_list.append(immunoscore)

data["Immunoscore"]=immunoscore_list


data.to_csv("/home/zhouchi/HCC/Immunoscore/RF5Y/RF5Y_immunoscore.txt",header=1,sep='\t',index=0)
