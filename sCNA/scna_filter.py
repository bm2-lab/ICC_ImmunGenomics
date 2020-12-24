import pandas as pd 
import sys

sample_name=sys.argv[1]


data=pd.read_csv("cnv_result/"+sample_name+".cnv",header=0,sep='\t')

data_filter=data[(abs(data["seg.mean"])>0.25)]

data_filter.to_csv("final_result/"+sample_name+"_final.cnv",header=1,sep='\t',index=0)