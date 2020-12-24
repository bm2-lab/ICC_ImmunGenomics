import os
import pandas as pd

sample_list=os.listdir("/home/zhouchi/HCC/neoantigen/Relapse_free_5y")


hla_dic={}
hla_info=pd.read_csv("/home/zhouchi/HCC/LOH/HLALOHand_docker_imagefile/hla_gen.infor.maxlen",header=-1,sep='\t')
hla_normal=hla_info.iloc[:,5]
hla_target=hla_info.iloc[:,2]
for i in range(len(hla_normal)):
	hla_dic[hla_target[i]]=hla_normal[i]


HLA_dict={}
for line in open("RF5Y_HLA_consensus.txt"):
	record=line.strip().split('\t')
	HLA_dict[record[0]]=[record[1],record[2],record[3],record[4],record[5],record[6]]

f_o=open("HLA_summary.txt",'w')
f_o.write("sample_name\tLOH_status\tKeptAllele\tLossAlle\n")
for sample in sample_list:
	allele_all=HLA_dict[sample]
	if not os.path.exists("/home/zhouchi/HCC/HLALOH/"+sample+'/'+sample+".30.DNA.HLAlossPrediction_CI.xls"):
		allele_kept_str=','.join(allele_all)
		f_o.write(sample+'\t'+"No"+'\t'+allele_kept_str+'\t'+"NA"+'\n')
	else:
		loss_allele=[]
		data=pd.read_table("/home/zhouchi/HCC/HLALOH/"+sample+'/'+sample+".30.DNA.HLAlossPrediction_CI.xls")
		for i in range(len(data.HLA_A_type1)):
			if ((data.HLA_type1copyNum_withBAFBin[i] < 0.5 or data.HLA_type2copyNum_withBAFBin[i] < 0.5) and data.PVal_unique[i] < 0.05 and (data.LossAllele[i] not in loss_allele)):
				loss_allele.append(data.LossAllele[i])
			else:
				continue
		loss_allele_short=[]
		for allele in loss_allele:
			loss_allele_short.append(hla_dic[allele])
		allele_kept=allele_all
		for i in loss_allele_short:
			allele_kept.remove(i)
		allele_kept_str=','.join(allele_kept)
		allele_loss_str=','.join(loss_allele_short)
		LOH_status="No"
		if len(loss_allele_short)!=0:
			#print len(loss_allele_short)
			LOH_status="Yes"
		f_o.write(sample+'\t'+LOH_status+'\t'+allele_kept_str+'\t'+allele_loss_str+'\n')
f_o.close()


