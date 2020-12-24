import os
import sys

sample_name=sys.argv[1]


f_dNdS_input=open("/home/zhouchi/HCC/dNdS/Relapse/"+sample_name+'/'+sample_name+"_dNdS_input.txt",'w')
f_dNdS_input.write("sampleID\tchr\tpos\tref\tmut\n")
pos_list=[]
for line in open("/home/zhouchi/HCC/dNdS/Relapse/"+sample_name+'/'+sample_name+"_vep_ann_all.txt"):
	if line.startswith("#"):
		continue
	else:
		record=line.strip().split('\t')
		if record[0] not in pos_list:
			f_dNdS_input.write(sample_name+"\t"+record[0].split('_')[0]+"\t"+record[0].split('_')[1]+"\t"+record[0].split('_')[2].split('/')[0]+"\t"+record[0].split('_')[2].split('/')[1]+"\n")
			pos_list.append(record[0])
		else:
			continue
f_dNdS_input.close()
