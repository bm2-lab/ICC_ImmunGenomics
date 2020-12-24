import sys
import pandas as pd

hla_file=sys.argv[1]
sample_name=sys.argv[2]
out_dir=sys.argv[3]

HLA_dict={}
for line in open(hla_file):
	record=line.strip().split('\t')
	HLA_dict[record[0]]=[record[1],record[2],record[3],record[4],record[5],record[6]]

hla_dic={}
hla_info=pd.read_csv("/home/zhouchi/HCC/LOH/HLALOHand_docker_imagefile/hla_gen.infor.maxlen",header=-1,sep='\t')
hla_normal=hla_info.iloc[:,5]
hla_target=hla_info.iloc[:,2]
for i in range(len(hla_normal)):
	hla_dic[hla_normal[i]]=hla_target[i]

#print hla_dic
#if not os.path.exists(out_dir+'/'+sample_name):
#	os.mkdir(out_dir+'/'+sample_name)


f_out=open(out_dir+'/'+sample_name+"_hlas",'w')
for line in HLA_dict[sample_name]:
	f_out.write(hla_dic[line]+'\n')
f_out.close()


purity=open("/home/zhouchi/HCC/script/mutation_clonity_result/"+sample_name+'/'+sample_name+"_cellularity.txt").readline().strip()
ploidy=open("/home/zhouchi/HCC/script/mutation_clonity_result/"+sample_name+'/'+sample_name+"_ploidy.txt").readline().strip()

f_out_1=open(out_dir+'/'+sample_name+"_copyNumsolution.txt",'w')
f_out_1.write("Ploidy\ttumorPurity\ttumorPloidy\t\n")
f_out_1.write(sample_name+"_tumor_sorted"+"\t"+ploidy+"\t"+purity+"\t"+ploidy+"\t\n")


f_out_1.close()