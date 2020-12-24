import os
import sys

sample_name=sys.argv[1]



if not os.path.exists("/home/zhouchi/HCC/dNdS/Relapse/"+sample_name):
	os.mkdir("/home/zhouchi/HCC/dNdS/Relapse/"+sample_name)


f_dnds=open("/home/zhouchi/HCC/dNdS/Relapse/"+sample_name+'/'+sample_name+"_pass.vcf","w")
for line in open("/home_2/zhouchi/HCC_only_intrahepatic/mutation_result/"+sample_name+'/'+sample_name+"_filter.vcf"):
	if line.startswith("#"):
		f_dnds.write(line)
	elif line.strip().split('\t')[6]=="PASS":
		f_dnds.write(line)
	else:
		continue
f_dnds.close()