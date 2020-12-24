###filter somatic mutation
import sys
import os

mutation_file_in=sys.argv[1]
sample_name=sys.argv[2]
out_dir=sys.argv[3]


flag=0
for line in open(mutation_file_in):
	tumor_sample=""
	#head=""
	if line.startswith("##tumor_sample"):
		tumor_sample=line.strip().split("=")[1]
		break
	else:
		continue
#print tumor_sample

for line in open(mutation_file_in):
	head=""
	#head=""
	if line.startswith("#CHROM"):
		head=line.strip().split("\t")[10]
		break
	else:
		continue
#print head


if tumor_sample==head:
	flag=1
	#print flag
if flag==1:
	cmd_filter="bcftools filter -m + -i \'(FMT/AF[1:0] >= 0.01 && SUM(FMT/AD[1:])>=30 && FMT/AD[1:1] >=3 && FMT/AD[0:1] <=3)\' " + mutation_file_in + " -o " + out_dir + "/" + sample_name + "_filter.vcf -O v"
	#print cmd_filter
	os.system(cmd_filter)
else:
	cmd_filter="bcftools filter -m + -i \'(FMT/AF[0:0] >= 0.01 && SUM(FMT/AD[0:])>=30 && FMT/AD[0:1] >=3 && FMT/AD[1:1] <=3)\' " + mutation_file_in + " -o " + out_dir + "/" + sample_name + "_filter.vcf -O v"
	#print cmd_filter
	os.system(cmd_filter)