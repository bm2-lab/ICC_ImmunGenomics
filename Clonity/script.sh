sample_name="M20180316"
if [ ! -f /home/zhouchi/HCC/data/5year_sampleout/${sample_name}/somatic_oncefiltered.vcf ];then
	gunzip /home/zhouchi/HCC/data/5year_sampleout/${sample_name}/somatic_oncefiltered.vcf.gz
fi
python mutation_filter.py /home/zhouchi/HCC/data/5year_sampleout/${sample_name}/somatic_oncefiltered.vcf ${sample_name} /home/zhouchi/HCC/data/5year_sampleout/${sample_name}
#bcftools filter -m + -i '(FMT/AF[0:0] >= 0.05 && SUM(FMT/AD[0:])>=30 && FMT/AF[1:0] <= 0.05)' /home/zhouchi/HCC/data/5year_sampleout/${sample_name}/somatic_oncefiltered.vcf -o /home/zhouchi/HCC/data/5year_sampleout/${sample_name}_fiter.vcf -O v
#gunzip /home/zhouchi/HCC/data/5year_sampleout/${sample_name}/somatic_oncefiltered.vcf.gz
bash mutation_clonity_pipeline.sh \
/home/zhouchi/HCC/data/5year_sampleout/${sample_name}/tumor_realigned.bam \
/home/zhouchi/HCC/data/5year_sampleout/${sample_name}/normal_realigned.bam \
/home/zhouchi/HCC/data/5year_sampleout/${sample_name}/${sample_name}_filter.vcf \
/home/zhouchi/HCC/database/human_g1k_v37.fasta \
/home/zhouchi/HCC/database/hg37.gc50Base.wig.gz \
${sample_name} \
/home/zhouchi/HCC/script/mutation_clonity_result/${sample_name}
