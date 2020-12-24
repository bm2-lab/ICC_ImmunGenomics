sample_name="M20180364"

normal_sample_name="M20180365"

bwa mem -t 64 /home/zhouchi/HCC/database/ucsc.hg19.fasta /home/zhouchi/HCC/data/5year_wufufa_66_cleandata/${sample_name}/${sample_name}_R1.trim.fastq.gz /home/zhouchi/HCC/data/5year_wufufa_66_cleandata/${sample_name}/${sample_name}_R2.trim.fastq.gz | samtools sort -@ 4 -m 2G - -o ${sample_name}_tumor_sorted.bam



bwa mem -t 64 /home/zhouchi/HCC/database/ucsc.hg19.fasta /home/zhouchi/HCC/data/5year_wufufa_66_cleandata/${normal_sample_name}/${normal_sample_name}_R1.trim.fastq.gz /home/zhouchi/HCC/data/5year_wufufa_66_cleandata/${normal_sample_name}/${normal_sample_name}_R2.trim.fastq.gz | samtools sort -@ 4 -m 2G - -o ${sample_name}_normal_sorted.bam


samtools index ${sample_name}_tumor_sorted.bam
samtools index ${sample_name}_normal_sorted.bam

docker run -v /home/zhouchi/HCC/sCNA/:/data quay.io/jeltje/varscan2 -c ${sample_name}_normal_sorted.bam -t ${sample_name}_tumor_sorted.bam -q ${sample_name} -i ucsc.hg19.fasta -b centromeres.bed -w targets.bed -s tmpdir > cnv_result/${sample_name}.cnv 

python scna_filter.py ${sample_name}


if [ -f final_result/${sample_name}_final.cnv ];then
	rm -rf ${sample_name}_tumor_sorted.bam*
	rm -rf ${sample_name}_normal_sorted.bam*
fi

docker rm $(docker ps -q -f status=exited)