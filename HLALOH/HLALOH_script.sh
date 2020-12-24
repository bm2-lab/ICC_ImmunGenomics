#bwa mem -t 64 /home/zhouchi/HCC/database/ucsc.hg19.fasta /home/zhouchi/HCC/test_data/M20180293/M20180293_R1.trim.fastq.gz /home/zhouchi/HCC/test_data/M20180293/M20180293_R2.trim.fastq.gz | samtools sort -@ 4 -m 2G - -o M20180293_tumor_sorted.bam



#bwa mem -t 64 /home/zhouchi/HCC/database/ucsc.hg19.fasta /home/zhouchi/HCC/test_data/M20180294/M20180294_R1.trim.fastq.gz /home/zhouchi/HCC/test_data/M20180294/M20180294_R2.trim.fastq.gz | samtools sort -@ 4 -m 2G - -o M20180293_GL_sorted.bam


#samtools index M20180293_tumor_sorted.bam
#samtools index M20180293_GL_sorted.bam

#samtools view -h /home/zhouchi/HCC/test_data/result/M20180293/normal_realigned.bam
#samtools view /home/zhouchi/HCC/test_data/result/M20180293/normal_realigned.bam | awk '$3="chr"$3' | awk -F ' ' '$7=($7=="=" || $7=="*"?$7:sprintf("chr%s",$7))' | tr ' ' '\t' | samtools view -b >> M20180293_GL_sorted.bam




sample_id="M20180336"

if [ ! -d /home/zhouchi/HCC/HLALOH/${sample_id} ];then
	mkdir -p /home/zhouchi/HCC/HLALOH/${sample_id}
fi

samtools view -H /home/zhouchi/HCC/data/5year_sampleout/${sample_id}/normal_realigned.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - /home/zhouchi/HCC/data/5year_sampleout/${sample_id}/normal_realigned.bam > /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_GL_sorted.bam

samtools view -H /home/zhouchi/HCC/data/5year_sampleout/${sample_id}/tumor_realigned.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - /home/zhouchi/HCC/data/5year_sampleout/${sample_id}/tumor_realigned.bam > /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_tumor_sorted.bam



samtools index /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_GL_sorted.bam
samtools index /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_tumor_sorted.bam

python generate_hla_cn.py RF5Y_HLA_consensus.txt ${sample_id} /home/zhouchi/HCC/HLALOH/${sample_id}

cd ${sample_id}

bash /home/zhouchi/HCC/LOH/HLALOHand_docker_imagefile/test.sh ${sample_id}_tumor_sorted.bam ${sample_id}_GL_sorted.bam ${sample_id}_copyNumsolution.txt ${sample_id}_hlas result

sudo cp result/${sample_id}.30.DNA.HLAlossPrediction_CI.xls .
sudo rm -rf result

cd ..
rm -rf /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_GL_sorted.bam /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_tumor_sorted.bam /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_GL_sorted.bam.bai /home/zhouchi/HCC/HLALOH/${sample_id}/${sample_id}_tumor_sorted.bam.bai

