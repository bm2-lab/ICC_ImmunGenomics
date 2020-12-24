#java -jar /home/zhouchi/pTuneos/software/trimmomatic-0.36.jar PE -threads 32 -phred33 /home_2/zhouchi/HCC_RNAseq/Norecurrencein5years_sample/R18027382-20180426RNA-POOL-01-M20180295-R_R1.fastq.gz /home_2/zhouchi/HCC_RNAseq/Norecurrencein5years_sample/R18027382-20180426RNA-POOL-01-M20180295-R_R2.fastq.gz -baseout  R18027382-20180426RNA-POOL-01-M20180295_clean.fq.gz ILLUMINACLIP:/home/zhouchi/iTunes/software/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#hisat2 -p 64 -x /home/zhouchi/HCC/DEG/database/genome/genome -1 ../R18027382-20180426RNA-POOL-01-M20180295_clean_1P.fq.gz -2 ../R18027382-20180426RNA-POOL-01-M20180295_clean_2P.fq.gz -S tumor.sam


#java -jar /home/zhouchi/pTuneos/software/trimmomatic-0.36.jar PE -threads 32 -phred33 /home_2/zhouchi/HCC_RNAseq/Norecurrencein5years_sample/R18027382-20180426RNA-POOL-01-M20180296-R_R1.fastq.gz /home_2/zhouchi/HCC_RNAseq/Norecurrencein5years_sample/R18027382-20180426RNA-POOL-01-M20180296-R_R2.fastq.gz -baseout  R18027382-20180426RNA-POOL-01-M20180296_clean.fq.gz ILLUMINACLIP:/home/zhouchi/iTunes/software/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#hisat2 -p 64 -x /home/zhouchi/HCC/DEG/database/genome/genome -1 ../R18027382-20180426RNA-POOL-01-M20180296_clean_1P.fq.gz -2 ../R18027382-20180426RNA-POOL-01-M20180296_clean_2P.fq.gz -S normal.sam


#gfold count -ann /home/zhouchi/HCC/DEG/database/genome.gtf -tag tumor.sam -o tumor.read_cnt
#gfold count -ann /home/zhouchi/HCC/DEG/database/genome.gtf -tag normal.sam -o normal.read_cnt
gfold diff -s1 tumor -s2 normal -suf .read_cnt -o tumorVSnormal.diff
