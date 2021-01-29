#!/usr/local/software/anaconda3/bin/python
import os
import logging
import glob
import time
from multiprocessing import Pool
import pandas as pd
from functools import wraps
import subprocess
import traceback
import pandas as pd
import opt


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwds):
        t0 = time.time()
        with open(log_path, 'a') as f:
            f.write('patient_{}:start {}\n'.format(args[0], func.__name__))
        result = func(*args, **kwds)
        with open(log_path, 'a') as f:
            t1 = time.time()
            f.write('patient_{}: {} finished! cost {:0.3f} seconds!\n'.format(
                args[0], func.__name__, t1-t0))
        return result
    return wrapper


@timer
def clean(patid, tmp_files):
    '''remove trash files to release disk space'''
    while tmp_files:
        trash = tmp_files.pop(0)
        os.remove(trash)
        with open(log_path, 'a') as f:
            f.write('files {} will be removed'.format(trash))


@timer
def trimmomatic(patid, settings):
    '''
    input: .fastq files(two files);
    output: accession_1U/1P/2U/2P(four files)
    '''
    cmd = 'java -jar {trimmomatic_path} PE -threads 5 {fastq_1} {fastq_2} -baseout ./tumor ILLUMINACLIP:{adapters}:2:30:12:1:true LEADING:3 TRAILING:3 MINLEN:36'.format(**settings)

    try:
        with open(log_path, 'a') as f:
            f.write(cmd)
            p = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
            p.wait()
    except Exception as e:
        print(e)
    # 查看Unpaired文件大小，如果发现>100MB就记录此异常
    unpaired_1 = glob.glob('./*_1U')[0]
    unpaired_2 = glob.glob('./*_2U')[0]
    with open(log_path, 'a') as f:
        if os.path.getsize(unpaired_1)/1024**2 > 100:
            f.write("Unpaired file {} too big!\n".format(unpaired_1))
        elif os.path.getsize(unpaired_2)/1024**2 > 100:
            f.write("Unpaired file {} too big!\n".format(unpaired_2))
        else:
            f.write("Pass trimmomatic")


@timer
def star(patid, settings):
    '''
    make index first; then mapping
    input: accession_1P/2P(two files)
    output: five files; only need is accessionAligned.sortedByCoord.out.bam
    '''
    # STAR --runThreadN 4 --genomeDir /home/mengflz/HPD/Gide_2019/index/ --outSAMtype BAM SortedByCoordinate --readFilesIn "/home/mengflz/HPD/Zhao_2019/test/SRR8281222_1P" "/home/mengflz/HPD/Zhao_2019/test/SRR8281222_2P" --outFileNamePrefix /home/mengflz/HPD/Zhao_2019/test/SRR8281222
    
    settings["filein_1"] = glob.glob('./*_1P')[0]
    settings["filein_2"] = glob.glob('./*_2P')[0]

    cmd = 'STAR --runThreadN {nt} --genomeDir {indexpath} --outSAMtype BAM SortedByCoordinate --readFilesIn {filein_1} {filein_2} --outFileNamePrefix ./{patid}'.format(
        **settings)
    with open(log_path, 'a') as f:
        f.write(cmd)
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()


@timer
def stringtie(patid, settings):
    '''
    input: accessionAligned.sortedByCoord.out.bam;
    output: .gtf file as transcripts, .tsv as expression matrix
    '''
    settings['filein'] = glob.glob('./*.bam')[0]
    cmd = '/home/mengflz/software/stringtie-2.1.4.Linux_x86_64/stringtie {filein} -p 10 -G {gtfpath} -o ./test/tumor.gtf -A ./test/tumor_expression.tsv'.format(
        **settings)
    with open(log_path, 'a') as f:
        f.write(cmd)
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()
    

@timer
def featureCounts(patid, settings):
    '''
    -t Specify the feature type
    -p isPairedEnd If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads.
    -B requireBothEndsMapped

    '''
    settings['filein'] = glob.glob('./*.bam')[0]
    # settings['filein_normal'] = glob.glob('./normal*.bam')[0]
    # cmd_baseline_normal = "{featureCounts_path} -t exon -g gene_id -Q 20 -p -B -T 10 -a {gtfpath} -o featureCounts_raw_normal.counts {filein_normal}".format(**settings)
    cmd = "{featureCounts_path} -t exon -g gene_id -Q 20 -p -B -T 10 -a {gtfpath} -o featureCounts_raw.counts {filein}".format(**settings)
    with open(log_path, 'a') as f:
        f.write(cmd)
        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()
    # cmd = "{featureCounts_path} --ignoreDup -Q 20 -B -p -T 16 -g gene_id -a {} -o STAR.raw.Counts  Aligned.sortedByCoord.out.bam.format(gtf)"


@timer
def cal_count(caseid, settings):
    counts_matrix = pd.read_csv("./featureCounts_raw.counts", sep='\t', header=1, index_col=0)
    gene_length = counts_matrix.loc[:, 'Length']/1000
    raw_count = counts_matrix.iloc[:, 5]
    FPKM = (raw_count*10**6/raw_count.sum())/gene_length
    TPM = (raw_count/gene_length)*10**6/(raw_count/gene_length).sum()
    output = pd.concat([TPM, FPKM], axis=1, keys=['TPM', 'FPKM'])
    output.to_csv('./cal_exp.tsv', sep='\t')


@timer
def trash_remove(caseid, settings):
    '''
    flag = 1, twice sequencing
    '''
    # need to return a list recording if got wrong
    # remain_file = ['output_filtered.vcf', 'output-tnhaplotyper.vcf.gz', 'output-tnhaplotyper.vcf','output-tnhaplotyper.vcf.gz.tbi',
    #                'normal_realigned.bam', 'tumor_realigned.bam', 'normal_realigned.bam.bai', 'tumor_realigned.bam.bai']
    remain_file = ['{}Aligned.sortedByCoord.out.bam'.format(caseid),
                   'featureCounts_raw.counts', 'cal_exp.tsv']
    file_total = glob.glob('*')
    with open(log_path, 'a') as f:
        for i in file_total:
            if i not in remain_file:
                os.remove(i)
                f.write('{} remove file {}\n'.format(caseid, i))
        for j in remain_file:
            if j not in file_total:
                f.write('{} WARNINGS: no needed file {}\n'.format(caseid, j))


def process(patid):
    # define a temporary list to record files need will be removed
    # 在执行完函数后删除它的输入文件
    try:
        with open(log_path, 'a') as f:
            f.write('***start handle {}***\n'.format(patid))
        if not os.path.exists('./{}'.format(patid)):
            os.mkdir('./{}'.format(patid))
        os.chdir('./{}'.format(patid))

        tissue_list = [item.split('/')[-1] for item in glob.glob(file_path + patid + '/WTS/*')]
        for type in tissue_list:
            if not os.path.exists('./{}'.format(type)):
                os.mkdir('./{}'.format(type))
            os.chdir('./{}'.format(type))

            settings = opt.wts_options(patid, type, nt=20).__dict__

            trimmomatic(patid, settings)
            star(patid, settings)
            featureCounts(patid, settings)  # stringtie(patid, settings)
            cal_count(patid, settings)
            trash_remove(patid, settings)
            os.chdir(os.path.dirname(os.getcwd()))

        os.chdir(os.path.dirname(os.getcwd()))
    except Exception:
        traceback.print_exc()


def TCGA_process(patid):
    with open(log_path, 'a') as f:
        f.write('***start handle {}***\n'.format(patid))
    if not os.path.exists('./{}'.format(patid)):
        os.mkdir('./{}'.format(patid))
    os.chdir('./{}'.format(patid))
    try:
        settings = opt.wts_options(patid, nt=20).__dict__
        featureCounts(patid, settings)  # stringtie(patid, settings)
        cal_count(patid, settings)
        trash_remove(patid, settings)
    except Exception as e:
        with open(log_path, 'a') as f:
            f.write(traceback.format_exc())
    finally:
        os.chdir(os.path.dirname(os.getcwd()))

@timer
def merge_multi_files(list_total):
    output_dataframe_fpkm = pd.DataFrame()
    output_dict = {}
    gene_id = output_dataframe_fpkm.index

    for i in list_total:
        tmp_frame = pd.read_csv(
            './{}/{}.tsv'.format(i, i), sep='\t', index_col=0, low_memory=False)
        tmp_frame = tmp_frame[~tmp_frame.index.str.match('STRG')]
        tmp_frame = tmp_frame.groupby(level=0).mean()
        output_dict[i] = tmp_frame
        gene_id = tmp_frame.index.union(gene_id)
    # gene_id,output_dict 获得了基因索引名的总和
    output_dataframe_fpkm = pd.DataFrame(index=gene_id)
    output_dataframe_tpm = pd.DataFrame(index=gene_id)
    for i in output_dict.keys():
        tmp_dataframe = output_dict[i].loc[gene_id, :]
        tmp_dataframe_seris_tpm = tmp_dataframe['TPM']
        tmp_dataframe_seris_fpkm = tmp_dataframe['FPKM']
        output_dataframe_tpm.insert(0, i, tmp_dataframe_seris_tpm)
        output_dataframe_fpkm.insert(0, i, tmp_dataframe_seris_fpkm)

    output_dataframe_fpkm.to_csv('./expression_FPKM.csv', sep=',')
    output_dataframe_tpm.to_csv('./expression_TPM.csv', sep=',')


def main():
    global log_path, file_path
    log_path = os.getcwd() + "/pipeline_metastasis.log"
    file_path = "/home_2/wzt/metastasis/"  # sys.argv[1]

    with open(log_path, 'a') as f:
        f.write("*********" + time.strftime("%a %b %d %H:%M:%S %Y",
                                            time.localtime()) + "**********" + '\n')
    # 创建结果路径
    output_path = os.path.normpath(os.getcwd()+'/output_metastasis')
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    os.chdir(output_path)

    # start build task queue
    file_list = glob.glob(file_path + '*')  # /home_2/wzt/rawData1/97/
    pool = Pool(processes=3)
    for case in file_list:
        pat_id = os.path.basename(case)
        pool.apply_async(process, (pat_id,))
    pool.close()
    pool.join()

    with open(log_path, 'a') as f:
        f.write('process stage Done!\n')


def test():
    featureCounts()


if __name__ == "__main__":
    main()

# 处理转移组的脚本
# TODO: try https://github.com/pachterlab/kallisto to quantify expression
