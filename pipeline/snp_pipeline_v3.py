#!/home/mengflz/miniconda3/bin/python

import glob
import os
import hashlib
import subprocess
import opt
import psutil
from sys import argv
import time
from multiprocessing import Pool
from functools import wraps
import traceback


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


def getBigFileMD5(filepath: str):
    md5obj = hashlib.md5()
    maxbuf = 1048576
    # maxbuf = 8192
    f = open(filepath, 'rb')
    while True:
        buf = f.read(maxbuf)
        if not buf:
            break
        md5obj.update(buf)
    f.close()
    hash = md5obj.hexdigest()
    return str(hash).lower()


def checkMD5(filepath: list):
    wrong_md5 = []
    for items in filepath:
        file_list = glob.glob(items+'*.fastq.gz')
        for file_name in file_list:
            with open(file_name + '.md5') as f:
                md5, _ = f.read().split('  ')
                checkmd5 = getBigFileMD5(file_name)
                if md5 != checkmd5:
                    wrong_md5.append(os.path.basename(file_name))
    return wrong_md5


@timer
def map_with_BWA(caseid, settings):
    cmd1 = r"""({SENTIEON_INSTALL_DIR}/bin/sentieon bwa mem -M -R "@RG\tID:rg_{normal_sample}\tSM:{normal_sample}\tPL:{platform}" -t {nt} -K 10000000 {fasta} {normal_fastq_1} {normal_fastq_2} || echo -n 'error' ) | {SENTIEON_INSTALL_DIR}/bin/sentieon util sort -o normal_sorted.bam -t {nt} --sam2bam -i -""".format(**settings)

    cmd2 = r"""({SENTIEON_INSTALL_DIR}/bin/sentieon bwa mem -M -R "@RG\tID:rg_{tumor_sample}\tSM:{tumor_sample}\tPL:{platform}" -t {nt} -K 10000000 {fasta} {tumor_fastq_1} {tumor_fastq_2} || echo -n 'error' ) | {SENTIEON_INSTALL_DIR}/bin/sentieon util sort -o tumor_sorted.bam -t {nt} --sam2bam -i -""".format(**settings)
    with open(log_path, 'a') as f:
        p = subprocess.Popen(cmd1 + '&&' + cmd2, shell=True,
                             stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()


@timer
def sample_matrics(caseid, settings):
    for tissue_type in ['normal', 'tumor']:
        settings['tissue_type'] = tissue_type
        cmd1 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i {tissue_type}_sorted.bam --algo MeanQualityByCycle {tissue_type}_mq_metrics.txt --algo QualDistribution {tissue_type}_qd_metrics.txt --algo GCBias --summary {tissue_type}_gc_summary.txt {tissue_type}_gc_metrics.txt --algo AlignmentStat --adapter_seq '' {tissue_type}_aln_metrics.txt --algo InsertSizeMetricAlgo {tissue_type}_is_metrics.txt".format(
            **settings)
        cmd2 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon plot GCBias -o {tissue_type}_gc-report.pdf {tissue_type}_gc_metrics.txt".format(
            **settings)
        cmd3 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon plot QualDistribution -o {tissue_type}_qd-report.pdf {tissue_type}_qd_metrics.txt".format(
            **settings)
        cmd4 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon plot MeanQualityByCycle -o {tissue_type}_mq-report.pdf {tissue_type}_mq_metrics.txt".format(
            **settings)
        cmd5 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon plot InsertSizeMetricAlgo -o {tissue_type}_is-report.pdf {tissue_type}_is_metrics.txt".format(
            **settings)
        with open(log_path, 'a') as f:
            p = subprocess.Popen(cmd1 + '&&' + cmd2 + '&&' + cmd3 + '&&' + cmd4 +
                                '&&' + cmd5, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
            p.wait()


@timer
def remove_duplicate(caseid, settings):
    for tissue_type in ['normal', 'tumor']:
        settings['tissue_type'] = tissue_type
        cmd1 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -t {nt} -i {tissue_type}_sorted.bam --algo LocusCollector --fun score_info {tissue_type}_score.txt".format(
            **settings)
        cmd2 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -t {nt} -i {tissue_type}_sorted.bam --algo Dedup --score_info {tissue_type}_score.txt --metrics {tissue_type}_dedup_metrics.txt {tissue_type}_deduped.bam".format(
            **settings)
        with open(log_path, 'a') as f:
            p = subprocess.Popen(cmd1 + '&&' + cmd2, shell=True,
                                stderr=subprocess.STDOUT, stdout=f, text=True)
            p.wait()


@timer
def coverage_metrics(caseid, settings):
    for tissue_type in ['normal', 'tumor']:
        settings['tissue_type'] = tissue_type
        cmd = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i {tissue_type}_deduped.bam --algo CoverageMetrics {tissue_type}_coverage_metrics".format(
            **settings)
        with open(log_path, 'a') as f:
            p = subprocess.Popen(
                cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
            p.wait()


@timer
def indel_realigner(caseid, settings):
    for tissue_type in ['normal', 'tumor']:
        settings['tissue_type'] = tissue_type
        cmd = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i {tissue_type}_deduped.bam --algo Realigner -k {known_Mills_indels} -k {known_1000G_indels} {tissue_type}_realigned.bam".format(**settings)
        with open(log_path, 'a') as f:
            p = subprocess.Popen(
                cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
            p.wait()


@timer
def base_recalibration(caseid, settings):
    for tissue_type in ['normal', 'tumor']:
        settings['tissue_type'] = tissue_type
        cmd1 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i {tissue_type}_deduped.bam --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {tissue_type}_recal_data.table".format(
            **settings)
        cmd2 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i {tissue_type}_deduped.bam -q {tissue_type}_recal_data.table --algo QualCal -k {dbsnp} -k {known_Mills_indels} -k {known_1000G_indels} {tissue_type}_recal_data.table.post".format(
            **settings)
        cmd3 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -t {nt} --algo QualCal --plot --before {tissue_type}_recal_data.table --after {tissue_type}_recal_data.table.post {tissue_type}_recal.csv".format(
            **settings)
        cmd4 = r"{SENTIEON_INSTALL_DIR}/bin/sentieon plot QualCal -o {tissue_type}_recal_plots.pdf {tissue_type}_recal.csv".format(
            **settings)
        with open(log_path, 'a') as f:
            p = subprocess.Popen(cmd1 + '&&' + cmd2 + '&&' + cmd3 + '&&' + cmd4, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
            p.wait()


@timer
def corealignment(caseid, settings):
    cmd = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i tumor_realigned.bam -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table --algo Realigner -k {known_Mills_indels} -k {known_1000G_indels} tn_corealigned.bam".format(**settings)
    with open(log_path, 'a') as f:
        p = subprocess.Popen(
            cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()


@timer
def variant_calling(caseid, settings):
    # cmd_ori = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} -i tn_corealigned.bam --algo TNhaplotyper --tumor_sample {tumor_sample} --normal_sample {normal_sample} --dbsnp {dbsnp} output-tnhaplotyper.vcf.gz".format(
    #     **settings)
    cmd = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} --interval {bed} -i tumor_realigned.bam  -i normal_realigned.bam  -q tumor_recal_data.table  -q normal_recal_data.table --algo TNhaplotyper2 --tumor_sample {tumor_sample} --normal_sample {normal_sample}  --germline_vcf {germline} --pon {pon} tmp_output_filter.vcf --algo OrientationBias --tumor_sample {tumor_sample} orientation_data_file".format(**settings)
    with open(log_path, 'a') as f:
        p = subprocess.Popen(
            cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()


@timer
def filtration(caseid, settings):
    """use bcftools/TNfilter to filter output-tnhaplotyper.vcf.gz"""
    # cmd_bcftools = r"bcftools filter -m + -i '(FMT/AF[1:0] >= 0.01 && SUM(FMT/AD[1:])>=30 && FMT/AD[1:1] >=3 && FMT/AD[0:1] <=3)' output-tnhaplotyper.vcf.gz -o output_filtered.vcf -O v"
    cmd = r"{SENTIEON_INSTALL_DIR}/bin/sentieon driver -r {fasta} -t {nt} --algo TNfilter --tumor_sample {tumor_sample} --normal_sample {normal_sample} -v  tmp_output_filter.vcf --max_fp_rate {max_fp_rate} --orientation_priors orientation_data_file output_filter.vcf".format(**settings)
    with open(log_path, 'a') as f:
        p = subprocess.Popen(
            cmd, shell=True, stderr=subprocess.STDOUT, stdout=f, text=True)
        p.wait()


def check_license(license_path):
    flag = 0
    for proc in psutil.process_iter():
        if 'licsrvr' in proc.name():
            flag = 1
    if flag == 0:
        print('start licsrvr!')
        subprocess.run(
            r'/home_2/mengflz/sentieon-genomics-201911/libexec/licsrvr --start {}'.format(license_path), shell=True)


@timer
def trash_remove(caseid, settings):
    # need to return a list recording if got wrong
    # remain_file = ['output_filtered.vcf', 'output-tnhaplotyper.vcf.gz', 'output-tnhaplotyper.vcf','output-tnhaplotyper.vcf.gz.tbi',
    # 'normal_realigned.bam', 'tumor_realigned.bam', 'normal_realigned.bam.bai', 'tumor_realigned.bam.bai']
    remain_file = ['tmp_output_filter.vcf', 'normal_realigned.bam', 'tumor_realigned.bam', 'normal_realigned.bam.bai', 'tumor_realigned.bam.bai', 'output_filter.vcf']
    file_total = glob.glob('*')
    with open(log_path, 'a') as f:
        for i in file_total:
            if i not in remain_file:
                os.remove(i)
                f.write('{} remove file {}\n'.format(caseid, i))
        for j in remain_file:
            if j not in file_total:
                f.write('{} WARNINGS: no needed file {}\n'.format(caseid, j))


def TNseq(caseid, settings: dict):
    '''
    rewrited by pipeline-example-TNscope.sh
    wrong_md5 storage file names not passing md5 check
    '''

    os.environ['SENTIEON_LICENSE'] = settings['SENTIEON_LICENSE']
    check_license(settings['SENTIEON_LICENSE'])

    # Step1 Mapping reads with BWA-MEM, sorting for normal and tumor sample
    map_with_BWA(caseid, settings)
    # Step2 Metrics for normal and tumor sample
    sample_matrics(caseid, settings)
    # Step3 Remove Duplicate Reads for normal and tumor sample
    remove_duplicate(caseid, settings)
    # Step4 Coverage metrics for normal and tumor sample
    coverage_metrics(caseid, settings)
    # Step5 Indel realigner
    indel_realigner(caseid, settings)
    # Step6 Base recalibration for normal and tumor sample
    base_recalibration(caseid, settings)
    # Step7 Corealignment of tumor and normal
    # corealignment(caseid, settings)
    # Step8 Variant calling
    variant_calling(caseid, settings)
    # Step9 Filteration
    filtration(caseid, settings)
    # Step final trash remove
    trash_remove(caseid, settings)
    

def process(case_id):
    '''
    total process func
    '''
    # mkdir patid
    if not os.path.exists(os.getcwd()+'/'+case_id):
        os.mkdir(os.getcwd()+'/'+case_id)
    os.chdir(os.getcwd()+'/'+case_id)

    try:
        settings = opt.snp_options(case_id, nt=30).__dict__
        TNseq(case_id, settings)
    except Exception:
        with open(log_path, 'a') as f:
            f.write(traceback.format_exc())
    finally:
        os.chdir(os.path.dirname(os.getcwd()))


def main():
    '''
    input: sample_path; data_path
        e.g. data_path = "/home_2/wzt/outcome/" process_num=3

    output: ./output/patient_id/[tumor_sample]+[normal_sample]/xxx
            ./output/pipeline.log
    '''
    global data_path
    data_path = argv[1]

    # make output dir and log file
    output_path = os.getcwd() + '/output_v3'
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    os.chdir(output_path)
    global log_path
    log_path = os.getcwd()+'/pipeline.log'

    with open(log_path, 'a') as f:
        f.write("*********" + time.strftime("%a %b %d %H:%M:%S %Y",
                                            time.localtime()) + "**********" + '\n')
    try:
        pool = Pool(processes=2)
        file_list = glob.glob(data_path + '*')
        for case in file_list:
            pat_id = os.path.basename(case)
            if os.path.exists(os.path.normpath(case + '/WES')):
                pool.apply_async(process, (pat_id,))
        pool.close()
        pool.join()
    except Exception as e:
        traceback.print_exc()


def test():
    settings = opt.snp_options(nt=60, normal_folder="/home_2/ICC_outcome/WES/Sample_RB20002146-T200387ND1L1-200813/", tumor_folder="/home_2/ICC_outcome/WES/Sample_RB20002147-T200386TD1L1-200813/").__dict__
    os.environ['SENTIEON_LICENSE'] = settings['SENTIEON_LICENSE']
    check_license(settings['SENTIEON_LICENSE'])
    global log_path
    log_path = './0108.log'
    corealignment("164", settings)
    variant_calling("164", settings)


if __name__ == '__main__':
    main()
    # test()


# 相对于v2，修改了读入文件的方法，直接从/home_2/wzt/outcome/路径中返回列表并读入文件
# 修改了opt.py
