###sequenza preprocess
#####test command


tumor_bam_path=$1
normal_bam_path=$2
mutation_vcf=$3
reference_path=$4
gc_file_path=$5
sample_name=$6
out=$7

if [ ! -d ${out} ];then
	mkdir -p ${out}
fi


#sequenza-utils bam2seqz -n ${normal_bam_path} -t ${tumor_bam_path} --fasta ${reference_path} -gc ${gc_file_path} -o ${out}/${sample_name}.out.seqz.gz

#sequenza-utils seqz_binning --seqz ${out}/${sample_name}.out.seqz.gz -w 50 -o ${out}/${sample_name}.small.seqz.gz

#Rscript /home/zhouchi/HCC/script/sequenza_process.R ${out}/${sample_name}.small.seqz.gz ${out} ${sample_name}


python /home/zhouchi/HCC/script/sequenza2pyclone.py ${mutation_vcf} ${out}/${sample_name}_segments.txt ${sample_name} ${out}


TUMOR_CONTENT=`cat ${out}/${sample_name}_cellularity.txt`

PyClone setup_analysis --in_files ${out}/${sample_name}_sequenza2pyclone.txt --tumour_contents ${TUMOR_CONTENT} --prior major_copy_number --working_dir ${out}/${sample_name}_pyclone_result
PyClone run_analysis --config_file ${out}/${sample_name}_pyclone_result/config.yaml 
PyClone build_table --config_file ${out}/${sample_name}_pyclone_result/config.yaml --out_file ${out}/${sample_name}_pyclone_result/loci.tsv --table_type loci
PyClone plot_clusters --config_file ${out}/${sample_name}_pyclone_result/config.yaml --plot_file ${out}/${sample_name}_pyclone_result/coordinates_cluster_plot.pdf --plot_type parallel_coordinates
PyClone plot_clusters --config_file ${out}/${sample_name}_pyclone_result/config.yaml --plot_file ${out}/${sample_name}_pyclone_result/density_cluster_plot.pdf --plot_type density

#if [ -f ${out}/${sample_name}_sequenza2pyclone.txt ];then
#	rm -rf ${out}/${sample_name}.out.seqz.gz ${out}/${sample_name}.small.seqz.gz
#fi