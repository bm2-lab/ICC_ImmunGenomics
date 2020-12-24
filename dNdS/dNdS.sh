for sample_name in `ls /home_2/zhouchi/HCC_only_intrahepatic/mutation_result`
do
	python dNdS_process_input.py ${sample_name}
	vep -i /home/zhouchi/HCC/dNdS/Relapse/${sample_name}/${sample_name}_pass.vcf --cache --dir /home/zhouchi/vep_data_hg19 --dir_cache /home/zhouchi/vep_data_hg19 --force_overwrite --canonical --symbol -o STDOUT --offline -o /home/zhouchi/HCC/dNdS/Relapse/${sample_name}/${sample_name}_vep_ann_all.txt
	python dNdS_input.py ${sample_name}
	Rscript dNdS.R /home/zhouchi/HCC/dNdS/Relapse/${sample_name}/${sample_name}_dNdS_input.txt /home/zhouchi/HCC/dNdS/Relapse/${sample_name}/${sample_name}_dNdS_result.txt
done