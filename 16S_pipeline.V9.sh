#/bin/sh -S
#########
#Please address any bugs to Cheng.
#Date 2018.9.13
#########

#echo $(readlink -f $1)

sample_metadata=$(readlink -f $1)
depth=$2

taxa_filtered="${3}"

category_sum=${4//,/ }
reference_trained=$(readlink -f $5)
close_reference_trained=$(readlink -f $6)
manifest_file=$(readlink -f $7)
not_rda=${8//\;/ }
classifier_type=$9
if_picrust=${10}


echo "Check wheather your categories are the following:"
for i in $category_sum;do echo $i;done

echo "Check wheather the group of enviromental factors excluded from rda are the following:"
for i in $not_rda;do echo $i;done


if [[ $if_picrust == 'yes' ]];then
	echo "

	Will run PICRUST"
else echo "

	Won't run PICRUST"
fi

declare -A tax_aa;
tax_aa=([k]=Kingdom [p]=Phylum [c]=Class [o]=Order [f]=Family [g]=Genus [s]=Species);

tax_levels["1"]="Kingdom"
tax_levels["2"]="Phylum"
tax_levels["3"]="Class"
tax_levels["4"]="Order"
tax_levels["5"]="Family"
tax_levels["6"]="Genus"
tax_levels["7"]="Species"

#tax_levels1["k"]="Kingdom"
#tax_levels1["p"]="Phylum"
#tax_levels1["c"]="Class"
#tax_levels1["o"]="Order"
#tax_levels1["f"]="Family"
#tax_levels1["g"]="Genus"
#tax_levels1["s"]="Species"


if [ -z "$9" ]; then
	echo "##########

		  Please prepare the directory structure before starting the program like below:
		  raw/fastq_files ...
		  mapping_file
		  manifest_file
		  \n\n"

	echo "Please provide following input parameters
		1) Path of the mapping file. (Accept both .csv or txt format.)
		2) Depth of the for subsampleing. (Suggested value: 4000)
		3) Taxa excluded from OTU table, seprated by commas.
		4) The name of categories in the mapping file seprated by commas.
		5) Path of the classifier for alignment.
		6) Path of the reference sequences for close reference alignment.
		7) Path of the manifest file.
		8) Specify numeric variables excluded from rda seprated by commas,use 'none' if all numeric variables is expected. use 'all' to skip correlation analysis.
		9) Specify the type of classifier, either silva or gg
		10) wheather to run PICRUST, 'yes' or 'no'
		Sample Usage:
		bash ~/github/Bayegy/16S_pipeline.V9.sh ../data/sample-metadata.tsv auto 1000 Group ~/database_16S/GG/338-806/gg_13_8_99_338_806_classifier.qza ~/database_16S/GG/338-806/gg_13_5_97_338_806_ref_seqs.qza ../data/manifest.txt  all gg yes
		"
	exit 0
else
	echo "################
	Running: sh $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} "
fi

check_file() {
	echo "Checking file for $1 ..."
	file_name=$1
	if [ -f $file_name ]; then
		echo "File $file_name exists"
	else
		echo "File $file_name does not exist"
		exit
	fi
}


	##Activate Qiime2 Version
	echo "##############################################################\n#Initiate directory name and set up the directory structure"

	SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

	source $SCRIPTPATH/path/clean_path.sh

	python3 $SCRIPTPATH/piputils/write_colors_plan.py -i $sample_metadata -c $4 -p $SCRIPTPATH/piputils/group_color.list -o colors_plan.json
	if [ -f colors_plan.json ];then
		source $SCRIPTPATH/piputils/put_colors_plan_path.sh $(readlink -f colors_plan.json)
		echo "COLORS_PLAN_PATH = ${COLORS_PLAN_PATH}"
	fi
	source activate qiime2-2018.11

	export TMPDIR='/home/admin1/TMPDIR'
	echo "Classifier TMPDIR is ${TMPDIR}"
	export JOBLIB_TEMP_FOLDER='/home/admin1/JOBLIB_TEMP_FOLDER'
	echo "Classifier JOBLIB_TEMP_FOLDER is ${JOBLIB_TEMP_FOLDER}"


	echo "##############################################################\n#paired end analysis using DADA2"

	qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path $manifest_file --output-path demux.qza --input-format PairedEndFastqManifestPhred33
	qiime demux summarize --i-data demux.qza --o-visualization demux.qzv


	#qiime dada2 denoise-paired  --i-demultiplexed-seqs demux.qza --p-trunc-len-f 290 --p-trunc-len-r 256 --p-trim-left-f 26 --p-trim-left-r 26 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza  --p-n-threads 0 --o-denoising-stats stats-dada2.qza --verbose
	#qiime dada2 denoise-paired --p-max-ee 100 --p-trunc-q 0  --i-demultiplexed-seqs demux.qza --p-trunc-len-f 290 --p-trunc-len-r 220 --p-trim-left-f 26 --p-trim-left-r 26 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza  --p-n-threads 0 --o-denoising-stats stats-dada2.qza --verbose
	qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-trunc-len-f 290 --p-trunc-len-r 256 --p-trim-left-f 26 --p-trim-left-r 26 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza  --p-n-threads 0 --o-denoising-stats stats-dada2.qza --verbose
	#qiime dada2 denoise-paired  --i-demultiplexed-seqs demux.qza --p-trunc-len-f 290 --p-trunc-len-r 220 --p-trim-left-f 26 --p-trim-left-r 26 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza  --p-n-threads 0 --o-denoising-stats stats-dada2.qza --verbose
	qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv
<<comment2
	echo "##############################################################\n#Single end analysis using DADA2"
	qiime tools import   --type 'SampleData[SequencesWithQuality]'   --input-path $manifest_file --output-path demux.qza --input-format SingleEndFastqManifestPhred33
	qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
	qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-max-ee 50 --p-trunc-q 0 --p-trunc-len 240 --p-trim-left 0 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza  --p-n-threads 0 --o-denoising-stats stats-dada2.qza --verbose
	qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv
comment2



<<comment1
	echo "############################################################\nAlternative methods of read-joining using deblur"
	qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path $manifest_file --output-path demux.qza --input-format PairedEndFastqManifestPhred33
	#qiime cutadapt trim-paired --i-demultiplexed-sequences demux.qza --p-front-f CCAGCASCYGCGGTAATTCC --p-front-r ACTTTCGTTCTTGATYRA --p-overlap 10 --o-trimmed-sequences trimmed-demux.qza
	#qiime tools import   --type 'SampleData[JoinedSequencesWithQuality]'  --input-path $manifest_file --output-path demux-joined.qza --input-format SingleEndFastqManifestPhred33
	qiime vsearch join-pairs --p-maxdiffs 5 --p-minovlen 15 --p-truncqual 2 --i-demultiplexed-seqs demux.qza --o-joined-sequences demux-joined.qza

	qiime demux summarize --i-data demux-joined.qza --o-visualization demux.qzv
	qiime quality-filter q-score-joined --i-demux demux-joined.qza --o-filtered-sequences demux-joined-filtered.qza --o-filter-stats demux-joined-filter-stats.qza
	qiime deblur denoise-16S --i-demultiplexed-seqs demux-joined-filtered.qza --p-trim-length 400  --p-sample-stats --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-stats stats-dada2.qza --p-jobs-to-start 16 --verbose
	qiime deblur visualize-stats --i-deblur-stats stats-dada2.qza --o-visualization stats-dada2.qzv
comment1



	mv rep-seqs-dada2.qza rep-seqs.withCandM.qza
	mv table-dada2.qza table.withCandM.qza

	# Output the OTUid_frequency-talbe txt file, this file is the raw talbe, without any filtering.
	echo "qiime tools export --input-path table.withCandM.qza --output-path table.withCandM.without_filtering"
	qiime tools export --input-path table.withCandM.qza --output-path table.withCandM.without_filtering
	echo "biom convert --to-tsv -i table.withCandM.without_filtering/feature-table.biom -o table.withCandM.without_filtering/feature-table.tsv"
	biom convert --to-tsv -i table.withCandM.without_filtering/feature-table.biom -o table.withCandM.without_filtering/feature-table.tsv

	echo "##############################################################\n#Classify the taxonomy"
	check_file $reference_trained
	qiime feature-classifier classify-sklearn --verbose --p-confidence 0.7 --p-n-jobs 1   --i-classifier $reference_trained  --i-reads rep-seqs.withCandM.qza  --o-classification taxonomy.withCandM.qza
	# qiime feature-classifier classify-sklearn --verbose --p-confidence 0.55 --p-n-jobs 1   --i-classifier $reference_trained  --i-query rep-seqs.withCandM.qza  --o-classification taxonomy.withCandM.qza
	#qiime feature-classifier  classify-consensus-vsearch  --i-query rep-seqs.withCandM.qza --i-reference-reads ~/database_ITS/UNITE_release/full-length/rep-set.qza  --i-reference-taxonomy ~/database_ITS/UNITE_release/full-length/ref-taxonomy.qza  --p-perc-identity 0.7 --o-classification taxonomy.withCandM.qza

	qiime metadata tabulate  --m-input-file taxonomy.withCandM.qza  --o-visualization taxonomy.withCandM.qzv

	# Output the OTUid_taxonomy txt file, without filtering the mitochondria,chloroplast,Unassigned
	echo "qiime tools export --input-path taxonomy.withCandM.qza --output-path taxonomy.withCandM.without_filtering"
	qiime tools export --input-path taxonomy.withCandM.qza --output-path taxonomy.withCandM.without_filtering

	#Archaea,
	# qiime taxa filter-table   --i-table table.withCandM.qza  --i-taxonomy taxonomy.withCandM.qza  --p-exclude ${taxa_filtered}mitochondria,chloroplast,Unassigned  --o-filtered-table table-no-mitochondria-no-chloroplast.qza
	# mv table-no-mitochondria-no-chloroplast.qza table.qza
	# qiime taxa filter-seqs   --i-sequences rep-seqs.withCandM.qza   --i-taxonomy taxonomy.withCandM.qza  --p-exclude ${taxa_filtered}mitochondria,chloroplast,Unassigned   --o-filtered-sequences rep-seqs-no-mitochondria-no-chloroplast.qza
	# mv rep-seqs-no-mitochondria-no-chloroplast.qza rep-seqs.qza
	# qiime feature-classifier classify-sklearn --verbose --p-confidence 0.4 --p-n-jobs 1   --i-classifier $reference_trained  --i-reads rep-seqs.qza  --o-classification taxonomy.qza
	echo "Filter parameter is ${taxa_filtered}"
	filter_source.py -i table.withCandM.qza -r rep-seqs.withCandM.qza -t taxonomy.withCandM.qza -f "${taxa_filtered}"
	if [[ $classifier_type == 'silva' ]];
		then python $SCRIPTPATH/format_silva_to_gg.py -i taxonomy.qza;
		else python $SCRIPTPATH/format_silva_to_gg.py -i taxonomy.qza -c;
	fi;


	echo "##############################################################\n#Generate tree";
	qiime alignment mafft   --i-sequences rep-seqs.qza   --o-alignment aligned-rep-seqs.qza
	qiime alignment mask   --i-alignment aligned-rep-seqs.qza   --o-masked-alignment masked-aligned-rep-seqs.qza
	qiime phylogeny fasttree   --i-alignment masked-aligned-rep-seqs.qza   --o-tree unrooted-tree.qza
	qiime phylogeny midpoint-root   --i-tree unrooted-tree.qza   --o-rooted-tree rooted-tree.qza
	qiime feature-table tabulate-seqs   --i-data rep-seqs.qza   --o-visualization rep-seqs.qzv


	echo "##############################################################\n#Generate OTUStats";
	for f in rep-seqs.qza table.qza taxonomy.qza ; do echo $f; qiime tools export --input-path $f --output-path exported; done
	#for f in alpha-rarefaction.qzv table.qzv taxa-bar-plots.qzv; do echo $f; qiime tools export --input-path $f --output-path exported_qzv; done
	qiime tools export --input-path rooted-tree.qza --output-path exported/
	mv exported/tree.nwk exported/tree.rooted.nwk 
	qiime tools export --input-path unrooted-tree.qza --output-path exported/
	mv exported/tree.nwk exported/tree.unrooted.nwk 
	biom add-metadata -i exported/feature-table.biom -o exported/feature-table.taxonomy.biom --observation-metadata-fp exported/taxonomy.tsv --observation-header OTUID,taxonomy,confidence
	biom convert -i exported/feature-table.taxonomy.biom -o exported/feature-table.taxonomy.txt --to-tsv --header-key taxonomy
	for f in $(find . -type f -name "*.qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); new=${dir}/${base}; qiime tools export --input-path $f --output-path ${new}.qzv.exported; done 

	source deactivate
	source activate qm2
	echo "##############################################################\n#Generate the figure for the percentage of annotated level"
	perl ${SCRIPTPATH}/stat_otu_tab.unspecifiedadded.pl -unif min exported/feature-table.taxonomy.txt -prefix exported/Relative/otu_table --even exported/Relative/otu_table.even.txt -spestat exported/Relative/classified_stat_relative.xls > /dev/null
	perl ${SCRIPTPATH}/bar_diagram.pl -table exported/Relative/classified_stat_relative.xls -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,7 > exported/Relative/Classified_stat_relative.svg
	for svg_file in exported/Relative/*svg; do echo $svg_file; n=$(basename "$svg_file" .svg); echo $n; rsvg-convert -h 3200 -b white $svg_file > exported/Relative/${n}.png; done

	if [ -d "./Result" ];then
		rm -r ./Result;
	fi;
	mkdir -p ./Result/1-OTUStats/3-RepresentiveSequence/;

	cp -r demux.qzv* stats-dada2.qzv* ./Result/1-OTUStats/
	cp -r ./rep-seqs.qzv* ./exported/*nwk ./exported/dna-sequences.fasta ./Result/1-OTUStats/3-RepresentiveSequence/

	for f in $(find ./Result/1-OTUStats/ -type f -name "*qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); mv $f ${f}.exported; mv ${f}.exported ${dir}/${base}; done
	for f in $(find ./Result/1-OTUStats/ -type f -name "index.html") ; do echo $f; base=$(basename $f .html); dir=$(dirname $f); new=${dir}/Summary_请点此文件查看.html; mv $f $new; done


	mv ./Result/1-OTUStats/stats-dada2 ./Result/1-OTUStats/2-Stats-dada2
	mv ./Result/1-OTUStats/demux/ ./Result/1-OTUStats/1-Stats-demux
	cp -r exported/feature-table.taxonomy.txt exported/feature-table.taxonomy.biom exported/Relative/Classified_stat_relative.png exported/Relative/classified_stat_relative.xls ./Result/1-OTUStats/
	cp -r exported/Relative/otu_table.even.txt ./Result/1-OTUStats/feature-table.taxonomy.even.txt


	echo "##############################################################\n#Generate the results of each group"
	for category_set in $category_sum;
		do echo $category_set;
		source deactivate;
		source activate qiime2-2018.11;
		#mkdir ${category_set}_Results;
		python $SCRIPTPATH/filter_source.py  -i table.qza -t taxonomy.qza -r rep-seqs.qza -m $sample_metadata -c $category_set -o ${category_set}_Results;
		cd ${category_set}_Results;
		mapping_file=$(readlink -f './mapping_file.txt');


		qiime metadata tabulate   --m-input-file taxonomy.qza   --o-visualization taxonomy.qzv;

		echo "##############################################################\n#Generate tree";
		qiime alignment mafft   --i-sequences rep-seqs.qza   --o-alignment aligned-rep-seqs.qza
		qiime alignment mask   --i-alignment aligned-rep-seqs.qza   --o-masked-alignment masked-aligned-rep-seqs.qza
		qiime phylogeny fasttree   --i-alignment masked-aligned-rep-seqs.qza   --o-tree unrooted-tree.qza
		qiime phylogeny midpoint-root   --i-tree unrooted-tree.qza   --o-rooted-tree rooted-tree.qza


		echo "##############################################################\n#Visulize of the table without Choloroplast and Mitochondira"
		qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file $mapping_file
		qiime taxa barplot   --i-table table.qza   --i-taxonomy taxonomy.qza   --m-metadata-file $mapping_file  --o-visualization taxa-bar-plots.qzv


		#########calculate the min sample depth
		#for f in rep-seqs.qza table.qza taxonomy.qza ; do echo $f; qiime tools export --input-path $f --output-path exported; done
		qiime tools export --input-path table.qzv --output-path exported_qzv
		if [[ $depth == 'auto' ]];
			then min_depth=$(echo \($(cut -f2 -d ',' exported_qzv/sample-frequency-detail.csv | sort -n | head -n1)/1000\)*1000 | bc);
			else min_depth=$depth;
		fi;



		echo "##############################################################The selected sample depth is $min_depth"

		echo "##############################################################\n#Core alpha and beta diversity analysis"
		qiime diversity core-metrics-phylogenetic   --i-phylogeny rooted-tree.qza   --i-table table.qza   --p-sampling-depth $min_depth   --m-metadata-file $mapping_file  --output-dir core-metrics-results

		#qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/faith_pd_vector.qza   --m-metadata-file $mapping_file  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
		#qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/evenness_vector.qza   --m-metadata-file $mapping_file  --o-visualization core-metrics-results/evenness-group-significance.qzv
		#qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/shannon_vector.qza   --m-metadata-file $mapping_file  --o-visualization core-metrics-results/shannon-group-significance.qzv
		#qiime diversity alpha-group-significance   --i-alpha-diversity core-metrics-results/observed_otus_vector.qza   --m-metadata-file $mapping_file  --o-visualization core-metrics-results/observed_otus-group-significance.qzv
		for category_1 in $category_set;
		do echo $category_1;
			qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method permanova --m-metadata-column $category_1   --o-visualization 'core-metrics-results/unweighted_unifrac-permanova-'$category_1'-significance.qzv'  --p-pairwise;
			#qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method permanova --m-metadata-column $category_2   --o-visualization 'core-metrics-results/unweighted_unifrac-permanova-'$category_2'-significance.qzv'  --p-pairwise
			qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method permanova --m-metadata-column $category_1   --o-visualization 'core-metrics-results/weighted_unifrac-permanova-'$category_1'-significance.qzv'  --p-pairwise;
			#qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method permanova --m-metadata-column $category_2   --o-visualization 'core-metrics-results/weighted_unifrac-permanova-'$category_2'-significance.qzv'  --p-pairwise
			qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method permanova --m-metadata-column $category_1   --o-visualization 'core-metrics-results/bray_curtis-permanova-'$category_1'-significance.qzv'  --p-pairwise;
			#qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method permanova --m-metadata-column $category_2   --o-visualization 'core-metrics-results/bray_curtis-permanova-'$category_2'-significance.qzv'  --p-pairwise
			qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method anosim --m-metadata-column $category_1   --o-visualization 'core-metrics-results/unweighted_unifrac-anosim-'$category_1'-significance.qzv'  --p-pairwise;
			#qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method anosim --m-metadata-column $category_2   --o-visualization 'core-metrics-results/unweighted_unifrac-anosim-'$category_2'-significance.qzv'  --p-pairwise
			qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method anosim --m-metadata-column $category_1   --o-visualization 'core-metrics-results/weighted_unifrac-anosim-'$category_1'-significance.qzv'  --p-pairwise;
			#qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method anosim --m-metadata-column $category_2   --o-visualization 'core-metrics-results/weighted_unifrac-anosim-'$category_2'-significance.qzv'  --p-pairwise
			qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method anosim --m-metadata-column $category_1   --o-visualization 'core-metrics-results/bray_curtis-anosim-'$category_1'-significance.qzv'  --p-pairwise;
		done;

		#qiime diversity beta-group-significance   --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza   --m-metadata-file $mapping_file  --p-method anosim --m-metadata-column $category_2   --o-visualization 'core-metrics-results/bray_curtis-anosim-'$category_2'-significance.qzv'  --p-pairwise
		qiime diversity alpha-rarefaction   --i-table table.qza   --i-phylogeny rooted-tree.qza   --p-max-depth $min_depth  --m-metadata-file $mapping_file  --o-visualization alpha-rarefaction.qzv   --p-steps 50
		##These following two commands work only for column with numeric values:
		##qiime emperor plot   --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza   --m-metadata-file $mapping_file --p-custom-axis $category_2   --o-visualization 'core-metrics-results/unweighted_unifrac-emperor-'$category_2'.qzv'
		##qiime emperor plot   --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza   --m-metadata-file $mapping_file   --p-custom-axis $category_2   --o-visualization 'core-metrics-results/bray-curtis-emperor-'$category_2'.qzv'


		echo "##############################################################\n#alpha dviersity summary"
		mkdir alpha alpha_groupsignificance
		#qiime diversity alpha --i-table table.qza --p-metric chao1 --output-dir alpha/chao1
		#qiime diversity alpha --i-table table.qza --p-metric shannon --output-dir alpha/shannon
		#qiime diversity alpha --i-table table.qza --p-metric observed_otus --output-dir alpha/observed_otus
		qiime diversity alpha-phylogenetic --i-table table.qza --i-phylogeny rooted-tree.qza --p-metric faith_pd --output-dir alpha/faith_pd
		#qiime diversity alpha --i-table table.qza --p-metric simpson --o-alpha-diversity alpha/simpson.qza
		#qiime diversity alpha --i-table table.qza --p-metric ace --o-alpha-diversity alpha/ace.qza
		#qiime diversity alpha --i-table table.qza --p-metric goods_coverage --o-alpha-diversity alpha/goods_coverage.qza

		for alpha_index in chao1 shannon observed_otus faith_pd simpson;
			do echo $alpha_index;
			if [ ! $alpha_index == faith_pd ];then
			qiime diversity alpha --i-table table.qza --p-metric ${alpha_index} --output-dir alpha/${alpha_index}/
			fi;
			qiime tools export --input-path alpha/${alpha_index}/alpha_diversity.qza --output-path alpha/${alpha_index}/;
			#if [ ! $alpha_index == goods_coverage ];then
			qiime diversity alpha-group-significance   --i-alpha-diversity alpha/${alpha_index}/alpha_diversity.qza   --m-metadata-file $mapping_file  --o-visualization alpha_groupsignificance/${alpha_index}-group-significance.qzv;
			#fi;
		done;



	 	#qiime tools export --input-path alpha/chao1/alpha_diversity.qza --output-path alpha/chao1/
	 	#qiime tools export --input-path alpha/shannon/alpha_diversity.qza --output-path alpha/shannon/
	 	#qiime tools export --input-path alpha/observed_otus/alpha_diversity.qza --output-path alpha/observed_otus/
	 	#qiime tools export --input-path alpha/faith_pd/alpha_diversity.qza --output-path alpha/faith_pd/
	 	#qiime tools export --input-path alpha/simpson.qza --output-path alpha/simpson/
	 	#qiime tools export --input-path alpha/ace.qza --output-path alpha/ace/
	 	#qiime tools export --input-path alpha/goods_coverage.qza --output-path alpha/goods_coverage/
	 	#paste alpha/observed_otus/alpha-diversity.tsv alpha/chao1/alpha-diversity.tsv alpha/shannon/alpha-diversity.tsv alpha/faith_pd/alpha-diversity.tsv alpha/simpson/alpha-diversity.tsv alpha/ace/alpha-diversity.tsv alpha/goods_coverage/alpha-diversity.tsv | awk -F'\t' 'BEGIN{OFS="\t"}{print $1, $2, $4, $6, $8, $10, $12, $14}' >  alpha/alpha-summary.tsv
	 	paste alpha/observed_otus/alpha-diversity.tsv alpha/chao1/alpha-diversity.tsv alpha/shannon/alpha-diversity.tsv alpha/faith_pd/alpha-diversity.tsv alpha/simpson/alpha-diversity.tsv | awk -F'\t' 'BEGIN{OFS="\t"}{print $1, $2, $4, $6, $8, $10}' >  alpha/alpha-summary.tsv

		echo "##############################################################\n#Export necessary files for future analysis"
		for f in rep-seqs.qza table.qza taxonomy.qza ; do echo $f; qiime tools export --input-path $f --output-path exported; done
		#for f in alpha-rarefaction.qzv table.qzv taxa-bar-plots.qzv; do echo $f; qiime tools export --input-path $f --output-path exported_qzv; done
		qiime tools export --input-path rooted-tree.qza --output-path exported/
		mv exported/tree.nwk exported/tree.rooted.nwk 
		qiime tools export --input-path unrooted-tree.qza --output-path exported/
		mv exported/tree.nwk exported/tree.unrooted.nwk 
		biom add-metadata -i exported/feature-table.biom -o exported/feature-table.taxonomy.biom --observation-metadata-fp exported/taxonomy.tsv --observation-header OTUID,taxonomy,confidence
		biom convert -i exported/feature-table.taxonomy.biom -o exported/feature-table.taxonomy.txt --to-tsv --header-key taxonomy
		biom convert -i exported/feature-table.taxonomy.biom -o exported/feature-table.txt --to-tsv

		echo "##############################################################\n#Generate heatmaps for top OTUs with different levels with minimum frequence reads supported"
		mkdir exported/collapsed
		# mkdir exported/${min_freq}
		for n in 2 3 4 5 6 7;
			do echo $n;
			qiime taxa collapse   --i-table table.qza   --i-taxonomy taxonomy.qza   --p-level $n   --o-collapsed-table exported/collapsed/collapsed-${tax_levels[${n}]}.qza;
			qiime feature-table summarize --i-table exported/collapsed/collapsed-${tax_levels[${n}]}.qza --o-visualization exported/collapsed/collapsed-${tax_levels[${n}]}.qzv '--m-sample-metadata-file' $mapping_file;
			# qiime feature-table filter-features   --i-table exported/collapsed/collapsed-${tax_levels[${n}]}.qza --p-min-frequency $min_freq  --o-filtered-table exported/${min_freq}/table-${tax_levels[${n}]}.${min_freq}.qza; 
<<COMMENT
			for category_1 in $category_set;
				do echo $category_1;
					Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -t exported/${min_freq}/table-${tax_levels[${n}]}.${min_freq}.qza -o media_files
					qiime feature-table heatmap --i-table media_files/filtered_feature_table.qza  --m-metadata-file media_files/cleaned_map.txt --m-metadata-column $category_1 --o-visualization exported/${min_freq}/${category_1}-table-${tax_levels[${n}]}.${min_freq}.qzv;
				done;
COMMENT
		done;

		source deactivate
		source activate qm2
		echo "##############################################################\n#Generate the figure for the percentage of annotated level"
		perl ${SCRIPTPATH}/stat_otu_tab.unspecifiedadded.pl -unif min exported/feature-table.taxonomy.txt -prefix exported/Relative/otu_table --even exported/Relative/otu_table.even.txt -spestat exported/Relative/classified_stat_relative.xls
		perl ${SCRIPTPATH}/bar_diagram.pl -table exported/Relative/classified_stat_relative.xls -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,7 > exported/Relative/Classified_stat_relative.svg
		# Rscript ${SCRIPTPATH}/taxonomy_stats.R -i exported/feature-table.ConsensusLineage.txt -o TaxaStatistic
		sed 's/taxonomy/Consensus Lineage/' < exported/feature-table.taxonomy.txt | sed 's/ConsensusLineage/Consensus Lineage/' > exported/feature-table.ConsensusLineage-raw.txt
		sed 's/taxonomy/Consensus Lineage/' < exported/Relative/otu_table.even.txt | sed 's/ConsensusLineage/Consensus Lineage/' > exported/feature-table.ConsensusLineage.txt
		sed -i '1i\# Constructed from biom file'  exported/feature-table.ConsensusLineage.txt 


		Rscript ${SCRIPTPATH}/taxonomy_stats.R -i exported/feature-table.ConsensusLineage-raw.txt -o TaxaStatistic

		Rscript ${SCRIPTPATH}/collapse_table_with_group_mean.R -m $mapping_file -c $category_set -t exported/Relative/classified_stat_relative.xls -o ./
		perl ${SCRIPTPATH}/bar_diagram.pl -table ${category_set}_classified_stat_relative.xls -style 1 -x_title "" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,7 > ${category_set}_classified_stat_relative.svg
		for svg_file in *svg; do echo $svg_file; n=$(basename "$svg_file" .svg); echo $n; rsvg-convert -h 3200 -b white $svg_file > ${n}.png; done

		for key in ${!tax_aa[*]};do mv exported/Relative/otu_table.${key}.relative.mat exported/Relative/otu_table.${tax_aa[$key]}.relative.txt;done;
			#mv exported/Relative/otu_table.p.relative.mat exported/Relative/otu_table.Phylum.relative.txt
			#mv exported/Relative/otu_table.c.relative.mat exported/Relative/otu_table.Class.relative.txt
			#mv exported/Relative/otu_table.o.relative.mat exported/Relative/otu_table.Order.relative.txt
			#mv exported/Relative/otu_table.f.relative.mat exported/Relative/otu_table.Family.relative.txt
			#mv exported/Relative/otu_table.g.relative.mat exported/Relative/otu_table.Genus.relative.txt
			#mv exported/Relative/otu_table.s.relative.mat exported/Relative/otu_table.Species.relative.txt


#		for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species"; 
#			do echo $n7; 
			#echo "mv exported/Relative/otu_table.${n7}.relative.mat exported/Relative/otu_table.${tax_levels[${n7}]}.relative.txt"
			#mv exported/Relative/otu_table.${n7}.relative.mat exported/Relative/otu_table.${tax_levels[${n7}]}.relative.txt
			#echo ${tax_levels[$n7]}
			#echo ${tax_levels[${n7}]}
#			perl -lane '$,="\t";pop(@F);print(@F)' exported/Relative/otu_table.${n7}.relative.txt > exported/Relative/otu_table.${n7}.relative.lastcolumn.txt; 
#			perl ${SCRIPTPATH}/get_table_head2.pl exported/Relative/otu_table.${n7}.relative.lastcolumn.txt 20 -trantab > exported/Relative/otu_table.${n7}.relative.lastcolumn.trans; 
			#perl ${SCRIPTPATH}/bar_diagram.pl -table exported/Relative/otu_table.${n7}.relative.lastcolumn.trans -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" -right -textup -rotate='-45' --y_mun 1,1 > exported/Relative/otu_table.${n7}.relative.svg;
#			perl ${SCRIPTPATH}/bar_diagram.pl -table exported/Relative/otu_table.${n7}.relative.lastcolumn.trans -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent (%)" -right -textup -rotate='-45' --y_mun 0.2,5 --micro_scale --percentage > exported/Relative/otu_table.${n7}.relative.svg
#		done;
		for svg_file in exported/Relative/*svg; do echo $svg_file; n=$(basename "$svg_file" .svg); echo $n; rsvg-convert -h 3200 -b white $svg_file > exported/Relative/${n}.png; done


		source deactivate
		source activate qiime2-2018.11
		echo "ANCOM analaysis for differential OTU"
		mkdir exported/ANCOM
		for n2 in 2 3 4 5 6 7;
			do echo $n2;
			for category_1 in $category_set;
				do echo $category_1;
					Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -t exported/collapsed/collapsed-${tax_levels[${n2}]}.qza -o media_files
					qiime composition add-pseudocount   --i-table media_files/filtered_feature_table.qza  --o-composition-table exported/ANCOM/composition.${tax_levels[${n2}]}.qza;
					qiime composition ancom  --i-table exported/ANCOM/composition.${tax_levels[${n2}]}.qza --m-metadata-file media_files/cleaned_map.txt --m-metadata-column $category_1 --o-visualization exported/ANCOM/${category_1}.ANCOM.${tax_levels[${n2}]}.qzv;
				done;
				#qiime composition ancom  --i-table exported/ANCOM/composition.${tax_levels[${n2}]}.qza --m-metadata-file $mapping_file --m-metadata-column $category_2 --o-visualization exported/ANCOM/SecondaryGroup/ANCOM.${tax_levels[${n2}]}.qzv;
		done;

		if [[ $if_picrust == 'yes' ]];then
			echo "##############################################################\n#Run for PICRUST analysis and STAMP visulization";
			qiime vsearch cluster-features-closed-reference --i-sequences rep-seqs.qza --i-table table.qza --i-reference-sequences $close_reference_trained --p-perc-identity 0.97 --p-threads 0 --output-dir closedRef_forPICRUSt
			qiime feature-table summarize --i-table closedRef_forPICRUSt/clustered_table.qza --o-visualization closedRef_forPICRUSt/clustered_table.qzv --m-sample-metadata-file $mapping_file
			qiime feature-table tabulate-seqs   --i-data closedRef_forPICRUSt/unmatched_sequences.qza   --o-visualization closedRef_forPICRUSt/unmatched_sequences.qzv
			qiime tools export --input-path closedRef_forPICRUSt/clustered_table.qza --output-path closedRef_forPICRUSt/
			biom convert -i closedRef_forPICRUSt/feature-table.biom -o closedRef_forPICRUSt/feature-table.txt --to-tsv

			export PYTHONPATH=/home/admin1/miniconda3/lib/python2.7/site-packages/
			normalize_by_copy_number.py -i closedRef_forPICRUSt/feature-table.biom -o closedRef_forPICRUSt/feature-table.normalized.biom
			predict_metagenomes.py -i closedRef_forPICRUSt/feature-table.normalized.biom -o closedRef_forPICRUSt/feature-table.metagenome.biom
			categorize_by_function.py -i closedRef_forPICRUSt/feature-table.metagenome.biom -o closedRef_forPICRUSt/feature-table.metagenome.L1.txt -c KEGG_Pathways -l 1 -f
			categorize_by_function.py -i closedRef_forPICRUSt/feature-table.metagenome.biom -o closedRef_forPICRUSt/feature-table.metagenome.L2.txt -c KEGG_Pathways -l 2 -f
			categorize_by_function.py -i closedRef_forPICRUSt/feature-table.metagenome.biom -o closedRef_forPICRUSt/feature-table.metagenome.L3.txt -c KEGG_Pathways -l 3 -f
			export PYTHONPATH=

			cd closedRef_forPICRUSt

			source deactivate
			source activate qm2

			for n3 in 1 2 3;
				do echo $n3;
				python ${SCRIPTPATH}/convert_percent.py -i feature-table.metagenome.L${n3}.txt;
				#perl ${SCRIPTPATH}/get_table_head2.pl percent.feature-table.metagenome.L${n3}.txt 35 -trantab > percent.feature-table.metagenome.L${n3}.tab
				#perl ${SCRIPTPATH}/top10_bar_diagram.pl  -right -grid -rotate='-45' -x_title 'Sample Name' -y_title 'Relative Abundance' --y_mun 0.25,4 --height 350 -table percent.feature-table.metagenome.L${n3}.tab > percent.feature-table.metagenome.L${n3}.svg
				Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category_set -i feature-table.metagenome.L${n3}.txt -o function-bar-plots-top20/ -p L${n3}_${category_set}_ -b F -s T;
				Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category_set -i feature-table.metagenome.L${n3}.txt -o function-bar-plots-top20-group-mean/ -p ${category_set}_L${n3}_mean_ -b T -s T;

				perl ${SCRIPTPATH}/cluster.pl  -BC -Z -x percent.feature-table.metagenome.L${n3}.txt > level1.relative.tree
				perl ${SCRIPTPATH}/draw_tree.pl -bun 0.25,4 -bline -type 4  level1.relative.tree  percent.feature-table.metagenome.L${n3}.tab --flank_x 100 >  tree.feature-table.metagenome.L${n3}.svg
				#${SCRIPTPATH}/PCA.R.pl ${PWD}/percent.feature-table.metagenome.L${n3}.txt 0.2 ${PWD}/PCA_L${n3}
				cp $mapping_file ./sample-metadata.PCA.txt
				perl -p -i.bak -e 's/#//' ./sample-metadata.PCA.txt
				tail -n +2 feature-table.metagenome.L${n3}.txt | grep -v "disease"> feature-table.metagenome.L${n3}.PCA.txt
				perl -p -i.bak -e 's/#OTU ID/KEGG_function/' feature-table.metagenome.L${n3}.PCA.txt
			done;
			for svg_file in *svg; do echo $svg_file; base=$(basename $svg_file .svg); rsvg-convert -h 3200 -b white $svg_file > ${base}.png; done

			cd ..
			for n5 in 1 2 3;
				do echo $n5;
				for category_1 in $category_set;do echo $category_1;Rscript ${SCRIPTPATH}/Function_PCA.r -i ${PWD}/closedRef_forPICRUSt/feature-table.metagenome.L${n5}.PCA.txt -m ${PWD}/closedRef_forPICRUSt//sample-metadata.PCA.txt -g $category_1;done;
				for category_1 in $category_set;do echo $category_1; Rscript ${SCRIPTPATH}/Function_DunnTest.r -i ${PWD}/closedRef_forPICRUSt/feature-table.metagenome.L${n5}.PCA.txt -m ${PWD}/closedRef_forPICRUSt/sample-metadata.PCA.txt -c $category_1 > /dev/null; done;
			done;
			source deactivate
			for category_1 in $category_set;
				do echo $category_1;
				function_barplot.R -i ./closedRef_forPICRUSt/feature-table.metagenome.L3.txt -m $mapping_file -c $category_1 -j T -a 0.05 -b T -o ./2-ANOVA_And_Duncan
				function_barplot.R -i ./closedRef_forPICRUSt/feature-table.metagenome.L3.txt -m $mapping_file -c $category_1 -j T -a 0.05 -b F -o ./2-ANOVA_And_Duncan
			done;

		fi;

		source deactivate
		source activate qiime2-2018.11
		echo "##############################################################\n#Make phylogenetic trees for ITOL"
<<COMMENT
		mkdir phylogeny
		qiime feature-table filter-features --i-table table.qza --p-min-frequency $min_freq --o-filtered-table phylogeny/table.${min_freq}.qza
		qiime tools export --input-path phylogeny/table.${min_freq}.qza --output-path phylogeny
		biom convert -i phylogeny/feature-table.biom -o phylogeny/feature-table.txt --to-tsv
		cut -f1 phylogeny/feature-table.txt | tail -n +3 > phylogeny/feature-table.list
		seqtk subseq exported/dna-sequences.fasta phylogeny/feature-table.list > phylogeny/dna-sequences.${min_freq}.fasta
		qiime tools import   --input-path phylogeny/dna-sequences.${min_freq}.fasta  --output-path phylogeny/dna-sequences.${min_freq}.qza   --type 'FeatureData[Sequence]'
		qiime alignment mafft   --i-sequences phylogeny/dna-sequences.${min_freq}.qza  --o-alignment phylogeny/dna-sequences.${min_freq}.aligned.qza
		qiime alignment mask   --i-alignment phylogeny/dna-sequences.${min_freq}.aligned.qza   --o-masked-alignment phylogeny/dna-sequences.${min_freq}.aligned.masked.qza
		qiime phylogeny fasttree   --i-alignment phylogeny/dna-sequences.${min_freq}.aligned.masked.qza  --o-tree phylogeny/dna-sequences.${min_freq}.unrooted-tree.qza
		qiime phylogeny midpoint-root   --i-tree phylogeny/dna-sequences.${min_freq}.unrooted-tree.qza   --o-rooted-tree phylogeny/dna-sequences.${min_freq}.rooted-tree.qza
		qiime feature-classifier classify-sklearn   --i-classifier  $reference_trained  --i-reads phylogeny/dna-sequences.${min_freq}.qza  --o-classification phylogeny/taxonomy.${min_freq}.qza

		biom add-metadata -i phylogeny/feature-table.biom -o phylogeny/feature-table.taxonomy.biom --observation-metadata-fp exported/taxonomy.tsv --observation-header OTUID,taxonomy,confidence
		biom convert -i phylogeny/feature-table.taxonomy.biom -o phylogeny/feature-table.taxonomy.txt --to-tsv --header-key taxonomy
		qiime tools export --input-path phylogeny/dna-sequences.${min_freq}.rooted-tree.qza --output-path phylogeny/
		mv phylogeny/tree.nwk phylogeny/tree.rooted.nwk
		perl ${SCRIPTPATH}/generate_file_Itol.pl phylogeny/feature-table.taxonomy.txt 
COMMENT

		qiime tools export --input-path masked-aligned-rep-seqs.qza --output-path ./
		qiime tools export --input-path rep-seqs.qza --output-path ./
		qiime tools export --input-path rooted-tree.qza --output-path ./

		echo "##############################################################\n#export all qzv files into clickable folders"
		#for f in $(find . -type f -name "*.qzv"); do echo $f; qiime tools export $f --output-dir ${f}.exported; done
		for f in $(find . -type f -name "*.qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); new=${dir}/${base}; qiime tools export --input-path $f --output-path ${new}.qzv.exported; done 


		echo "##############################################################\n#Run Qiime1 for differOTU analysis"
		source deactivate
		source activate qm1
		mkdir exported/DiffAbundance
		biom convert -i exported/Relative/otu_table.even.txt -o exported/DiffAbundance/otu_table.even.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
		summarize_taxa.py -i exported/DiffAbundance/otu_table.even.biom -a -o exported/DiffAbundance/tax
		summarize_taxa.py -i exported/DiffAbundance/otu_table.even.biom -a -L 7 -o exported/DiffAbundance/tax
		source ~/.bash_profile

		min_observation=$(echo \(`wc -l $mapping_file | sed 's/ .*//g'`-1\)/4 | bc)
		echo "###############min observation of otu in samples is $min_observation"
		mkdir tables_for_deseq_anova_kruskal
		for n4 in 2 3 4 5 6 7;
			do echo $n4;
			#the biom file should include taxonomy information for group_significance.py script
			cut -f1 exported/DiffAbundance/tax/otu_table.even_L${n4}.txt > exported/DiffAbundance/tax/otu_table.even_L${n4}.1stColumn.txt
			perl -p -i.bak -e 's/#OTU ID/taxonomy/' exported/DiffAbundance/tax/otu_table.even_L${n4}.1stColumn.txt
			paste exported/DiffAbundance/tax/otu_table.even_L${n4}.txt exported/DiffAbundance/tax/otu_table.even_L${n4}.1stColumn.txt > exported/DiffAbundance/tax/otu_table.even_${tax_levels[${n4}]}.taxonomy.txt
			biom convert -i exported/DiffAbundance/tax/otu_table.even_${tax_levels[${n4}]}.taxonomy.txt -o exported/DiffAbundance/tax/otu_table.even_${tax_levels[${n4}]}.taxonomy.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

			cp exported/DiffAbundance/tax/otu_table.even_${tax_levels[${n4}]}.taxonomy.txt tables_for_deseq_anova_kruskal/
			filter_otus_from_otu_table.py -i exported/DiffAbundance/tax/otu_table.even_${tax_levels[${n4}]}.taxonomy.biom -s $min_observation -o filtered_otu_table.biom

			for category_1 in $category_set;
				do echo $category_1;
				#Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -o media_files
				group_significance.py -i filtered_otu_table.biom -m $mapping_file -c $category_1 -s kruskal_wallis -o exported/DiffAbundance/kruskal_wallis_${category_1}_DiffAbundance_${tax_levels[${n4}]}.txt --biom_samples_are_superset --print_non_overlap;
				group_significance.py -i filtered_otu_table.biom -m $mapping_file -c $category_1 -s ANOVA -o exported/DiffAbundance/ANOVA_${category_1}_DiffAbundance_${tax_levels[${n4}]}.txt --biom_samples_are_superset --print_non_overlap;

	#			python ${SCRIPTPATH}/auto_DESeq.py -m $mapping_file -g $category_1 -l ${tax_levels[${n4}]};
				done;


			for category_1 in $category_set;
				do echo $category_1;
				#Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -o media_files
	#			group_significance.py -i filtered_otu_table.biom -m $mapping_file -c $category_1 -s kruskal_wallis -o exported/DiffAbundance/kruskal_wallis_${category_1}_DiffAbundance_${tax_levels[${n4}]}.txt --biom_samples_are_superset --print_non_overlap;
	#			group_significance.py -i filtered_otu_table.biom -m $mapping_file -c $category_1 -s ANOVA -o exported/DiffAbundance/ANOVA_${category_1}_DiffAbundance_${tax_levels[${n4}]}.txt --biom_samples_are_superset --print_non_overlap;

				python3 ${SCRIPTPATH}/auto_DESeq.py -m $mapping_file -g $category_1 -l ${tax_levels[${n4}]};
				done;

			done;



		echo "##############################################################\n#Run R script for additional R related figure generation"
		source deactivate
		source activate qm2
		mkdir R_output
		#Change format of meta-data file for Rscript of PLSDA analysis
		#cp $mapping_file ./R_output/sample-metadata.txt
		#tail -n +2 exported/feature-table.txt > ./R_output/feature-table.PLSDA.txt
		#perl -p -i.bak -e 's/#OTU ID//' ./R_output/feature-table.PLSDA.txt
		#sort ./R_output/sample-metadata.txt > ./R_output/sample-metadata.PLSDA.txt
		#Change format of meta-data file for Rscript of alpha diversity analysis
		#cp $mapping_file ./alpha/sample-metadata_alphadiversity.txt
		#perl -p -i.bak -e 's/#SampleID//' ./alpha/sample-metadata_alphadiversity.txt

		source deactivate
		source activate qm2
		for category_1 in $category_set;
			do echo $category_1;
				Rscript ${SCRIPTPATH}/clean_na_of_inputs.R -m $mapping_file --group $category_1 -o media_files
				map=$(readlink -f ./media_files/cleaned_map.txt)
				#otu=$(readlink -f ./media_files/cleaned_feature_table.txt)
				Rscript ${SCRIPTPATH}/RRelatedOutput.R $map $category_1 T;
				Rscript ${SCRIPTPATH}/alphaboxplotwitSig.R -m $map -c $category_1 -i ./alpha/alpha-summary.tsv -o ./alpha/;
			done;


		source deactivate
		source activate qm2
		Rscript ${SCRIPTPATH}/beta_heatmap.R -i exported/feature-table.ConsensusLineage.txt -m $mapping_file -t exported/tree.rooted.nwk -r exported/dna-sequences.fasta -o R_output -c $category_set -p 'unclustered_';
		Rscript ${SCRIPTPATH}/beta_heatmap.R -i exported/feature-table.ConsensusLineage.txt  -t exported/tree.rooted.nwk -r exported/dna-sequences.fasta -o R_output;

		perl ${SCRIPTPATH}/table_data_svg.pl --colors cyan-orange R_output/bray_curtis_matrix.xls R_output/weighted_unifrac_matrix.xls R_output/unweighted_unifrac_matrix.xls --symbol 'Beta Diversity' > R_output/BetaDiversity_heatmap.svg


		rsvg-convert -h 3200 -b white R_output/BetaDiversity_heatmap.svg > R_output/BetaDiversity_heatmap.png

		#python2 ${SCRIPTPATH}/biom_to_stamp.py -m KEGG_Pathways closedRef_forPICRUSt/feature-table.metagenome.biom > closedRef_forPICRUSt/feature-table.metagenome.KEGG_Pathways.STAMP.txt


		echo "##############################################################\n#Generate the absolute directory for enviromental factors relational analysis"

	#	cd exported/
		perl ${SCRIPTPATH}/stat_otu_tab.unspecifiedadded.pl -unif min exported/feature-table.taxonomy.txt --prefix exported/Absolute/otu_table -nomat -abs -spestat exported/Absolute/classified_stat.xls

	#	cd exported/Absolute/
		for key in ${!tax_aa[*]};do mv exported/Absolute/otu_table.${key}.absolute.mat exported/Absolute/otu_table.${tax_aa[$key]}.absolute.txt;done;
		#mv otu_table.k.absolute.mat otu_table.Kingdom.absolute.txt
		#mv otu_table.p.absolute.mat otu_table.Pylumn.absolute.txt
		#mv otu_table.c.absolute.mat otu_table.Class.absolute.txt
		#mv otu_table.o.absolute.mat otu_table.Order.absolute.txt
		#mv otu_table.f.absolute.mat otu_table.Family.absolute.txt
		#mv otu_table.g.absolute.mat otu_table.Genus.absolute.txt
		#mv otu_table.s.absolute.mat otu_table.Species.absolute.txt
	#	mkdir RDA
	#	for n6 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
	#		do echo $n6;
	#		mkdir RDA/${n6}
	#		cp otu_table.${n6}.absolute.txt RDA/${n6}
	#		cd RDA/${n6}
	#		for category_1 in $category_set;do echo $category_1;python ${SCRIPTPATH}/RDA.py -i otu_table.${n6}.absolute.txt -m $mapping_file -g $category_1 -o ./ -n 25 -e $not_rda;done;
	#		cd ../../
	#	done;
	#	cd ../../




		perl ${SCRIPTPATH}/stat_otu_tab.pl -unif min exported/feature-table.taxonomy.txt -prefix otu_table_forlefse/otu_table
		for key in ${!tax_aa[*]};do mv otu_table_forlefse/otu_table.${key}.relative.mat otu_table_forlefse/otu_table.${tax_aa[$key]}.relative.txt;done;

		source deactivate
		source activate qm2
<<skip
		test=${not_rda// */}
		if [ ! $test == "all" ];then
			echo "##############################################################\nCorrelation heatmap analysis"
			for nrda in $not_rda;
				do echo $nrda;
				prefix=${nrda//,/_}_excluded_;
				prefix=${prefix//none_excluded_/};
				prefix=${prefix//\//-};
				prefix=${prefix//\\/-};
				prefix=${prefix//\(/};
				prefix=${prefix//\)/};
				prefix=${prefix//\%/};
				for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
					do echo $n7;
					Rscript ${SCRIPTPATH}/cor_heatmap.R -i otu_table_forlefse/otu_table.${n7}.relative.txt -o 2-CorrelationHeatmap/${n7}/ -n 25 -m $mapping_file -e $nrda -p "$prefix";
					for category_1 in $category_set;do echo $category_1;Rscript ${SCRIPTPATH}/RDA.R -i exported/Relative/otu_table.${n7}.relative.txt -m $mapping_file -c $category_1 -o exported/Absolute/RDA/${n7} -n 25 -e $nrda -p "$prefix";done;
				done;
			done;
		fi;
skip

	    test=${not_rda// */}
	    if [ ! $test == "all" ];then
	        echo "##############################################################\nCorrelation heatmap analysis"
	        for nrda in $not_rda;
	            do echo $nrda;
	            arr=(${nrda//:/ });
	            nrda=${arr[0]};
	            prefix=${arr[1]};
	            for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
	                do echo $n7;
					Rscript ${SCRIPTPATH}/cor_heatmap.R -i otu_table_forlefse/otu_table.${n7}.relative.txt -o 2-CorrelationHeatmap/${n7}/ -n 25 -m $mapping_file -e $nrda -p "$prefix";
					for category_1 in $category_set;do echo $category_1;Rscript ${SCRIPTPATH}/RDA.R -i exported/Relative/otu_table.${n7}.relative.txt -m $mapping_file -c $category_1 -o exported/Absolute/RDA/${n7} -n 25 -e $nrda -p "$prefix";done;
	            done;
	        done;
	    fi;



		source deactivate
		source activate qm2
		echo "##############################################################network and abundance heatmap" 
		for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
			do echo $n7;
			Rscript ${SCRIPTPATH}/network.R -c 0.5 -i otu_table_forlefse/otu_table.${n7}.relative.txt -o 3-NetworkAnalysis/${n7}/;
			#Rscript ${SCRIPTPATH}/abundance_heatmap.R -n 20 -i exported/Relative/otu_table.${n7}.relative.txt -o Heatmap_top20/${n7}/;
			Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category_set -n 20 -i exported/Absolute/otu_table.${n7}.absolute.txt -o Heatmap_top20/${n7}/ -l T -t F;
			Rscript ${SCRIPTPATH}/abundance_heatmap.R -m $mapping_file -c $category_set  -n 20 -i exported/Absolute/otu_table.${n7}.absolute.txt -o Heatmap_top20_clustered/${n7}/ -l T -t F -u T;

			Rscript ${SCRIPTPATH}/abundance_heatmap.R  -m $mapping_file -c $category_set -n 20 -i exported/Absolute/otu_table.${n7}.absolute.txt -o Heatmap_top20/${n7}/ -b T -l T -p 'group_mean_' -t T;

			done;


		echo "###############################################################\nAdditional plot"
		mkdir 4-VennAndFlower
		for category_1 in $category_set;
			do echo $category_1;
			Rscript ${SCRIPTPATH}/venn_and_flower_plot.R  -i ./exported/feature-table.taxonomy.txt -m $mapping_file -c $category_1 -o ./4-VennAndFlower;
			Rscript ${SCRIPTPATH}/phylotree_and_heatmap.R -i ./exported/feature-table.ConsensusLineage.txt -t ./exported/tree.rooted.nwk -m $mapping_file -c $category_1 -o AdditionalPhylogeny/ -n 50
		done;

		source deactivate
		source activate qm2
		echo "##############################################################\n#Barplot and RDA according to group mean"
		for category_1 in $category_set;
			do echo $category_1;
			for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species"; 
				do echo $n7;
				Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category_1 -i exported/Relative/otu_table.${n7}.relative.txt -o taxa-bar-plots-top20-group-ordered/ -p ${n7}_${category_1}_ -b F;
				Rscript ${SCRIPTPATH}/abundance_barplot.R -n 20 -m $mapping_file -c $category_1 -i exported/Relative/otu_table.${n7}.relative.txt -o Barplot-of-Group-Mean/ -p ${category_1}_${n7}_mean_ -b T;
			done;
		done;

<<COMMENT
		for n6 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
			do echo $n6;
			Rscript ${SCRIPTPATH}/collapse_table_with_group_mean.R -i exported/Relative/otu_table.${n6}.relative.txt -m $mapping_file -c Group4 -o media_files -s F
			for category_1 in Group3;do echo $category_1;python ${SCRIPTPATH}/RDA.py -i media_files/abundance_table_collapsed_with_group_mean.txt -m media_files/map_collapsed_with_group_mean.txt -g $category_1 -o exported/Absolute/RDA/${n6} -n 25 -e $not_rda;done;
		done;
COMMENT

		##########alpha rarefacation
		Rscript ${SCRIPTPATH}/alphararefaction.R -i alpha-rarefaction.qzv.exported -o alpha-rarefaction-ggplot2

		Rscript ${SCRIPTPATH}/alphararefaction.R -i alpha-rarefaction.qzv.exported -o alpha-rarefaction-ggplot2 -m $mapping_file -c $category_set -p "group_mean_"

		echo "##############################################################\n#Run LEFSE for Group"

		#perl ${SCRIPTPATH}/stat_otu_tab.pl -unif min exported/feature-table.taxonomy.txt -prefix otu_table_forlefse/otu_table
		#for key in ${!tax_aa[*]};do mv otu_table_forlefse/otu_table.${key}.relative.mat otu_table_forlefse/otu_table.${tax_aa[$key]}.relative.txt;done;

		for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
			do echo $n7;
			for category_1 in $category_set;
				do echo $category_1;
				Rscript ${SCRIPTPATH}/write_data_for_lefse.R -i  otu_table_forlefse/otu_table.${n7}.relative.txt -m  $mapping_file -c  $category_1 -o  exported/Relative/Lefse/${n7}/${category_1}_${n7}_lefse.txt -u l;
			done;
		done;


		source deactivate
		source deactivate
		source activate lefse
		cd exported/Relative
		for n7 in "Phylum" "Class" "Order" "Family" "Genus" "Species";
			do echo $n7;
			cd Lefse/${n7}
			for category_1 in $category_set;
				do echo $category_1;
					base="${category_1}_${n7}_lefse_LDA2"; format_input.py ${category_1}_${n7}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 2;  
					${SCRIPTPATH}/mod_lefse-plot_res.py --category $category_1 --map $mapping_file  --max_feature_len 200 --orientation h --format pdf ${base}.LDA.txt ${base}.pdf; ${SCRIPTPATH}/mod_lefse-plot_cladogram.py ${base}.LDA.txt --map $mapping_file  --category $category_1 ${base}.cladogram.pdf --format pdf;
					base="${category_1}_${n7}_lefse_LDA4"; format_input.py ${category_1}_${n7}_lefse.txt ${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py ${base}.lefseinput.txt ${base}.LDA.txt -l 4;  
					${SCRIPTPATH}/mod_lefse-plot_res.py --category $category_1 --map $mapping_file  --max_feature_len 200 --orientation h --format pdf ${base}.LDA.txt ${base}.pdf; ${SCRIPTPATH}/mod_lefse-plot_cladogram.py ${base}.LDA.txt --map $mapping_file --category $category_1 ${base}.cladogram.pdf --format pdf;
			done;
			cd ../../
		done;
		cd ../../

	#	source deactivate
	#	source activate lefse
	#	mkdir -p LEfSe/Genus/
	#	for category_1 in $category_set;
	#		do echo $category_1;
	#			Rscript ${SCRIPTPATH}/write_data_for_lefse.R  exported/Absolute/otu_table.Genus.absolute.txt  $mapping_file  $category_1  LEfSe/Genus/${category_1}_table_for_lefse.txt F;
	#			base="${category_1}_Genus_LEfSe_LDA2"; format_input.py LEfSe/Genus/${category_1}_table_for_lefse.txt LEfSe/Genus/${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py LEfSe/Genus/${base}.lefseinput.txt LEfSe/Genus/${base}.LDA.txt -l 2;  ${SCRIPTPATH}/mod_lefse-plot_res.py --dpi 300 LEfSe/Genus/${base}.LDA.txt LEfSe/Genus/${base}.png; ${SCRIPTPATH}/mod_lefse-plot_cladogram.py LEfSe/Genus/${base}.LDA.txt --dpi 300 LEfSe/Genus/${base}.cladogram.png --format png   --labeled_stop_lev 4;
	#			${SCRIPTPATH}/mod_lefse-plot_res.py  --max_feature_len 200 --orientation h --format pdf --dpi 300 LEfSe/Genus/${base}.LDA.txt LEfSe/Genus/${base}.pdf; ${SCRIPTPATH}/mod_lefse-plot_cladogram.py LEfSe/Genus/${base}.LDA.txt --dpi 300 LEfSe/Genus/${base}.cladogram.pdf  --format pdf   --labeled_stop_lev 4;
	#		done;

	#	mkdir -p LEfSe/OTU/
	#	for category_1 in $category_set;
	#		do echo $category_1;
	#			Rscript ${SCRIPTPATH}/write_data_for_lefse.R  exported/feature-table.taxonomy.txt  $mapping_file  $category_1  LEfSe/OTU/${category_1}_table_for_lefse.txt T;
				#Rscript ${SCRIPTPATH}/write_data_for_lefse.R  exported/Absolute/otu_table.Genus.absolute.txt  $mapping_file  $category_1  LEfSe/OTU/${category_1}_table_for_lefse.txt F;
	#			base="${category_1}_OTU_LEfSe_LDA4"; format_input.py LEfSe/OTU/${category_1}_table_for_lefse.txt LEfSe/OTU/${base}.lefseinput.txt -c 2 -u 1 -o 1000000; run_lefse.py LEfSe/OTU/${base}.lefseinput.txt LEfSe/OTU/${base}.LDA.txt -l 4;  ${SCRIPTPATH}/mod_lefse-plot_res.py --dpi 300 LEfSe/OTU/${base}.LDA.txt LEfSe/OTU/${base}.png; ${SCRIPTPATH}/mod_lefse-plot_cladogram.py LEfSe/OTU/${base}.LDA.txt --dpi 300 LEfSe/OTU/${base}.cladogram.png --format png  ;
	#			${SCRIPTPATH}/mod_lefse-plot_res.py  --max_feature_len 200 --orientation h --format pdf --dpi 300 LEfSe/OTU/${base}.LDA.txt LEfSe/OTU/${base}.pdf; ${SCRIPTPATH}/mod_lefse-plot_cladogram.py LEfSe/OTU/${base}.LDA.txt --dpi 300 LEfSe/OTU/${base}.cladogram.pdf  --format pdf  ;
	#		done;



		source deactivate
		source activate qm2
		category_report=($category_set)
		echo "##############################################################\n#Organize the result files";
		#cp -r ${SCRIPTPATH}/Result ./
		bash ${SCRIPTPATH}/organize_dir_structure_V2.sh $mapping_file $category_report ${SCRIPTPATH} wu $if_picrust $prefix;
		if [ -d "../Result/${category_report}" ];then
			rm -r ../Result/${category_report};
		fi;
		mv Result ../Result/${category_report}


		cd ../;
	done;
