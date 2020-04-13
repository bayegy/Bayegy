#/usr/bin/env bash
#########
#Please address any bugs to Cheng. 
#Date 2017.12.19
#########
mapping_file=$1
category_1=$2
SCRIPTPATH=$3
number=$4
if_picrust=$5
prefix=$6


if [ -z "$3" ]; then
	echo "##########

		  Please prepare the directory structure before starting the program like below:
		  raw/fastq_files ...
		  mapping_file
		  manifest_file
		  \n\n"

	echo "Please provide following input parameters
		1) Full path of the mapping file. (Accept both .csv or txt format.)
		2) The name of the first category in the mapping file. 

		Sample Usage:
		sh $0 M231_Mapping_2.tsv Group1 readme.pdf
		"
	exit 0
else
	echo "################
	Running: sh $0 $1 $2 $3 $4 $5 $6"
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

	echo "##############################################################\n#Organize the Result folder"
	SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
	if [ -d "./Result" ];then
		rm -r ./Result;
	fi;
	mkdir -p ./Result/ \
	./Result/2-AbundanceAnalysis/ ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/ \
	./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/ \
	./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/3-CollapsedStats/ \
	./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/ \
	./Result/2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/ \
	./Result/2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/ \
	./Result/2-AbundanceAnalysis/2-AbundanceComparison/2-ANOVA/ \
	./Result/2-AbundanceAnalysis/2-AbundanceComparison/3-KruskalWallis/ \
	./Result/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/ \
	./Result/2-AbundanceAnalysis/2-AbundanceComparison/5-DESeq2/ \
	./Result/3-AlphaDiversity/1-AlphaDiversitySummary/ \
	./Result/3-AlphaDiversity/2-AlphaRarefaction/ \
	./Result/3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/ \
	./Result/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/ \
	./Result/4-BetaDiversity/1-BetaDiversitySummary/ \
	./Result/4-BetaDiversity/2-PCoA/ \
	./Result/4-BetaDiversity/3-NMDS/ \
	./Result/4-BetaDiversity/4-PLS-DA/ \
	./Result/4-BetaDiversity/5-GroupSignificance/ \
	./Result/4-BetaDiversity/6-PCA/ \
	./Result/5-Phylogenetics/ \
	./Result/6-AssociationAnalysis/1-RDA/ \
	./Result/FiguresTablesForReport \
 	./Result/4-BetaDiversity/2-PCoA/PCoA-Phyloseq \
	./Result/4-BetaDiversity/2-PCoA/PCoA-Qiime2


	if [[ $if_picrust == 'yes' ]]; then
		mkdir mkdir -p ./Result/7-FunctionAnalysis/1-KEGG_Pathway/ \
		./Result/7-FunctionAnalysis/2-PCAPlots/ \
		./Result/7-FunctionAnalysis/3-SignifcanceAnalysis/ \
		./Result/7-FunctionAnalysis/3-SignifcanceAnalysis/1-DunnTest
	fi

	echo "Start organize the files for deliverables ..."
	cp ${SCRIPTPATH}/Result_README.pdf ./Result/
	cp $mapping_file ./Result/mapping_file.txt

	cp ${category_1}_classified_stat_relative.xls ${category_1}_classified_stat_relative.png ./Result/2-AbundanceAnalysis/

	cp exported/Absolute/*absolute.txt ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/
	#cp exported/Absolute/otu_table.p.absolute.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Phylum.absolute.txt
	#cp exported/Absolute/otu_table.c.absolute.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Class.absolute.txt
	#cp exported/Absolute/otu_table.o.absolute.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Order.absolute.txt
	#cp exported/Absolute/otu_table.f.absolute.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Family.absolute.txt
	#cp exported/Absolute/otu_table.g.absolute.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Genus.absolute.txt
	#cp exported/Absolute/otu_table.s.absolute.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Species.absolute.txt
	cp exported/Relative/*relative.txt ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/
	#cp exported/Relative/otu_table.p.relative.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Phylum.relative.txt
	#cp exported/Relative/otu_table.c.relative.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Class.relative.txt
	#cp exported/Relative/otu_table.o.relative.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Order.relative.txt
	#cp exported/Relative/otu_table.f.relative.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Family.relative.txt
	#cp exported/Relative/otu_table.g.relative.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Genus.relative.txt
	#cp exported/Relative/otu_table.s.relative.mat ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Species.relative.txt

	cp -r exported/collapsed/*qzv* ./Result/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/3-CollapsedStats/

	cp -r taxa-bar-plots.qzv* ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/
	#cp -r exported/Relative/*relative.txt exported/Relative/otu*png ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots-top20


	# cp -r exported/${number}/*.qzv* ./Result/2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/Heatmap-Qiime2/
	cp -r Heatmap_top20 ./Result/2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/
	cp -r Heatmap_top20_clustered ./Result/2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/

	cp -r exported/ANCOM/*.qzv* ./Result/2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/

	cp TaxaStatistic/* ./Result/2-AbundanceAnalysis/

	cp -r exported/DiffAbundance/ANOVA_*txt ./Result/2-AbundanceAnalysis/2-AbundanceComparison/2-ANOVA/

	cp -r exported/DiffAbundance/kruskal_wallis_*txt ./Result/2-AbundanceAnalysis/2-AbundanceComparison/3-KruskalWallis/

	#cp lefse result ./Result/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/
	cp -r exported/Relative/Lefse/*/ ./Result/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/
	cp -r exported/DiffAbundance/DESeq2_*txt ./Result/2-AbundanceAnalysis/2-AbundanceComparison/5-DESeq2/
	cp -rp tables_for_deseq_anova_kruskal ./Result/2-AbundanceAnalysis/2-AbundanceComparison/

	cp -r alpha/alpha-summary.tsv R_output/*alpha_diversity_* ./Result/3-AlphaDiversity/1-AlphaDiversitySummary/

	cp -r alpha-rarefaction.qzv* ./Result/3-AlphaDiversity/2-AlphaRarefaction/

	cp -r alpha/*wilcox* ./Result/3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/

	#cp -r core-metrics-results/observed*qzv* ./Result/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/

	#cp -r core-metrics-results/shannon*qzv*  ./Result/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/

	#cp -r core-metrics-results/faith*qzv* ./Result/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/
	cp -r alpha_groupsignificance/* ./Result/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/



	cp -r R_output/*matrix.xls R_output/BetaDiversity_heatmap.png ./Result/4-BetaDiversity/1-BetaDiversitySummary/

	cp -r R_output/*summary.pdf ./Result/4-BetaDiversity/1-BetaDiversitySummary/
	cp -r core-metrics-results/*_emperor.qzv* ./Result/4-BetaDiversity/2-PCoA/PCoA-Qiime2
	cp -r R_output/*PCoA* ./Result/4-BetaDiversity/2-PCoA/PCoA-Phyloseq
	cp -r R_output/*NMDS* ./Result/4-BetaDiversity/3-NMDS/
	cp -r core-metrics-results/*permanova*significance.qzv* ./Result/4-BetaDiversity/5-GroupSignificance/
	cp -r core-metrics-results/*anosim*significance.qzv* ./Result/4-BetaDiversity/5-GroupSignificance/
	#cp -r R_output/unifrac*summary.pdf ./Result/4-BetaDiversity/1-BetaDiversitySummary/
	#cp -r core-metrics-results/unweighted*_emperor.qzv* R_output/*unifrac*PCoA* ./Result/4-BetaDiversity/2-PCoA/
	#cp -r R_output/*unifrac*NMDS* ./Result/4-BetaDiversity/3-NMDS/
	#cp -r core-metrics-results/unweighted*significance.qzv* ./Result/4-BetaDiversity/5-GroupSignificance/

	#cp -r R_output/wunifrac*summary.pdf ./Result/4-BetaDiversity/1-BetaDiversitySummary/
	#cp -r core-metrics-results/weighted*_emperor.qzv* R_output/*wunifrac*PCoA* ./Result/4-BetaDiversity/2-PCoA/
	#cp -r R_output/*wunifrac*NMDS* ./Result/4-BetaDiversity/3-NMDS/
	cp -r R_output/*PLSDA* ./Result/4-BetaDiversity/4-PLS-DA/
	cp -r R_output/*PCA* ./Result/4-BetaDiversity/6-PCA/
	#cp -r core-metrics-results/weighted*significance.qzv* ./Result/4-BetaDiversity/5-GroupSignificance/
	cp -r exported/feature-table.taxonomy.txt ./Result/2-AbundanceAnalysis/${category_1}_feature-table.taxonomy.txt
	cp dna-sequences.fasta ./Result/5-Phylogenetics/${category_1}_representative-sequence.fasta
	cp tree.nwk ./Result/5-Phylogenetics/${category_1}_rooted-tree.nwk
	#cp -r R_output/Bacteria.phylogeny.pdf ./Result/5-Phylogenetics/1-MajorPhylums/
	#cp -r phylogeny/tol_* phylogeny/tree.rooted.nwk ./Result/5-Phylogenetics/2-MajorOTUs/

	#cp -r phylogeny/tol_* phylogeny/tree.rooted.nwk ./Result/5-Phylogenetics/


	cp -r exported/Absolute/RDA/* ./Result/6-AssociationAnalysis/1-RDA/
	#rm ./Result/6-AssociationAnalysis/1-RDA/*/data.txt ./Result/6-AssociationAnalysis/1-RDA/*/rda.R
	#cp -r exported/Absolute/p.rda.pdf ./Result/6-AssociationAnalysis/1-RDA/Phylum.rda.pdf
	#cp -r exported/Absolute/c.rda.pdf ./Result/6-AssociationAnalysis/1-RDA/Class.rda.pdf
	#cp -r exported/Absolute/o.rda.pdf ./Result/6-AssociationAnalysis/1-RDA/Order.rda.pdf
	#cp -r exported/Absolute/f.rda.pdf ./Result/6-AssociationAnalysis/1-RDA/Family.rda.pdf
	#cp -r exported/Absolute/g.rda.pdf ./Result/6-AssociationAnalysis/1-RDA/Genus.rda.pdf
	#cp -r exported/Absolute/s.rda.pdf ./Result/6-AssociationAnalysis/1-RDA/Species.rda.pdf

	#cp -r /Result/6-AssociationAnalysis
	if [[ $if_picrust == 'yes' ]]; then
		cp -rp 2-ANOVA_And_Duncan/ ./Result/7-FunctionAnalysis/3-SignifcanceAnalysis/
		cp -r closedRef_forPICRUSt/feature-table.metagenome.L1.txt closedRef_forPICRUSt/feature-table.metagenome.L2.txt closedRef_forPICRUSt/feature-table.metagenome.L3.txt ./Result/7-FunctionAnalysis/1-KEGG_Pathway/
		#rm ./Result/7-FunctionAnalysis/1-KEGG_Pathway/*PCA* ./Result/7-FunctionAnalysis/1-KEGG_Pathway/*DunnTest*
		cp -rp closedRef_forPICRUSt/function-bar-plots-top20-group-mean closedRef_forPICRUSt/function-bar-plots-top20  ./Result/7-FunctionAnalysis/1-KEGG_Pathway/

		cp -r closedRef_forPICRUSt/*PCA*pdf ./Result/7-FunctionAnalysis/2-PCAPlots/
		cp -r closedRef_forPICRUSt/*DunnTest*txt ./Result/7-FunctionAnalysis/3-SignifcanceAnalysis/1-DunnTest/
	fi
	#mv ./Result/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.PCA.txt.PCA.pdf ./Result/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.PCA.pdf
	#mv ./Result/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L2.PCA.txt.PCA.pdf ./Result/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L2.PCA.pdf
	#mv ./Result/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L3.PCA.txt.PCA.pdf ./Result/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L3.PCA.pdf
	#rm ./Result/7-FunctionAnalysis/2-PCAPlots/PCA*/PCA.R ./Result/7-FunctionAnalysis/2-PCAPlots/PCA*/Rplots.pdf

	cp AdditionalPhylogeny/*.pdf AdditionalPhylogeny/*table.txt ./Result/5-Phylogenetics/
	cp -rp 4-VennAndFlower/ ./Result/1-VennAndFlower
	rm ./Result/1-VennAndFlower/*.log
	cp -rp 3-NetworkAnalysis/ ./Result/6-AssociationAnalysis/
	cp -rp 2-CorrelationHeatmap/ ./Result/6-AssociationAnalysis/

	cp -rp alpha-rarefaction-ggplot2/ ./Result/3-AlphaDiversity/2-AlphaRarefaction/
	cp -rp Barplot-of-Group-Mean/ ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots-top20-group-mean
	cp -rp taxa-bar-plots-top20-group-ordered/ ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots-top20
	#rm ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots-top20-group-mean/*svg
	#change index.html to a more obvious name, and organize the qzv.exported and qzv files.
	#cd ./Result/
	for f in $(find ./Result/ -type f -name "*qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); mv $f ${f}.exported; mv ${f}.exported ${dir}/${base}; done
	for f in $(find ./Result/ -type f -name "index.html") ; do echo $f; base=$(basename $f .html); dir=$(dirname $f); new=${dir}/Summary_请点此文件查看.html; mv $f $new; done
	#cd ../



	#minor adjustment of file structure
	mv ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots/ ./Result/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots_Qiime2
	mv ./Result/3-AlphaDiversity/2-AlphaRarefaction/alpha-rarefaction ./Result/3-AlphaDiversity/2-AlphaRarefaction/alpha-rarefaction-Qiime2
	rm -r ./Result/4-BetaDiversity/2-PCoA/PCoA-Qiime2/jaccard_emperor
	###rename the pcoa results
	#mv ./Result/4-BetaDiversity/*emperor* ./Result/4-BetaDiversity/2-PCoA/1-Plots-from-QIIME2
	#mv ./Result/4-BetaDiversity/*PCoA* ./Result/4-BetaDiversity/2-PCoA/2-Plots-from-R
	#######################For FiguresTablesForReport

	cp -rp ${SCRIPTPATH}/Report/src Result/FiguresTablesForReport/
	cp ${SCRIPTPATH}/Report/结题报告.html Result/

	cd ./Result/FiguresTablesForReport
	#cp -rp ../2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/Heatmap-Qiime2/${category_1}-table-Phylum.${number}/ page4
	cp ../2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/Heatmap_top20_clustered/Phylum/heatmap.pdf Figure4-3.pdf
	cp -rp ../2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/${category_1}.ANCOM.Genus/ page4-2
	cp -rp ../4-BetaDiversity/5-GroupSignificance/unweighted_unifrac-permanova-${category_1}-significance/ page6-2
	#cp ../1-QCStats/1-Stats-demux/demultiplex-summary.png Figure3-1.png
	#cp ../2-AbundanceAnalysis/Classified_stat_relative.png Figure4-1.png
	cp ../2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots-top20-group-mean/${category_1}_Phylum_mean_barplot.pdf Figure4-2.pdf
	cp ../2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/Genus/${category_1}_Genus_lefse_LDA2.pdf Figure4-4.pdf
	cp ../2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/Genus/${category_1}_Genus_lefse_LDA2.cladogram.pdf Figure4-5.pdf
	cp ../1-VennAndFlower/${category_1}_Venn_plot.png Figure4-6.png
	cp ../3-AlphaDiversity/1-AlphaDiversitySummary/${category_1}_alpha_diversity_shannon.boxplot.pdf Figure5-1.pdf
	cp ../3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/shannon_${category_1}_wilcox_compare_boxplot.png Figure5-2.png
	cp ../4-BetaDiversity/1-BetaDiversitySummary/BetaDiversity_heatmap.png Figure6-1.png
	cp ../4-BetaDiversity/3-NMDS/${category_1}_unweighted_unifrac_NMDS_without_labels.pdf Figure6-2.pdf
	cp ../5-Phylogenetics/${category_1}_phylogenetic_tree_heatmap.pdf Figure7-1.pdf
	cp ../6-AssociationAnalysis/1-RDA/Genus/${prefix}${category_1}*features_location_plot.pdf Figure8-1.pdf
	cp ../6-AssociationAnalysis/2-CorrelationHeatmap/Genus/${prefix}Correlation_heatmap.pdf Figure8-2.pdf
	cp ../6-AssociationAnalysis/3-NetworkAnalysis/Genus/Correlation_network.pdf Figure8-3.pdf

	if [[ $if_picrust == 'yes' ]]; then
		cp ../7-FunctionAnalysis/1-KEGG_Pathway/function-bar-plots-top20-group-mean/${category_1}_L1_mean_*.pdf Figure9-1.pdf
		cp ../7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.${category_1}.PCA.pdf Figure9-2.pdf
		cp ../7-FunctionAnalysis/3-SignifcanceAnalysis/2-ANOVA_And_Duncan/${category_1}_all_significant_pathway_barplot_of_duncan.pdf Figure9-3.pdf
	fi


	if [ -f Figure4-2.pdf ];then echo "Converting pdf to png"; for pdfs in *.pdf; do echo $pdfs; base=$(basename $pdfs .pdf); convert  -density 300 -quality 80 $pdfs ${base}.png; rm $pdfs;done;fi;

	python ~/github/Bayegy/convert_to_html_table.py -i ../../../Result/1-OTUStats/2-Stats-dada2/Summary_请点此文件查看.html -o src/pages/main_cleaned.html -t dada2html -k '无法加载表格表3-1'