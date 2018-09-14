#/bin/sh -S
#########
#Please address any bugs to Cheng. 
#Date 2017.12.19
#########
mapping_file=$1
category_1=$2
SCRIPTPATH=$3
number=$4

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
	Running: sh $0 $1 $2 $3"
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

organize_deliverable_structure() {

	mkdir -p ./Result_AmpliconSequencing/ ./Result_AmpliconSequencing/1-QCStats/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/ ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/3-CollapsedStats/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/2-RepresentiveSequence/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/4-Heatmaps/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/2-ANOVA/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/3-KruskalWallis/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/5-DESeq2/ \
	./Result_AmpliconSequencing/3-AlphaDiversity/1-AlphaDiversitySummary/ \
	./Result_AmpliconSequencing/3-AlphaDiversity/2-AlphaRarefaction/ \
	./Result_AmpliconSequencing/3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/ \
	./Result_AmpliconSequencing/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/ \
	./Result_AmpliconSequencing/4-BetaDiversity/1-BetaDiversitySummary/ \
	./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/ \
	./Result_AmpliconSequencing/4-BetaDiversity/3-NMDS/ \
	./Result_AmpliconSequencing/4-BetaDiversity/4-PLS-DA/ \
	./Result_AmpliconSequencing/4-BetaDiversity/5-GroupSignificance/ \
	./Result_AmpliconSequencing/5-Phylogenetics/ \
	./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/ \
	./Result_AmpliconSequencing/7-FunctionAnalysis/1-KEGG_Pathway/ \
	./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/ \
	./Result_AmpliconSequencing/7-FunctionAnalysis/3-TreeBasedPlots/ \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/ \
	./Result_AmpliconSequencing/8-FiguresTablesForReport \
 	./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/PCoA-Phyloseq \
	./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/PCoA-Qiime2 \
	./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/taxa-bar-plots-top20 


	echo "Start organize the files for deliverables ..."
	
	cp ${SCRIPTPATH}/Result_README.pdf ./Result_AmpliconSequencing/
	cp $mapping_file ./Result_AmpliconSequencing/mapping_file.txt

	cp -r demux.qzv* stats-dada2.qzv* ./Result_AmpliconSequencing/1-QCStats/

	cp -r exported/feature-table.taxonomy.txt exported/feature-table.taxonomy.biom exported/Relative/Classified_stat_relative.png ./Result_AmpliconSequencing/2-AbundanceAnalysis/
	cp -r exported/Relative/otu_table.even.txt ./Result_AmpliconSequencing/2-AbundanceAnalysis/feature-table.taxonomy.even.txt

	cp exported/Absolute/*absolute.txt ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/
	#cp exported/Absolute/otu_table.p.absolute.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Phylum.absolute.txt
	#cp exported/Absolute/otu_table.c.absolute.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Class.absolute.txt
	#cp exported/Absolute/otu_table.o.absolute.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Order.absolute.txt
	#cp exported/Absolute/otu_table.f.absolute.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Family.absolute.txt
	#cp exported/Absolute/otu_table.g.absolute.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Genus.absolute.txt
	#cp exported/Absolute/otu_table.s.absolute.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/1-Absolute/feature-table.Species.absolute.txt
	cp exported/Relative/*relative.txt ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/
	#cp exported/Relative/otu_table.p.relative.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Phylum.relative.txt
	#cp exported/Relative/otu_table.c.relative.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Class.relative.txt
	#cp exported/Relative/otu_table.o.relative.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Order.relative.txt
	#cp exported/Relative/otu_table.f.relative.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Family.relative.txt
	#cp exported/Relative/otu_table.g.relative.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Genus.relative.txt
	#cp exported/Relative/otu_table.s.relative.mat ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/2-Relative/feature-table.Species.relative.txt

	cp -r exported/collapsed/*qzv* ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/3-CollapsedStats/

	cp -r ./rep-seqs.qzv* ./exported/*nwk ./exported/dna-sequences.fasta ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/2-RepresentiveSequence/

	cp -r taxa-bar-plots.qzv* ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/
	cp -r exported/Relative/*relative.txt exported/Relative/otu*png ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/taxa-bar-plots-top20

	cp -r exported/${number}/*.qzv* ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/4-Heatmaps/

	cp -r exported/ANCOM/*.qzv* ./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/

	cp -r exported/DiffAbundance/ANOVA_*txt ./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/2-ANOVA/

	cp -r exported/DiffAbundance/kruskal_wallis_*txt ./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/3-KruskalWallis/

	#cp lefse result ./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/
	cp -r exported/Relative/Lefse/*/ ./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/
	cp -r exported/DiffAbundance/DESeq2_*txt ./Result_AmpliconSequencing/2-AbundanceAnalysis/2-AbundanceComparison/5-DESeq2/

	cp -r alpha/alpha-summary.tsv R_output/*alpha_diversity_* ./Result_AmpliconSequencing/3-AlphaDiversity/1-AlphaDiversitySummary/

	cp -r alpha-rarefaction.qzv* ./Result_AmpliconSequencing/3-AlphaDiversity/2-AlphaRarefaction/

	cp -r alpha/*wilcox* ./Result_AmpliconSequencing/3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/

	cp -r core-metrics-results/observed*qzv* ./Result_AmpliconSequencing/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/

	cp -r core-metrics-results/shannon*qzv*  ./Result_AmpliconSequencing/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/

	cp -r core-metrics-results/faith*qzv* ./Result_AmpliconSequencing/3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/

	cp -r R_output/*matrix.txt R_output/BetaDiversity_heatmap.png ./Result_AmpliconSequencing/4-BetaDiversity/1-BetaDiversitySummary/

	cp -r R_output/*summary.pdf ./Result_AmpliconSequencing/4-BetaDiversity/1-BetaDiversitySummary/
	cp -r core-metrics-results/*_emperor.qzv* ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/PCoA-Qiime2
	cp -r R_output/*PCoA* ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/PCoA-Phyloseq
	cp -r R_output/*NMDS* ./Result_AmpliconSequencing/4-BetaDiversity/3-NMDS/
	cp -r core-metrics-results/*permanova*significance.qzv* ./Result_AmpliconSequencing/4-BetaDiversity/5-GroupSignificance/
	cp -r core-metrics-results/*anosim*significance.qzv* ./Result_AmpliconSequencing/4-BetaDiversity/5-GroupSignificance/

	#cp -r R_output/unifrac*summary.pdf ./Result_AmpliconSequencing/4-BetaDiversity/1-BetaDiversitySummary/
	#cp -r core-metrics-results/unweighted*_emperor.qzv* R_output/*unifrac*PCoA* ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/
	#cp -r R_output/*unifrac*NMDS* ./Result_AmpliconSequencing/4-BetaDiversity/3-NMDS/
	#cp -r core-metrics-results/unweighted*significance.qzv* ./Result_AmpliconSequencing/4-BetaDiversity/5-GroupSignificance/

	#cp -r R_output/wunifrac*summary.pdf ./Result_AmpliconSequencing/4-BetaDiversity/1-BetaDiversitySummary/
	#cp -r core-metrics-results/weighted*_emperor.qzv* R_output/*wunifrac*PCoA* ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/
	#cp -r R_output/*wunifrac*NMDS* ./Result_AmpliconSequencing/4-BetaDiversity/3-NMDS/
	cp -r R_output/*PCA* R_output/*PLSDA* ./Result_AmpliconSequencing/4-BetaDiversity/4-PLS-DA/
	#cp -r core-metrics-results/weighted*significance.qzv* ./Result_AmpliconSequencing/4-BetaDiversity/5-GroupSignificance/
	
	#cp -r R_output/Bacteria.phylogeny.pdf ./Result_AmpliconSequencing/5-Phylogenetics/1-MajorPhylums/
	#cp -r phylogeny/tol_* phylogeny/tree.rooted.nwk ./Result_AmpliconSequencing/5-Phylogenetics/2-MajorOTUs/
	cp -r phylogeny/tol_* phylogeny/tree.rooted.nwk ./Result_AmpliconSequencing/5-Phylogenetics/


	cp -r exported/Absolute/RDA/* ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/
	#rm ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/*/data.txt ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/*/rda.R
	#cp -r exported/Absolute/p.rda.pdf ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/Phylum.rda.pdf
	#cp -r exported/Absolute/c.rda.pdf ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/Class.rda.pdf
	#cp -r exported/Absolute/o.rda.pdf ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/Order.rda.pdf
	#cp -r exported/Absolute/f.rda.pdf ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/Family.rda.pdf
	#cp -r exported/Absolute/g.rda.pdf ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/Genus.rda.pdf
	#cp -r exported/Absolute/s.rda.pdf ./Result_AmpliconSequencing/6-AssociationAnalysis/1-RDA/Species.rda.pdf

	#cp -r /Result_AmpliconSequencing/6-AssociationAnalysis

	cp -r closedRef_forPICRUSt/feature-table.metagenome* closedRef_forPICRUSt/percent.feature-table.metagenome*png ./Result_AmpliconSequencing/7-FunctionAnalysis/1-KEGG_Pathway/
	rm ./Result_AmpliconSequencing/7-FunctionAnalysis/1-KEGG_Pathway/*PCA*
	cp -r closedRef_forPICRUSt/*PCA*pdf ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/
	#mv ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.PCA.txt.PCA.pdf ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.PCA.pdf
	#mv ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L2.PCA.txt.PCA.pdf ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L2.PCA.pdf
	#mv ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L3.PCA.txt.PCA.pdf ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L3.PCA.pdf
	#rm ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/PCA*/PCA.R ./Result_AmpliconSequencing/7-FunctionAnalysis/2-PCAPlots/PCA*/Rplots.pdf

	cp -r closedRef_forPICRUSt/tree*png ./Result_AmpliconSequencing/7-FunctionAnalysis/3-TreeBasedPlots/

	cp AdditionalPhylogeny/*.pdf ./Result_AmpliconSequencing/5-Phylogenetics/
	cp -rp 4-VennAndFlower/ ./Result_AmpliconSequencing/2-AbundanceAnalysis/3-VennAndFlower
	rm ./Result_AmpliconSequencing/2-AbundanceAnalysis/3-VennAndFlower/*.log
	cp -rp 3-NetworkAnalysis/ ./Result_AmpliconSequencing/6-AssociationAnalysis/
	cp -rp 2-CorrelationHeatmap/ ./Result_AmpliconSequencing/6-AssociationAnalysis/
	cp -rp alpha-rarefaction-ggplot2/ ./Result_AmpliconSequencing/3-AlphaDiversity/2-AlphaRarefaction/
	
	#change index.html to a more obvious name, and organize the qzv.exported and qzv files.
	#cd ./Result_AmpliconSequencing/
	for f in $(find ./Result_AmpliconSequencing/ -type f -name "*qzv"); do echo $f; base=$(basename $f .qzv); dir=$(dirname $f); mv $f ${f}.exported; mv ${f}.exported ${dir}/${base}; done
	for f in $(find ./Result_AmpliconSequencing/ -type f -name "index.html") ; do echo $f; base=$(basename $f .html); dir=$(dirname $f); new=${dir}/Summary_请点此文件查看.html; mv $f $new; done
	#cd ../



	#minor adjustment of file structure
	mv ./Result_AmpliconSequencing/1-QCStats/demux/ ./Result_AmpliconSequencing/1-QCStats/1-Stats-demux
	mv ./Result_AmpliconSequencing/1-QCStats/stats-dada2 ./Result_AmpliconSequencing/1-QCStats/2-Stats-dada2
	mv ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/taxa-bar-plots/ ./Result_AmpliconSequencing/2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/taxa-bar-plots_Qiime2
	mv ./Result_AmpliconSequencing/3-AlphaDiversity/2-AlphaRarefaction/alpha-rarefaction ./Result_AmpliconSequencing/3-AlphaDiversity/2-AlphaRarefaction/alpha-rarefaction-Qiime2
	rm -r ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/PCoA-Qiime2/jaccard_emperor
	###rename the pcoa results
	#mv ./Result_AmpliconSequencing/4-BetaDiversity/*emperor* ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/1-Plots-from-QIIME2
	#mv ./Result_AmpliconSequencing/4-BetaDiversity/*PCoA* ./Result_AmpliconSequencing/4-BetaDiversity/2-PCoA/2-Plots-from-R
	#######################For 8-FiguresTablesForReport

	cp -rp ${SCRIPTPATH}/Report/src Result_AmpliconSequencing/8-FiguresTablesForReport/
	cp ${SCRIPTPATH}/Report/结题报告.html Result_AmpliconSequencing/



	cd ./Result_AmpliconSequencing/8-FiguresTablesForReport
	cp -rp ../2-AbundanceAnalysis/1-AbundanceSummary/4-Heatmaps/${category_1}-table-Phylum.${number}/ page4
	cp -rp ../2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/${category_1}.ANCOM.Genus/ page4-2
	cp -rp ../4-BetaDiversity/5-GroupSignificance/unweighted_unifrac-permanova-${category_1}-significance/ page6-2
	cp ../1-QCStats/1-Stats-demux/demultiplex-summary.png Figure3-1.png
	cp ../2-AbundanceAnalysis/Classified_stat_relative.png Figure4-1.png
	cp ../2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/taxa-bar-plots-top20/otu_table.Phylum.relative.png Figure4-2.png
	cp ../2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/Phylum/${category_1}_Phylum_lefse.png Figure4-3.png
	cp ../2-AbundanceAnalysis/3-VennAndFlower/${category_1}_Venn_plot.png Figure4-4.png
	cp ../3-AlphaDiversity/1-AlphaDiversitySummary/${category_1}_alpha_diversity_shannon.boxplot.pdf Figure5-1.pdf
	cp ../3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/shannon_${category_1}_wilcox_compare_boxplot.png Figure5-2.png
	cp ../4-BetaDiversity/1-BetaDiversitySummary/BetaDiversity_heatmap.png Figure6-1.png
	cp ../4-BetaDiversity/3-NMDS/${category_1}_unweighted_unifrac_NMDS.pdf Figure6-2.pdf
	cp ../5-Phylogenetics/${category_1}_phylogenetic_tree_heatmap.pdf Figure7-1.pdf
	cp ../6-AssociationAnalysis/1-RDA/Genus/*_bacteria_location_plot.png Figure8-1.png
	cp ../6-AssociationAnalysis/2-CorrelationHeatmap/Genus/Correlation_heatmap.pdf Figure8-2.pdf
	cp ../6-AssociationAnalysis/3-NetworkAnalysis/Genus/Correlation_network.pdf Figure8-3.pdf
	cp ../7-FunctionAnalysis/1-KEGG_Pathway/percent.feature-table.metagenome.L1.png Figure9-1.png
	cp ../7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.${category_1}.PCA.pdf Figure9-2.pdf
	cp ../7-FunctionAnalysis/3-TreeBasedPlots/tree.feature-table.metagenome.L1.png Figure9-3.png



	#cp ../2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/${category_1}.ANCOM.Genus/${category_1}.ANCOM.Genus.qzv Table3-1.qzv
	#cp ../2-AbundanceAnalysis/1-AbundanceSummary/1-AbundanceTable/3-CollapsedStats/collapsed-Species/collapsed-Species.qzv  Table3-2.qzv
	#cp ../3-AlphaDiversity/1-AlphaDiversitySummary/alpha-summary.tsv Table3-3.txt
	#cp ../3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/shannon-group-significance/shannon-group-significance.qzv Table3-4.qzv
	#cp ../4-BetaDiversity/5-GroupSignificance/unweighted_unifrac-permanova-${category_1}-significance/unweighted_unifrac-permanova-${category_1}-significance.qzv Table3-5.qzv
	
	#cp ../1-QCStats/1-Stats-demux/demux.qzv Figure2-1.qzv
	
	#cp ../2-AbundanceAnalysis/1-AbundanceSummary/3-Barplots/taxa-bar-plots_Qiime2/taxa-bar-plots.qzv Figure3-2.qzv
	#cp ../2-AbundanceAnalysis/1-AbundanceSummary/4-Heatmaps/${category_1}-table-Phylum.${number}/${category_1}-table-Phylum.${number}.qzv Figure3-3.qzv
	#cp ../2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/${category_1}.ANCOM.Genus/${category_1}.ANCOM.Genus.qzv Figure3-4.qzv
	#cp ../3-AlphaDiversity/2-AlphaRarefaction/alpha-rarefaction/alpha-rarefaction.qzv Figure3-6.qzv
	#cp ../3-AlphaDiversity/3-SignificanceAnalysis/2-Kruskal_Wallis/shannon-group-significance/shannon-group-significance.qzv Figure3-7.qzv
	
	#cp ../4-BetaDiversity/2-PCoA/unweighted_unifrac_emperor/unweighted_unifrac_emperor.qzv Figure3-9.qzv
	#cp ../4-BetaDiversity/3-NMDS/${category_1}_unifrac_NMDS.pdf Figure3-10.pdf
	#cp ../4-BetaDiversity/5-GroupSignificance/unweighted_unifrac-permanova-${category_1}-significance/unweighted_unifrac-permanova-${category_1}-significance.qzv Figure3-11.qzv
	#cp ../5-Phylogenetics/1-MajorPhylums/Bacteria.phylogeny.pdf Figure3-12.pdf
	#cp ../6-AssociationAnalysis/1-RDA/Genus/*_bacteria_location_plot.png Figure3-14.png
	#cp ../6-AssociationAnalysis/permanova.pdf Figure3-15.pdf
	#cp ../7-FunctionAnalysis/1-KEGG_Pathway/percent.feature-table.metagenome.L1.png Figure3-19.png
	#cp ../7-FunctionAnalysis/2-PCAPlots/feature-table.metagenome.L1.${category_1}.PCA.pdf Figure3-20.pdf
	#cp ../7-FunctionAnalysis/3-TreeBasedPlots/tree.feature-table.metagenome.L1.png Figure3-21.png

}

MAIN() {

	echo "##############################################################\n#Organize the Result folder"
	organize_deliverable_structure
<<COMMENT1
	echo "##############################################################\n#Organize the Essential folder ----- part1"
	mkdir Essential
	mkdir Essential/1-Demux/ Essential/2-AbundanceAnalysis/ Essential/3-AlphaDiversity Essential/4-BetaDiversity
	mkdir Essential/2-AbundanceAnalysis/OTUSummary Essential/2-AbundanceAnalysis/OTUDifferentialAnalysis Essential/2-AbundanceAnalysis/OTUSummary/Heatmap Essential/2-AbundanceAnalysis/OTUDifferentialAnalysis/ANCOM
	cp -r demux.qzv.exported Essential/1-Demux/
	cp -r taxa-bar-plots.qzv.exported Essential/2-AbundanceAnalysis/OTUSummary/
	cp exported/feature-table.taxonomy.txt Essential/2-AbundanceAnalysis/OTUSummary/
	cp -r exported/Relative Essential/2-AbundanceAnalysis/OTUSummary/
	cp -r exported/exported/${number}/table-l2.${number} exported/exported/${number}/table-l6.${number} Essential/2-OTUTable/Heatmap
	cp -r exported/ANCOM/ANCOM.l2.qzv.exported Essential/2-AbundanceAnalysis/OTUDifferentialAnalysis/ANCOM/phylum
	cp -r exported/ANCOM/ANCOM.l6.qzv.exported Essential/2-AbundanceAnalysis/OTUDifferentialAnalysis/ANCOM/genus
	cp exported/kruskal_wallis* exported/ANOVA* Essential/2-AbundanceAnalysis/OTUDifferentialAnalysis
	cp alpha/alpha-summary.tsv Essential/3-AlphaDiversity
	cp R_output/bray_matrix.txt R_output/wunifrac_matrix.txt R_output/unifrac_matrix.txt R_output/*unifrac_NMDS.pdf R_output/bray_NMDS.pdf R_output/P*_plot.pdf Essential/4-BetaDiversity/

	cp -r Essential Result
COMMENT1
}

MAIN;
