from .libShell import qiime2, qiime1, lefse
from . import settings
from .systemMixin import SystemMixin
from .utils import generate_tree, soar_outpath, picrust2_pipline
import os
import shutil
import pandas as pd
import logging
logger = logging.getLogger(__name__)


class VisualizeAmplicon(SystemMixin):

    def __init__(self, otu_table, rep_seqs, mapping_file,
                 category, out_dir="./", taxonomy=None, exclude="none",
                 unspecified=True, asv_type="16s", processors=1, report=True,
                 colors=False, min_align=False, deseq2=True, orders=[]):
        asv_type = asv_type.lower()
        arguments = locals()
        if taxonomy:
            self.set_path(force=False, taxonomy=taxonomy)
        else:
            self.set_attr(taxonomy="")
        self.set_path(
            force=False,
            # base_dir=os.path.dirname(__file__),
            otu_table=otu_table,
            rep_seqs=rep_seqs,
            mapping_file=mapping_file
        )
        self.load_settings_path()
        self.map = pd.read_csv(mapping_file, sep="\t", index_col=0)
        groups = self.map[category]
        groups = groups.loc[groups.notna()]
        self.groups = list(set(groups))
        self.groups.sort()
        abc_dir = os.path.join(out_dir, "2-AbundanceAnalysis")
        abc_sum = os.path.join(abc_dir, "1-AbundanceSummary")
        abc_com = os.path.join(abc_dir, "2-AbundanceComparison")
        abc_tab = os.path.join(abc_sum, "1-AbundanceTable")
        alpha_dir = os.path.join(out_dir, "3-AlphaDiversity")
        beta_dir = os.path.join(out_dir, "4-BetaDiversity")
        cor_dir = os.path.join(out_dir, "6-AssociationAnalysis")
        tmp_dir = os.path.join(out_dir, "_tmp_dir")

        self.set_path(
            force=True,
            out_dir=out_dir,
            tmp_dir=tmp_dir,
            abc_rel=os.path.join(tmp_dir, "relative"),
            venn_dir=os.path.join(out_dir, "1-VennAndFlower"),
            tree_dir=os.path.join(out_dir, "5-Phylogenetics"),
            abc_dir=abc_dir,
            abc_sum=abc_sum,
            abc_com=abc_com,
            abc_tab=abc_tab,
            abc_abs=os.path.join(abc_tab, "1-Absolute"),
            abc_rel_unspecified=os.path.join(abc_tab, "2-Relative"),
            abc_even=os.path.join(abc_com, "tables_for_deseq_anova_kruskal"),
            abc_bar=os.path.join(abc_sum, "2-Barplots"),
            abc_heatmap=os.path.join(abc_sum, "3-Heatmaps"),
            alpha_dir=alpha_dir,
            alpha_sum=os.path.join(alpha_dir, "1-AlphaDiversitySummary"),
            alpha_rare=os.path.join(alpha_dir, "2-AlphaRarefaction"),
            alpha_sig=os.path.join(alpha_dir, "3-SignificanceAnalysis"),
            beta_dir=beta_dir,
            beta_pcoa=os.path.join(beta_dir, "2-PCoA"),
            beta_sig=os.path.join(beta_dir, "5-GroupSignificance"),
            cor_dir=cor_dir,
            cor_rda=os.path.join(cor_dir, "1-RDA"),
            cor_heatmap=os.path.join(cor_dir, "2-CorrelationHeatmap"),
            cor_network=os.path.join(cor_dir, "3-NetworkAnalysis"),
            func_dir=os.path.join(out_dir, "7-FunctionAnalysis")
        )
        self.levels_map = {
            "k": "Kingdom",
            "p": "Phylum",
            "c": "Class",
            "o": "Order",
            "f": "Family",
            "g": "Genus",
            "s": "Species"
        }
        self.levels = self.levels_map.values()
        self.set_attr(
            category=category,
            exclude=exclude,
            asv_type=asv_type,
            processors=processors,
            report=report,
            colors=colors,
            min_align=min_align,
            deseq2=deseq2
        )
        self.unspecified = unspecified
        self.__tables = {}
        self.table_patterns = dict(
            rel="{abc_rel}/otu_table.{level}.relative.xls",
            rel_unspecified="{abc_rel_unspecified}/otu_table.{level}.relative.xls",
            abs="{abc_abs}/otu_table.{level}.absolute.xls",
            biom="{abc_even}/otu_table.{level}.even.filter.biom"
        )
        self.exclude_list = []
        exclude_list = [] if self.exclude == "all" else self.exclude.split(';')
        self.orders = orders
        for el in exclude_list:
            el, *aprefix = el.split(":") + [""]
            prefix = "".join(aprefix)
            self.exclude_list.append((el, prefix))

    def load_settings_path(self):
        self.set_path(force=False, **settings.path)

    def tax_levels(self, strip=True):
        levels = []
        iter_level = self.levels_map.items()
        for i, ol in enumerate(iter_level):
            level_no = i + 1
            level_abbr, level = ol
            if (level_no == 1 or level_no > 7) and strip:
                continue
            levels.append((level_no, level_abbr, level))
        return levels

    @property
    def all_levels(self):
        return self.tax_levels(strip=False)

    @property
    def main_levels(self):
        return self.tax_levels(strip=True)

    def prepare_dirs(self, dirs):
        for d in dirs:
            d = d.format(**self.context)
            if os.path.exists(d):
                self.system("rm -r {d}", d=d)
            os.makedirs(d)

    @soar_outpath(
        otu_table="{tmp_dir}/table.qza",
        rep_seqs="{tmp_dir}/rep-seqs.qza",
        taxonomy="{tmp_dir}/taxonomy.qza",
        otu_abs_tsv="{abc_dir}/feature-table.taxonomy.xls",
        consensus="{tmp_dir}/consensus.txt",
        mapping_file="{out_dir}/mapping_file.xls",
        depth_file="{tmp_dir}/depth.txt"
    )
    def process_source(self):
        with qiime2():
            self.system("""
{bayegy_home}/filter_source.py -i {otu_table} -r {rep_seqs} \
 -t '{taxonomy}' -o {tmp_dir} -m {mapping_file} -c {category} \
 --tsv {abc_dir}/feature-table.taxonomy.xls \
 --fmap {out_dir}/mapping_file.xls \
 --odepth {tmp_dir}/depth.txt \
 --consensus {tmp_dir}/consensus.txt""")

    def set_depth(self):
        with open(self.depth_file) as f:
            self.set_attr(depth=int(f.read().strip()))

    @soar_outpath(
        evenconsensus="{tmp_dir}/evenconsensus.txt",
        classified_stat="{tmp_dir}/classified_stat_relative.xls"
    )
    def calculate_otu_tables(self):
        stat_script = "stat_otu_tab.unspecifiedadded.pl" if self.unspecified else "stat_otu_tab.pl"
        min_sam = int(self.map.shape[0] / 4)
        with qiime1():
            self.system("""
{perl_path} {bayegy_home}/{stat_script} -unif min {otu_abs_tsv} -prefix {abc_rel_unspecified}/otu_table \
 --even {abc_dir}/feature-table.taxonomy.even.xls -spestat {tmp_dir}/classified_stat_relative.xls > /dev/null && \
biom convert -i {abc_dir}/feature-table.taxonomy.even.xls  -o {tmp_dir}/otu_table.even.biom --to-hdf5 \
 --table-type="OTU table" --process-obs-metadata taxonomy && \
summarize_taxa.py -i {tmp_dir}/otu_table.even.biom -a -L 1,2,3,4,5,6,7 -o {abc_even};
{perl_path} {bayegy_home}/{stat_script} -unif min {otu_abs_tsv} --prefix {abc_abs}/otu_table -nomat -abs  > /dev/null;
{perl_path} {bayegy_home}/stat_otu_tab.pl -unif min {otu_abs_tsv} -prefix {abc_rel}/otu_table \
 -spestat {tmp_dir}/classified_stat_relative.xls;""", **locals())
            for level_no, level_abbr, level in self.all_levels:
                self.system("""
filter_otus_from_otu_table.py -i {abc_even}/otu_table.even_L{level_no}.biom -s {min_sam} \
 -o {abc_even}/otu_table.{level}.even.filter.biom;
mv {abc_even}/otu_table.even_L{level_no}.txt {abc_even}/otu_table.{level}.even.xls;
mv {abc_abs}/otu_table.{level_abbr}.absolute.mat {abc_abs}/otu_table.{level}.absolute.xls;
mv {abc_rel}/otu_table.{level_abbr}.relative.mat {abc_rel}/otu_table.{level}.relative.xls;
mv {abc_rel_unspecified}/otu_table.{level_abbr}.relative.mat \
 {abc_rel_unspecified}/otu_table.{level}.relative.xls""", **locals())
        self.system("""
sed 's/taxonomy/Consensus Lineage/' < {abc_dir}/feature-table.taxonomy.even.xls > {tmp_dir}/evenconsensus.txt;
rm {abc_tab}/*/*.mat""")

    def taxa_levels_stat(self):
        self.system("""
{R_path} {bayegy_home}/collapse_table_with_group_mean.R -m {mapping_file} -c {category} \
 -t {classified_stat} -o {abc_dir}/ && \
{perl_path} {bayegy_home}/bar_diagram.pl -table {abc_dir}/{category}_classified_stat_relative.xls -style 1 \
 -x_title "" -y_title "Sequence Number Percent" -right -textup -rotate='-45' \
 --y_mun 1,7 > {tmp_dir}/{category}_classified_stat_relative.svg && \
convert -density 300 -quality 100 {tmp_dir}/{category}_classified_stat_relative.svg \
 {abc_dir}/{category}_classified_stat_relative.png;
{R_path} {bayegy_home}/taxonomy_stats.R -i {consensus} -o {abc_dir}""")

    def tables(self, tt):
        if not self.__tables.get(tt):
            tbs = []
            self.__tables[tt] = tbs
            ptn = self.table_patterns[tt]
            for level_no, level_abbr, level in self.main_levels:
                tb = ptn.format(**self.context, **locals())
                if os.path.exists(tb):
                    tbs.append((level, tb))
        return self.__tables[tt]

    def venn_and_flower(self):
        self.system("""
{R_path} {bayegy_home}/venn_and_flower_plot.R  -i {otu_abs_tsv} --skip F \
 -m {mapping_file} -c {category} -o {venn_dir}""")

    @soar_outpath(
        rooted_tree="{tmp_dir}/rooted-tree.qza",
        # figure_tree="{tree_dir}/{category}_phylogenetic_tree_heatmap.pdf"
        tree_nwk="{tree_dir}/tree.rooted.nwk"
    )
    def tree(self):
        generate_tree(
            rep_seqs_qza=self.rep_seqs,
            nwk_out=self.tree_dir,
            qza_out=self.tmp_dir
        )
        self.system("""
{R_path} {bayegy_home}/phylotree_and_heatmap.R -i {consensus} -t {tree_dir}/tree.rooted.nwk \
 -m {mapping_file} -c {category} -o {tree_dir} -n 50""")

    def qiime_bar(self):
        with qiime2():
            self.system("""
qiime taxa barplot --i-table {otu_table} --i-taxonomy {taxonomy} \
 --m-metadata-file {mapping_file}  --o-visualization {abc_bar}/taxa-bar-plots_Qiime2.qzv""")

    def run_deseq2(self, table, level):
        gp_num = len(self.groups)
        for i in range(gp_num):
            for j in range(gp_num):
                if i < j:
                    group1, group2 = self.groups[i], self.groups[j]
                    self.system("""
differential_abundance.py -i {table} \
 -o {abc_com}/5-DESeq2/DESeq2_{category}_Between_{group1}_and_{group2}_DiffAbundance_{level}.txt \
 -a DESeq2_nbinom -m {mapping_file} -c {category} -x '{group1}' -y '{group2}' -d""", **locals())

    def qiime1_levels_analyze(self):
        dirs = ['{abc_com}/3-KruskalWallis', '{abc_com}/2-ANOVA']
        if self.deseq2:
            dirs.append('{abc_com}/5-DESeq2')
        self.prepare_dirs(dirs)
        with qiime1():
            for level, tb in self.tables('biom'):
                self.system("""
group_significance.py -i {tb} -m {mapping_file} -c {category} -s kruskal_wallis \
 -o {abc_com}/3-KruskalWallis/kruskal_wallis_{category}_DiffAbundance_{level}.txt \
 --biom_samples_are_superset --print_non_overlap && \
group_significance.py -i {tb} -m {mapping_file} -c {category} -s ANOVA \
 -o {abc_com}/2-ANOVA/ANOVA_{category}_DiffAbundance_{level}.txt \
 --biom_samples_are_superset --print_non_overlap""", **locals())
                if self.deseq2:
                    self.run_deseq2(tb, level)

    def qiime2_levels_analyze(self):
        self.prepare_dirs(["{tmp_dir}/collapsed", "{abc_com}/1-ANCOM", "{abc_tab}/3-CollapsedStats"])
        with qiime2():
            for level_no, level_abbr, level in self.main_levels:
                self.system("""
qiime taxa collapse  --i-table {otu_table} --i-taxonomy {taxonomy} \
 --p-level {level_no} --o-collapsed-table {tmp_dir}/collapsed/collapsed-{level}.qza && \
qiime feature-table summarize --i-table {tmp_dir}/collapsed/collapsed-{level}.qza \
 --m-sample-metadata-file {mapping_file}\
 --o-visualization {abc_tab}/3-CollapsedStats/collapsed-{level}.qzv;
qiime composition add-pseudocount   --i-table {tmp_dir}/collapsed/collapsed-{level}.qza \
 --o-composition-table {tmp_dir}/collapsed/composition.{level}.qza&&\
qiime composition ancom  --i-table {tmp_dir}/collapsed/composition.{level}.qza \
 --m-metadata-file {mapping_file} --m-metadata-column {category} \
 --o-visualization {abc_com}/1-ANCOM/{category}.ANCOM.{level}.qzv""", **locals())

    @soar_outpath(alpha_table="{alpha_sum}/alpha-summary.xls")
    def qiime_alpha(self,
                    alpha_indexes=["chao1", "shannon", "observed_otus", "faith_pd", "simpson"],
                    steps=50):
        self.prepare_dirs(["{tmp_dir}/alpha", "{alpha_sig}/2-Kruskal_Wallis"])
        with qiime2():
            max_steps = self.depth - 1
            self.system("""
qiime diversity alpha-rarefaction --i-table {otu_table} --i-phylogeny {rooted_tree} \
 --p-max-depth {depth} --m-metadata-file {mapping_file} --p-steps {steps} \
 --o-visualization {alpha_rare}/alpha-rarefaction-Qiime2.qzv""",
                        steps=steps if steps < max_steps else max_steps)
            for alpha_index in alpha_indexes:
                if alpha_index == "faith_pd":
                    self.system("""
qiime diversity alpha-phylogenetic --i-table {otu_table} --i-phylogeny {rooted_tree} \
 --p-metric faith_pd --output-dir {tmp_dir}/alpha/faith_pd""")
                else:
                    self.system("""
qiime diversity alpha --i-table {otu_table} --p-metric {alpha_index} \
 --output-dir {tmp_dir}/alpha/{alpha_index}/""", alpha_index=alpha_index)

                self.system("""
qiime tools export --input-path {tmp_dir}/alpha/{alpha_index}/alpha_diversity.qza \
 --output-path {tmp_dir}/alpha/{alpha_index}/&& \
qiime diversity alpha-group-significance  \
 --i-alpha-diversity {tmp_dir}/alpha/{alpha_index}/alpha_diversity.qza  \
 --m-metadata-file {mapping_file} \
 --o-visualization {alpha_sig}/2-Kruskal_Wallis/{alpha_index}-group-significance.qzv""", alpha_index=alpha_index)
        self.system("""{python3_path} {bayegy_home}/merge_tables.py {tmp_dir}/alpha/*/*.tsv {alpha_sum}/alpha-summary.xls""")

    def qiime_beta(self, methods=["permanova", "anosim"]):
        with qiime2():
            beta_tmp = "{tmp_dir}/core-metrics-results".format(**self.context)
            if os.path.exists(beta_tmp):
                shutil.rmtree(beta_tmp)
            self.prepare_dirs(["{beta_pcoa}/PCoA-Qiime2"])
            self.system("""
qiime diversity core-metrics-phylogenetic --i-phylogeny {rooted_tree} \
 --i-table {otu_table}   --p-sampling-depth {depth}  \
 --m-metadata-file {mapping_file}  --output-dir {beta_tmp}&&\
mv {beta_tmp}/*emperor.qzv {beta_pcoa}/PCoA-Qiime2/&&\
rm {beta_pcoa}/PCoA-Qiime2/jaccard_emperor.qzv""", **locals())
            for beta_index in ["unweighted_unifrac", "weighted_unifrac", "bray_curtis"]:
                for method in methods:
                    self.system("""
qiime diversity beta-group-significance \
 --i-distance-matrix {beta_tmp}/{beta_index}_distance_matrix.qza  \
 --m-metadata-file {mapping_file}  --p-method {method} --m-metadata-column {category}  \
 --o-visualization '{beta_sig}/{beta_index}-{method}-{category}-significance.qzv' \
 --p-pairwise""", **locals())

    def export_qzv(self):
        with qiime2():
            self.system("""
for f in $(find {out_dir} -type f -name "*.qzv");
    do echo $f && \
    base=$(basename $f .qzv) && \
    dir=$(dirname $f) && \
    new=${{dir}}/${{base}} && \
    qiime tools export --input-path $f --output-path ${{new}}/;
done;
for f in $(find {out_dir} -type f -name "index.html");
    do echo $f; dir=$(dirname $f);
    mv $f ${{dir}}/Summary_请点此文件查看.html;
done;""")

    @soar_outpath(
        rep_seqs_fasta="{tree_dir}/representative-sequences.fasta"
    )
    def export_qza(self):
        with qiime2():
            self.system("""
qiime tools export --input-path {rep_seqs} --output-path {tree_dir} &&
mv {tree_dir}/dna-sequences.fasta {tree_dir}/representative-sequences.fasta""")

    def base_alpha_beta(self):
        self.set_path(
            force=False,
            qiime_alpha_res="{alpha_rare}/alpha-rarefaction-Qiime2"
        )
        self.system("""
{R_path} {bayegy_home}/base_alpha_beta.R -i {evenconsensus} -a {alpha_table} -t {tree_nwk} \
 -m {mapping_file} -c {category} --output-pcoa {beta_pcoa}/PCoA-Phyloseq \
 --output-nmds {beta_dir}/3-NMDS --output-plsda {beta_dir}/4-PLS-DA \
 --output-pca {beta_dir}/6-PCA  --output-alpha-heatmap {alpha_sum};
{R_path} {bayegy_home}/beta_heatmap.R -i {evenconsensus} -m {mapping_file} -t {tree_nwk} \
 -o {beta_dir}/1-BetaDiversitySummary -c {category} -p 'unclustered_';
{R_path} {bayegy_home}/beta_heatmap.R -i {evenconsensus}  -t {tree_nwk} \
 -o {beta_dir}/1-BetaDiversitySummary;
{perl_path} {bayegy_home}/table_data_svg.pl --colors cyan-orange {beta_dir}/1-BetaDiversitySummary/*matrix.xls \
 --symbol 'Beta Diversity' > {tmp_dir}/BetaDiversity_heatmap.svg;
rsvg-convert -h 3200 -b white {tmp_dir}/BetaDiversity_heatmap.svg \
 > {beta_dir}/1-BetaDiversitySummary/BetaDiversity_heatmap.png;
{R_path} {bayegy_home}/alphararefaction.R -i {qiime_alpha_res} -o {alpha_rare}/alpha-rarefaction-ggplot2;
{R_path} {bayegy_home}/alphararefaction.R -i {qiime_alpha_res} -o {alpha_rare}/alpha-rarefaction-ggplot2 \
 -m {mapping_file} -c {category} -p "group_mean_";
{R_path} {bayegy_home}/alphaboxplotwitSig.R -m {mapping_file} -c {category} \
 -i {alpha_table} -o {alpha_sig}/1-Wilcox_Test""")

    def levels_analyze(self):
        for level, tb in self.tables('abs'):
            self.system("""
{R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {tb} \
 -o {abc_heatmap}/Heatmap_top20/{level}/ -l T -t F;
{R_path} {bayegy_home}/abundance_heatmap.R -m {mapping_file} -c {category}  -n 20 -i {tb} \
 -o {abc_heatmap}/Heatmap_top20_clustered/{level}/ -l T -t F -u T;
{R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {tb} \
 -o {abc_heatmap}/Heatmap_top20/{level}/ -b T -l T -p 'group_mean_' -t T""", **locals())
            for order in self.orders:
                self.system("""
{R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {tb} \
 -o {abc_heatmap}/Heatmap_top20/{level}/ -l T -t F -O {order} -p in_{order}_;
{R_path} {bayegy_home}/abundance_heatmap.R  -m {mapping_file} -c {category} -n 20 -i {tb} \
 -o {abc_heatmap}/Heatmap_top20/{level}/ -b T -l T -p 'group_mean_in_{order}_' \
 -t T -O {order};""", **locals())

        for level, tb in self.tables('rel_unspecified'):
            self.system("""
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {tb} \
 -o {abc_bar}/taxa-bar-plots-top20 -p {level}_{category}_ -b F;
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {tb} \
 -o {abc_bar}/taxa-bar-plots-top20-group-mean -p {category}_{level}_mean_ -b T""", **locals())
            for order in self.orders:
                self.system("""
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {tb} \
 -o {abc_bar}/taxa-bar-plots-top20 -p {level}_in_{order}_ -b F -O {order};
 {R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {tb} \
 -o {abc_bar}/taxa-bar-plots-top20-group-mean -p {level}_in_{order}_ \
 -b T -O {order};""", **locals())
            for el, prefix in self.exclude_list:
                self.system("""
{R_path} {bayegy_home}/RDA.R -i {tb} -m {mapping_file} -c {category} \
 -o {cor_rda}/{level} -n 25 -e {el} -p '{prefix}'""", **locals())
        for level, tb in self.tables('rel'):
            for el, prefix in self.exclude_list:
                self.system("""
{R_path} {bayegy_home}/cor_heatmap.R -i {tb} -o {cor_heatmap}/{level}/ -n 25 \
 -m {mapping_file} -e {el} -p '{prefix}'""", **locals())
            self.system("""
{R_path} {bayegy_home}/network.R -c 0.5 -i {tb} -o {cor_network}/{level}/;
{R_path} {bayegy_home}/write_data_for_lefse.R -i  {tb} -m  {mapping_file} -c  {category} \
 -o  {abc_com}/4-LEfSe/{level}/{category}_{level}_lefse.txt -u l""", **locals())
            with lefse():
                self.system("""
base={abc_com}/4-LEfSe/{level}/{category}_{level}_lefse_LDA2 && \
{bayegy_home}/mod_format_input.py {abc_com}/4-LEfSe/{level}/{category}_{level}_lefse.txt \
 ${{base}}.lefseinput.txt -c 2 -u 1 -o 1000000 && \
{bayegy_home}/mod_run_lefse.py ${{base}}.lefseinput.txt ${{base}}.LDA.txt -l 2 && \
{bayegy_home}/mod_lefse-plot_res.py --category {category} --map {mapping_file} \
 --max_feature_len 200 --orientation h --format pdf ${{base}}.LDA.txt ${{base}}.pdf;
{bayegy_home}/mod_lefse-plot_cladogram.py ${{base}}.LDA.txt --map {mapping_file} \
 --category {category} ${{base}}.cladogram.pdf --format pdf;
base4={abc_com}/4-LEfSe/{level}/{category}_{level}_lefse_LDA4 && \
{metagenome_home}/lda22ldamt.py ${{base}}.LDA.txt ${{base4}}.LDA.txt 4 && \
{bayegy_home}/mod_lefse-plot_res.py --category {category} --map {mapping_file} \
 --max_feature_len 200 --orientation h --format pdf ${{base4}}.LDA.txt ${{base4}}.pdf;
{bayegy_home}/mod_lefse-plot_cladogram.py ${{base4}}.LDA.txt --map {mapping_file} \
 --category {category} ${{base4}}.cladogram.pdf --format pdf;""", **locals())

    def function(self):
        if self.asv_type not in ["its", "16s", "18s"]:
            print("{} ASV type is currently not supported for function analysis.".format(self.asv_type))
            return
        picrust2_pipline(self.otu_abs_tsv, self.rep_seqs_fasta,
                         out_dir=self.func_dir, asv_type=self.asv_type,
                         processors=self.processors, min_align=self.min_align)
        if self.asv_type == "16s":
            kegg_pathways = [
                ("KEGG.L%s_" % (ln), "{func_dir}/KEGG.pathway.L%s.abundance.xls" % (ln)) for ln in (1, 2, 3)]
        else:
            kegg_pathways = []
        for prefix, func_abc_file in [
            ("MetaCyc_", "{func_dir}/MetaCyc.pathway.abundance.xls"),
            *kegg_pathways
        ]:
            func_abc_file = func_abc_file.format(**self.context)
            if not os.path.exists(func_abc_file):
                logger.warning(
                    "{func_abc_file} does not exists! function prediction may get an empty result".format(**locals())
                )
                return
            self.system("""
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {func_abc_file} \
 -o {func_dir}/1-Barplots/bar-plots-top20/ -p {prefix} -b F -s F;
{R_path} {bayegy_home}/abundance_barplot.R -n 20 -m {mapping_file} -c {category} -i {func_abc_file} \
 -o {func_dir}/1-Barplots/bar-plots-top20-group-mean/ -p {prefix} -b T -s F;
{R_path} {bayegy_home}/PCA.R --input {func_abc_file} --map {mapping_file} --category {category} \
 -o {func_dir}/2-PCAPlots -p {prefix};
{R_path} {bayegy_home}/Function_DunnTest.r -i {func_abc_file}  -m {mapping_file} \
 -p {prefix} -c {category} -o {func_dir}/3-SignifcanceAnalysis/1-DunnTest > /dev/null;
{R_path} {bayegy_home}/function_barplot.R -i {func_abc_file} -m {mapping_file} -c {category} --skip F \
 -j T -a 0.05 -b F --feature-col 1 -o {func_dir}/3-SignifcanceAnalysis/2-ANOVA_And_Duncan -p {prefix}""", **locals())
            if func_abc_file.endswith("KEGG.pathway.L3.abundance.xls"):
                self.system("""
{R_path} {bayegy_home}/function_barplot.R -i {func_abc_file} \
 -m {mapping_file} -c {category} --skip F -j T -a 0.05 -b T \
 --feature-col -1 -o {func_dir}/3-SignifcanceAnalysis/2-ANOVA_And_Duncan""", **locals())

    def generate_report(self):
        self.set_path(force=True, report_dir=os.path.join(self.out_dir, "FiguresTablesForReport"))
        self.system("""
cp -rp {bayegy_home}/Report/src {out_dir}/FiguresTablesForReport/
cp {bayegy_home}/Report/结题报告.html {out_dir}/
cp {out_dir}/2-AbundanceAnalysis/1-AbundanceSummary/3-Heatmaps/Heatmap_top20_clustered/Phylum/heatmap.pdf {report_dir}/Figure4-3.pdf
cp -rp {out_dir}/2-AbundanceAnalysis/2-AbundanceComparison/1-ANCOM/{category}.ANCOM.Genus/ {report_dir}/page4-2
cp -rp {out_dir}/4-BetaDiversity/5-GroupSignificance/unweighted_unifrac-permanova-{category}-significance/ {report_dir}/page6-2
cp {out_dir}/2-AbundanceAnalysis/1-AbundanceSummary/2-Barplots/taxa-bar-plots-top20-group-mean/{category}_Phylum_mean_barplot.pdf {report_dir}/Figure4-2.pdf
cp {out_dir}/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/Genus/{category}_Genus_lefse_LDA2.pdf {report_dir}/Figure4-4.pdf
cp {out_dir}/2-AbundanceAnalysis/2-AbundanceComparison/4-LEfSe/Genus/{category}_Genus_lefse_LDA2.cladogram.pdf {report_dir}/Figure4-5.pdf
cp {out_dir}/1-VennAndFlower/{category}_Venn_plot.png {report_dir}/Figure4-6.png
cp {out_dir}/3-AlphaDiversity/1-AlphaDiversitySummary/{category}_alpha_diversity_shannon.boxplot.pdf {report_dir}/Figure5-1.pdf
cp {out_dir}/3-AlphaDiversity/3-SignificanceAnalysis/1-Wilcox_Test/shannon_{category}_wilcox_compare_boxplot.png {report_dir}/Figure5-2.png
cp {out_dir}/4-BetaDiversity/1-BetaDiversitySummary/BetaDiversity_heatmap.png {report_dir}/Figure6-1.png
cp {out_dir}/4-BetaDiversity/3-NMDS/{category}_unweighted_unifrac_NMDS_without_labels.pdf {report_dir}/Figure6-2.pdf
cp {out_dir}/5-Phylogenetics/{category}_phylogenetic_tree_heatmap.pdf {report_dir}/Figure7-1.pdf
{python3_path} {metagenome_home}/cpfirst {out_dir}/6-AssociationAnalysis/1-RDA/Genus/features_location_plot.pdf {report_dir}/Figure8-1.pdf
{python3_path} {metagenome_home}/cpfirst {out_dir}/6-AssociationAnalysis/2-CorrelationHeatmap/Genus/Correlation_heatmap.pdf {report_dir}/Figure8-2.pdf
cp {out_dir}/6-AssociationAnalysis/3-NetworkAnalysis/Genus/Correlation_network.pdf {report_dir}/Figure8-3.pdf
cp {out_dir}/7-FunctionAnalysis/1-Barplots/bar-plots-top20-group-mean/MetaCyc_barplot.pdf {report_dir}/Figure9-1.pdf
cp {out_dir}/7-FunctionAnalysis/2-PCAPlots/MetaCyc_PCA.pdf {report_dir}/Figure9-2.pdf
cp {out_dir}/7-FunctionAnalysis/3-SignifcanceAnalysis/2-ANOVA_And_Duncan/MetaCyc_{category}_all_significant_pathway_barplot_of_duncan.pdf {report_dir}/Figure9-3.pdf
for pdf in {report_dir}/*.pdf; do echo $pdf; dir=$(dirname $pdf); base=$(basename $pdf .pdf); convert  -density 300 -quality 80 $pdf ${{dir}}/${{base}}.png; rm $pdf;done;
{python3_path} {bayegy_home}/convert_to_html_table.py -i {out_dir}/../1-OTUStats/2-Stats-OTU/Summary_请点此文件查看.html \
 -o {report_dir}/src/pages/main_cleaned.html -t dada2html -k '无法加载表格表3-1'""")

    def clean_up(self):
        self.system("""
for f in $(find {out_dir} -type f -name "*.qzv");
    do rm $f;
done;
rm -r {tmp_dir};
if [ -d {abc_com}/5-DESeq2 ]; then rm {abc_com}/5-DESeq2/*pdf; fi;
{python3_path} {metagenome_home}/change_suffix.py {out_dir} -s '1-Stats-demux,2-Stats-dada2,3-RepresentiveSequence,taxa-bar-plots_Qiime2,\
1-ANCOM,alpha-rarefaction-Qiime2,2-Kruskal_Wallis,PCoA-Qiime2,5-GroupSignificance,FiguresTablesForReport'""")

    def set_colors(self, colors=False):
        if colors:
            colors_list = '\n'.join([colors[k] for k in sorted(colors.keys(), key=str.lower)])
            colors_list_file = os.path.join(self.out_dir, 'group_color.list')
            with open(colors_list_file, 'w') as f:
                f.write(colors_list)
        else:
            colors_list_file = "{bayegy_home}/piputils/group_color.list".format(**self.context)
        self.system(
            "{python3_path} {bayegy_home}/piputils/write_colors_plan.py -i {mapping_file} -c {category} \
            -p {colors_list_file} -o {out_dir}/colors_plan.json", colors_list_file=colors_list_file)
        os.environ['COLORS_PLAN_PATH'] = os.path.join(self.out_dir, 'colors_plan.json')

    def visualize(self):
        if self.colors:
            self.set_colors(self.colors)
        self.process_source()  # not skippable
        self.calculate_otu_tables()  # not skippable
        self.taxa_levels_stat()
        self.qiime1_levels_analyze()
        self.qiime2_levels_analyze()
        self.venn_and_flower()
        self.tree()  # not skippable
        self.set_depth()  # not skippable
        self.qiime_alpha()  # not skippable
        self.qiime_beta()
        self.qiime_bar()
        # self.qiime2_export()
        self.export_qzv()  # not skippable
        self.export_qza()  # not skippable
        self.base_alpha_beta()
        self.levels_analyze()
        self.function()
        if self.report:
            self.generate_report()
        self.clean_up()
