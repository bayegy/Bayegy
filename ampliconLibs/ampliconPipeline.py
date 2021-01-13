from .systemMixin import SystemMixin
import os
from . import settings
from .libShell import qiime2
from .commandArgs import CommandArgs
from .visualizeAmplicon import VisualizeAmplicon
import re
import json
from .utils import generate_tree, soar_outpath
import logging
logger = logging.Logger(__name__)

class AmpliconPipeline(SystemMixin):

    def __init__(self, raw_dir, pre_map,
                 out_dir="./",
                 sam_pattern=r'raw\.split\.(.+)\.[12]\.fq$',
                 fwd_pattern=r'\.1\.fq$',
                 rev_pattern=r'\.2\.fq$',
                 fwd_primer="",  # 338F ACTCCTACGGGAGGCAGCAG
                 rev_primer="",  # 806R GGACTACHVGGGTWTCTAAT
                 raw_type="16S_PE",
                 db_classifier="",
                 db_taxonomy="",
                 db_rep_set="",
                 threads=16,
                 otu_method="dada2",
                 classify_method="sklearn",
                 classify_args="",
                 taxa_filter="exclude:mitochondria,chloroplast,Unassigned",
                 otu_args="",
                 join_args="",
                 categories=False,
                 orders=False,
                 **visualize_kwargs
                 ):
        asv_type, end_type = raw_type.lower().split("_")
        db_taxonomy = db_taxonomy or settings.db_taxonomy.get(asv_type, "")
        db_rep_set = db_rep_set or settings.db_rep_set.get(asv_type, "")
        if not ((db_taxonomy and db_rep_set) or db_classifier):
            raise Exception("Database is not specified!")
        arguments = locals()
        dy_paths = {p: arguments[p] for p in ["db_classifier", "db_taxonomy", "db_rep_set"] if arguments[p]}
        self.set_path(
            force=False,
            # base_dir=os.path.dirname(__file__),
            raw_dir=raw_dir,
            pre_map=pre_map,
            **dy_paths,
            **settings.path
        )
        result_dir = os.path.join(out_dir, "result")
        data_dir = os.path.join(out_dir, "_DATA_DIR")
        tmp_dir = os.path.join(out_dir, "_TMP_DIR")
        stat_dir = os.path.join(result_dir, "1-OTUStats")
        self.set_path(
            force=True,
            out_dir=out_dir,
            tmp_dir=tmp_dir,
            data_dir=data_dir,
            result_dir=result_dir,
            stat_dir=stat_dir,
            demux_stat_dir=os.path.join(stat_dir, "1-Stats-demux"),
            otu_stat_dir=os.path.join(stat_dir, "2-Stats-OTU"),
            rep_seq_dir=os.path.join(stat_dir, "3-RepresentiveSequence/rep-seqs"),
            classifier_dir=os.path.join(data_dir, "classifier")
        )
        self.set_attr(
            sam_pattern=sam_pattern,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern if end_type == "pe" else 'NONE',
            fwd_primer=fwd_primer,
            rev_primer=rev_primer,
            threads=threads,
            end_type=end_type,
            asv_type=asv_type,
            otu_method=otu_method,
            classify_method=classify_method,
            classify_args=classify_args,
            taxa_filter=taxa_filter,
            otu_args=otu_args,
            join_args=join_args
        )

        with open(self.pre_map) as f:
            line = f.readline()
        if categories:
            categories = categories.split(',')
        else:
            categories = re.findall(r'Category\d*', line, flags=re.IGNORECASE)
        if orders:
            orders = orders.split(',')
        else:
            orders = re.findall(r'Order\d*', line, flags=re.IGNORECASE)
        logger.info("The categories are: {}".format(str(categories)))
        logger.info("The Orders are: {}".format(str(orders)))
        self.categories = categories
        self.orders = orders
        self.visualize_kwargs = visualize_kwargs

    @soar_outpath(
        manifest_path="{data_dir}/manifest.txt",
        mapping_file="{data_dir}/sample-metadata.tsv"
    )
    def create_manifest(self):
        self.system("""
{python3_path} {bayegy_home}/write_manifest.py \
  -i {raw_dir} -m {pre_map} -o {data_dir} -f '{fwd_pattern}' \
  -r '{rev_pattern}' -s '{sam_pattern}'""")

    @soar_outpath(demux_qza="{data_dir}/demux.qza")
    def import_raw(self):
        rtype = dict(
            pe="SampleData[PairedEndSequencesWithQuality]",
            se="SampleData[SequencesWithQuality]",
            jpe="SampleData[JoinedSequencesWithQuality]"
        )
        rformat = dict(
            pe="PairedEndFastqManifestPhred33",
            se="SingleEndFastqManifestPhred33",
            jpe="SingleEndFastqManifestPhred33"
        )
        with qiime2():
            self.system("""
qiime tools import --type '{type}' \
  --input-path {manifest_path} --output-path {data_dir}/demux.qza \
  --input-format {format}&& \
qiime demux summarize --i-data {data_dir}/demux.qza --o-visualization {demux_stat_dir}/demux.qzv&&\
qiime tools export --input-path {demux_stat_dir}/demux.qzv --output-path {demux_stat_dir}
""", type=rtype[self.end_type], format=rformat[self.end_type])

    @soar_outpath(
        otu_table="{data_dir}/table-dada2.qza",
        rep_seqs="{data_dir}/rep-seqs-dada2.qza",
        otu_status="{otu_stat_dir}/stats-dada2.qzv"
    )
    def dada2(self):
        data = self.get_demux_data('{demux_qza}')
        dada2_context = dict(
            trunc_f=self.get_trunc_len(data[1]),
            trim_f=len(self.fwd_primer) + 6 if self.fwd_primer else 26
        )
        with qiime2():
            if self.end_type == "pe":
                dada2_context.update(
                    trunc_r=self.get_trunc_len(data[2]),
                    trim_r=len(self.rev_primer) + 6 if self.rev_primer else 26
                )
                otu_args = """--i-demultiplexed-seqs {demux_qza} \
  --p-trunc-len-f {trunc_f} --p-trim-left-f {trim_f} \
  --p-trunc-len-r {trunc_r} --p-trim-left-r {trim_r} \
  --o-representative-sequences {data_dir}/rep-seqs-dada2.qza \
  --o-table {data_dir}/table-dada2.qza  \
  --o-denoising-stats {data_dir}/stats-dada2.qza \
  --p-n-threads {threads} --verbose"""
                if self.otu_args:
                    otu_args = CommandArgs.extend(otu_args, self.otu_args)
                self.system("qiime dada2 denoise-paired " + otu_args, **dada2_context)
            elif self.end_type == "se":
                otu_args = """--i-demultiplexed-seqs {demux_qza} \
  --p-max-ee 50 --p-trunc-len {trunc_f} \
  --p-trunc-q 0 --p-trim-left {trim_f} \
  --o-representative-sequences {data_dir}/rep-seqs-dada2.qza \
  --o-table {data_dir}/table-dada2.qza \
  --o-denoising-stats {data_dir}/stats-dada2.qza \
  --p-n-threads {threads} --verbose"""
                if self.otu_args:
                    otu_args = CommandArgs.extend(otu_args, self.otu_args)
                self.system("qiime dada2 denoise-single " + otu_args, **dada2_context)
            else:
                raise Exception(
                    "Dada2 can not handler this type ({}) of raw data".format(self.end_type)
                )
            self.system("""
qiime metadata tabulate --m-input-file {data_dir}/stats-dada2.qza \
 --o-visualization {otu_stat_dir}/stats-dada2.qzv""")

    @staticmethod
    def get_trunc_len(data, q="25%", t=20):
        for length, obj in data.items():
            if obj[q] < t:
                return int(length) - 1
        return int(length)

    @soar_outpath(
        otu_table="{data_dir}/table-deblur.qza",
        rep_seqs="{data_dir}/rep-seqs-deblur.qza",
        otu_status="{otu_stat_dir}/stats-deblur.qzv"
    )
    def deblur(self):
        with qiime2():
            if self.end_type == "pe":
                if self.fwd_primer and self.rev_primer:
                    self.system("""
qiime cutadapt trim-paired --i-demultiplexed-sequences {demux_qza} \
 --p-front-f {fwd_primer} --p-front-r {rev_primer} \
 --o-trimmed-sequences {data_dir}/trimmed-demux.qza""")
                    self.set_path(force=False, demux_qza="{data_dir}/trimmed-demux.qza".format(**self.context))
                join_args = " --i-demultiplexed-seqs {demux_qza} --verbose \
 --p-allowmergestagger --p-truncqual 5 --o-joined-sequences {data_dir}/demux-joined.qza"
                if self.join_args:
                    join_args = CommandArgs.extend(join_args, self.join_args)
                self.system("qiime vsearch join-pairs " + join_args)
                self.set_path(force=False, demux_qza="{data_dir}/demux-joined.qza".format(**self.context))
            if self.end_type == "pe" or self.end_type == "jpe":
                self.system("""
qiime quality-filter q-score-joined --i-demux {demux_qza} \
 --o-filtered-sequences {data_dir}/demux-joined-filtered.qza \
 --o-filter-stats {data_dir}/demux-joined-filter-stats.qza""")
                self.set_path(force=False, demux_qza="{data_dir}/demux-joined-filtered.qza".format(**self.context))
            elif self.end_type == "se":
                if self.fwd_primer:
                    self.system("""
qiime cutadapt trim-single --i-demultiplexed-sequences {demux_qza} \
 --p-front {fwd_primer} --o-trimmed-sequences {data_dir}/trimmed-demux.qza""")
                    self.set_path(force=False, demux_qza="{data_dir}/trimmed-demux.qza".format(**self.context))
                self.system("""
qiime quality-filter q-score --i-demux {demux_qza} \
 --o-filtered-sequences {data_dir}/demux-filtered.qza \
 --o-filter-stats {data_dir}/demux-filter-stats.qza """)
                self.set_path(force=False, demux_qza="{data_dir}/demux-filtered.qza".format(**self.context))
            else:
                raise Exception("Deblur can not handler this type of raw data")
            trim_length = self.get_deblur_trim_len(self.get_demux_data("{demux_qza}")[1])
            # trim_length = "-1"
            deblur_cmd = "qiime deblur denoise-16S" if self.asv_type == "16s" else "qiime \
 deblur denoise-other --i-reference-seqs {db_rep_set}"
            deblur_args = """--i-demultiplexed-seqs {demux_qza} --p-trim-length {trim_length} \
 --p-sample-stats --o-representative-sequences {data_dir}/rep-seqs-deblur.qza \
 --o-table {data_dir}/table-deblur.qza \
 --o-stats {data_dir}/stats-deblur.qza --p-jobs-to-start {threads} --verbose"""
            if self.otu_args:
                deblur_args = CommandArgs.extend(deblur_args, self.otu_args)
            self.system(" ".join([deblur_cmd, deblur_args]), trim_length=trim_length)
            self.system("""
qiime deblur visualize-stats --i-deblur-stats {data_dir}/stats-deblur.qza \
 --o-visualization {otu_stat_dir}/stats-deblur.qzv""")

    @staticmethod
    def get_deblur_trim_len(data):
        total = data['0']['count']
        for length, obj in data.items():
            if obj['count'] / total < 0.9:
                return int(length) - 1
        return int(length)

    def generate_otu(self):
        if self.otu_method == "dada2":
            self.dada2()
        else:
            self.deblur()

    def get_demux_data(self, demux_qza):
        demux = demux_qza.format(**self.context)
        with qiime2():
            self.system("""
mkdir -p {tmp_dir}/get_demux_data;
qiime demux summarize --i-data {demux} --o-visualization {tmp_dir}/get_demux_data/demux.qzv&&\
qiime tools export --input-path {tmp_dir}/get_demux_data/demux.qzv \
 --output-path {tmp_dir}/get_demux_data/""", demux=demux)
        data = self.parse_jsonp("{tmp_dir}/get_demux_data/data.jsonp".format(**self.context))
        self.system("rm -r {tmp_dir}/get_demux_data")
        return data

    def train_classifier(self):
        if self.classify_method == "sklearn" and self.context.get("db_classifier") and \
                (self.otu_method == "dada2" or self.asv_type == "16s"):
            return
        with qiime2():
            if not self.db_rep_set.endswith("qza") or not self.db_taxonomy.endswith('qza'):
                self.system("""
qiime tools import --type 'FeatureData[Sequence]'   --input-path {db_rep_set} \
 --output-path {classifier_dir}/rep-set.qza&&\
qiime tools import  --type 'FeatureData[Taxonomy]'  --input-format HeaderlessTSVTaxonomyFormat \
 --input-path {db_taxonomy} --output-path {classifier_dir}/ref-taxonomy.qza""")
                self.set_path(
                    force=False,
                    db_rep_set="{classifier_dir}/rep-set.qza".format(**self.context),
                    db_taxonomy="{classifier_dir}/ref-taxonomy.qza".format(**self.context)
                )
            if self.fwd_primer and self.rev_primer:
                self.system("""
qiime feature-classifier extract-reads  --i-sequences {db_rep_set} \
 --p-f-primer {fwd_primer}   --p-r-primer {rev_primer} --o-reads {classifier_dir}/ref-seqs.qza""")
                self.set_path(
                    force=False,
                    db_rep_set="{classifier_dir}/ref-seqs.qza".format(**self.context)
                )
            if self.classify_method == 'sklearn':
                self.system("""
qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads {db_rep_set} \
 --i-reference-taxonomy {db_taxonomy} \
 --o-classifier {classifier_dir}/classifier.qza""")
                self.set_path(
                    force=False,
                    db_classifier="{classifier_dir}/classifier.qza".format(**self.context)
                )

    @soar_outpath(
        taxonomy="{data_dir}/taxonomy.withCandM.qza"
    )
    def classify(self):
        with qiime2():
            if self.classify_method == "sklearn":
                classify_args = "--verbose --p-confidence 0.7 --p-n-jobs 1  \
 --i-classifier {db_classifier}  --i-reads {rep_seqs} \
 --o-classification {data_dir}/taxonomy.withCandM.qza"
                classify_cmd = "qiime feature-classifier classify-sklearn"
            else:
                classify_args = "--i-query {rep_seqs}  \
 --i-reference-reads {db_rep_set}  --i-reference-taxonomy {db_taxonomy} \
 --p-perc-identity 0.7 --o-classification {data_dir}/taxonomy.withCandM.qza"
                classify_cmd = "qiime feature-classifier  classify-consensus-vsearch"
            if self.classify_args:
                classify_args = CommandArgs.extend(classify_args, self.classify_args)
            self.system(" ".join([classify_cmd, classify_args]))

    @soar_outpath(
        otu_table="{data_dir}/table.qza",
        rep_seqs="{data_dir}/rep-seqs.qza",
        taxonomy="{data_dir}/taxonomy.qza"
    )
    def taxonomy_filter_stat(self):
        with qiime2():
            self.system("""
{bayegy_home}/filter_source.py -i {otu_table} -r {rep_seqs} \
 -t {taxonomy} -f '{taxa_filter}' -o {data_dir} --raw {data_dir} \
 --tsv {stat_dir}/feature-table.taxonomy.xls \
 --even {stat_dir}/feature-table.taxonomy.even.xls&&\
biom convert -i {stat_dir}/feature-table.taxonomy.xls -o {stat_dir}/feature-table.taxonomy.biom \
 --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy && \
{perl_path} {bayegy_home}/stat_otu_tab.pl {stat_dir}/feature-table.taxonomy.xls \
 -prefix {tmp_dir}/ -spestat {stat_dir}/classified_stat_relative.xls && \
{perl_path} {bayegy_home}/bar_diagram.pl -table {stat_dir}/classified_stat_relative.xls \
 -style 1 -x_title "Sample Name" -y_title "Sequence Number Percent" \
  -right -textup -rotate='-45' --y_mun 1,7 > {tmp_dir}/Classified_stat_relative.svg && \
convert {tmp_dir}/Classified_stat_relative.svg {stat_dir}/Classified_stat_relative.png && \
rm -r {tmp_dir}/*""")

    def tree(self):
        generate_tree(
            rep_seqs_qza=self.rep_seqs,
            nwk_out=os.path.join(self.rep_seq_dir, '../'),
            qza_out=self.data_dir
        )

    def put_otu_results(self):
        with qiime2():
            self.system("""
qiime tools export --input-path {otu_status} --output-path {otu_stat_dir}&&\
qiime tools export --input-path {rep_seqs} --output-path {rep_seq_dir}/../&&\
qiime feature-table tabulate-seqs --i-data {rep_seqs} \
 --o-visualization {rep_seq_dir}/rep-seqs.qzv &&\
qiime tools export --input-path {rep_seq_dir}/rep-seqs.qzv \
 --output-path {rep_seq_dir}/;
for f in $(find {stat_dir} -type f -name "index.html");
    do echo $f; dir=$(dirname $f);
    mv $f ${{dir}}/Summary_请点此文件查看.html;
done;""")

    @staticmethod
    def parse_jsonp(jsonp):
        with open(jsonp) as f:
            string = f.read()
        string = re.sub(r"^[^\(]*\(", "[", string)
        string = re.sub(r"\)[^\)]*$", "]", string)
        return json.loads(string)

    def set_colors(self, colors=False):
        if colors:
            colors_list = '\n'.join([colors[k] for k in sorted(colors.keys(), key=str.lower)])
            colors_list_file = os.path.join(self.out_dir, 'group_color.list')
            with open(colors_list_file, 'w') as f:
                f.write(colors_list)
        else:
            colors_list_file = "{bayegy_home}/piputils/group_color.list".format(**self.context)
        self.system(
            "{python3_path} {bayegy_home}/piputils/write_colors_plan.py -i {mapping_file} -c {cgs} \
            -p {colors_list_file} -o {out_dir}/colors_plan.json", cgs=','.join(self.categories),
            colors_list_file=colors_list_file)
        os.environ['COLORS_PLAN_PATH'] = os.path.join(self.out_dir, 'colors_plan.json')

    def visualize(self):
        self.set_colors()
        for category in self.categories:
            out_dir = os.path.join(self.result_dir, category)
            v = VisualizeAmplicon(
                otu_table=self.otu_table,
                rep_seqs=self.rep_seqs,
                taxonomy=self.taxonomy,
                mapping_file=self.mapping_file,
                category=category,
                out_dir=out_dir,
                processors=self.threads,
                asv_type=self.asv_type,
                orders=self.orders,
                **self.visualize_kwargs
            )
            v.visualize()

    def run(self):
        self.create_manifest()
        self.import_raw()
        self.train_classifier()
        self.generate_otu()
        self.classify()
        self.taxonomy_filter_stat()
        self.tree()
        self.put_otu_results()
        self.visualize()
        print("AmpliconPipeline is finished.")
        os.system("date")
