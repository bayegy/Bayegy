from .libShell import qiime2, picrust2
import os
from . import settings
import pandas as pd
import numpy as np
from MetaGenome.mapInfo import MapInfo
import logging
logger = logging.getLogger(__name__)


def generate_tree(rep_seqs_qza, nwk_out, qza_out):
    with qiime2():
        if not rep_seqs_qza.endswith('.qza'):
            cmd0 = "qiime tools import --input-path {rep_seqs_qza} \
 --type FeatureData[Sequence] --output-path {qza_out}/rep_seqs.qza".format(**locals())
            print(cmd0)
            os.system(cmd0)
            rep_seqs_qza = "{qza_out}/rep_seqs.qza".format(**locals())
            if not os.path.exists(rep_seqs_qza):
                raise Exception("Unable to import the seqs, please check the input format.")
        cmd = """
qiime alignment mafft --i-sequences {rep_seqs_qza} --o-alignment {qza_out}/aligned-rep-seqs.qza && \
qiime alignment mask --i-alignment {qza_out}/aligned-rep-seqs.qza --o-masked-alignment {qza_out}/masked-aligned-rep-seqs.qza && \
qiime phylogeny fasttree --i-alignment {qza_out}/masked-aligned-rep-seqs.qza   --o-tree {qza_out}/unrooted-tree.qza && \
qiime phylogeny midpoint-root   --i-tree {qza_out}/unrooted-tree.qza   --o-rooted-tree {qza_out}/rooted-tree.qza && \
qiime tools export --input-path {qza_out}/rooted-tree.qza --output-path {nwk_out}/ && \
mv {nwk_out}/tree.nwk {nwk_out}/tree.rooted.nwk && \
qiime tools export --input-path {qza_out}/unrooted-tree.qza --output-path {nwk_out}/ && \
mv {nwk_out}/tree.nwk {nwk_out}/tree.unrooted.nwk
        """.format(**locals())
        print(cmd)
        os.system(cmd)


def picrust2_pipline(table, seqs, out_dir="./", asv_type="16s", processors=1, min_align=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    database = settings.path.get("picrust2_database")
    if not database:
        raise Exception("picrust2_database not configed!")
    refs = {
        "16s": os.path.join(database, "prokaryotic/pro_ref"),
        "18s": os.path.join(database, "fungi/fungi_18S"),
        "its": os.path.join(database, "fungi/fungi_ITS")
    }
    ref = refs[asv_type]
    pathway_maps = {
        "16s": os.path.join(database, "pathway_mapfiles/metacyc_path2rxn_struc_filt_pro.txt"),
        "18s": os.path.join(database, "pathway_mapfiles/metacyc_path2rxn_struc_filt_euk.txt"),
        "its": os.path.join(database, "pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt"),
    }
    pathway_map = pathway_maps[asv_type]
    ec_numbers = {
        "16s": os.path.join(database, "prokaryotic/ec.txt.gz"),
        "18s": os.path.join(database, "fungi/ec_18S_counts.txt.gz"),
        "its": os.path.join(database, "fungi/ec_ITS_counts.txt.gz")
    }
    ec_number = ec_numbers[asv_type]
    copy_numbers = {
        "16s": os.path.join(database, "prokaryotic/16S.txt.gz"),
        "18s": os.path.join(database, "fungi/18S_counts.txt.gz"),
        "its": os.path.join(database, "fungi/ITS_counts.txt.gz")
    }
    copy_number = copy_numbers[asv_type]
    regroup_map = os.path.join(database, "pathway_mapfiles/ec_level4_to_metacyc_rxn.tsv")
    if asv_type == "16s":
        custom_args = "--min_align {}".format(min_align or 0.8)
    else:
        min_align = min_align or 0.6
        custom_args = "--ref_dir {ref} --pathway_map {pathway_map} --regroup_map {regroup_map}  \
 --custom_trait_tables {ec_number} --reaction_func {ec_number} \
 --marker_gene_table {copy_number} --min_align {min_align}".format(**locals())

    df = pd.read_csv(table, sep="\t", index_col=0)
    if np.dtype("O") == df.dtypes[-1]:
        df = df.drop(df.columns[-1], axis=1)
    table = os.path.join(out_dir, "input_table_preprocessed.tsv")
    df.to_csv(table, sep="\t")
    ec_prefix = "EC" if asv_type == "16s" else os.path.basename(ec_number).rstrip('.gz')
    cmd = """
picrust2_pipeline.py -s {seqs} -i {table} -o {out_dir}/initial_res -p {processors} \
{custom_args} --skip_nsti --skip_norm --remove_intermediate --verbose && \
gunzip {out_dir}/initial_res/*gz {out_dir}/initial_res/*/*gz && \
mv {out_dir}/initial_res/{ec_prefix}_predicted.tsv {out_dir}/OTU.EC.abundance.xls;
mv {out_dir}/initial_res/{ec_prefix}_metagenome_out/pred_metagenome_unstrat.tsv {out_dir}/Sample.EC.abundance.xls;
mv {out_dir}/initial_res/pathways_out/path_abun_unstrat.tsv {out_dir}/MetaCyc.pathway.abundance.xls;""".format(**locals())
    print(cmd)
    with picrust2():
        os.system(cmd)

    abc_maps = [
        ("{out_dir}/Sample.EC.abundance.xls", "{picrust2_database}/description_mapfiles/ec_level4_info.tsv.gz"),
        ("{out_dir}/MetaCyc.pathway.abundance.xls", "{picrust2_database}/description_mapfiles/metacyc_pathways_info.txt.gz")
    ]
    if asv_type == "16s":
        cmd1 = """
mv {out_dir}/initial_res/KO_predicted.tsv {out_dir}/OTU.KO.abundance.xls;
mv {out_dir}/initial_res/KO_metagenome_out/pred_metagenome_unstrat.tsv {out_dir}/Sample.KO.abundance.xls;
{perl_path} {metagenome_home}/ConvergeKO2Pathway.pl {out_dir}/Sample.KO.abundance.xls > {out_dir}/KEGG.pathway.abundance.xls;
    """.format(**locals(), **settings.path)
        print(cmd1)
        os.system(cmd1)
        abc_maps.extend([
            ("{out_dir}/KEGG.pathway.abundance.xls", "{mapfiles_home}/kegg/kegg_pathways_merge_levels.tsv"),
            ("{out_dir}/Sample.KO.abundance.xls", "{picrust2_database}/description_mapfiles/ko_info.tsv.gz")
        ])
    mi = MapInfo()
    for abc_table, mapfile in abc_maps:
        abc_table = abc_table.format(**locals())
        if not os.path.exists(abc_table):
            logger.warning("{abc_table} does not exists! function prediction may get an empty result.".format(**locals()))
            break
        mapfile = mapfile.format(**settings.path)
        mi.mapping(abc_table, mapfile)
    cmd2 = "rm -r {table} {out_dir}/initial_res"
    if asv_type == "16s":
        cmd2 += "; {python3_path} {bayegy_home}/summarize_levels.py -i {out_dir}/KEGG.pathway.abundance.xls \
 --output-abs {out_dir}/KEGG.pathway.L{{}}.abundance.xls"
    cmd2 = cmd2.format(**locals(), **settings.path)
    print(cmd2)
    os.system(cmd2)


def all_path_exists(paths):
    for path in paths:
        if not os.path.exists(path):
            return False
    return True


def soar_outpath(**paths):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            paths_dict = {k: p.format(**self.context) for k, p in paths.items()}
            rendered_paths = paths_dict.values()
            paths_exists = [os.path.exists(p) for p in rendered_paths]
            result = None
            if sum(paths_exists) == len(paths_exists):
                print("paths: {} exists, {} was passed!".format(
                    " and ".join(rendered_paths),
                    str(func)
                ))
            else:
                for p, exists in zip(rendered_paths, paths_exists):
                    if exists:
                        os.system("rm -r {}".format(p))
                result = func(self, *args, **kwargs)
            self.set_path(force=False, **paths_dict)
            return result
        return wrapper
    return decorator
