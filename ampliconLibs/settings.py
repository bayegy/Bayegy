from MetaGenome.pipconfig.settings import path as metagenomePath


ampliconPath = {
    # amplicon special path
    "qiime1_home": "/home/bayegy/pipelines/metagenome/miniconda2/envs/qiime1",
    "picrust2_home": "/home/bayegy/pipelines/metagenome/miniconda2/envs/picrust2",
    "picrust2_database": "/home/bayegy/pipelines/metagenome/miniconda2/envs/picrust2/lib/python3.6/site-packages/PICRUSt2-2.3.0b0-py3.6.egg/picrust2/default_files",
    "metagenome_home": "/home/bayegy/pipelines/metagenome/MetaGenome",
    "mapfiles_home": "/home/bayegy/Databases/mapfiles"
}


path = metagenomePath.copy()
path.update(**ampliconPath)


db_taxonomy = {
    "16s": "/home/bayegy/Databases/amplicon/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt",
    "18s": "/home/bayegy/Databases/amplicon/silva_18S_fungi_fixed/silva_taxonomy_fungi_fixed.txt",
    "its": "/home/bayegy/Databases/amplicon/UNITE_release/sh_qiime_release_s/sh_taxonomy_qiime_ver7_dynamic_s_01.12.2017.txt"
}
# "16s_silva": "/home/bayegy/Databases/amplicon/silva_16S/taxonomy_silvaContent_kpcofgsFormat_only16S_99.txt"

db_rep_set = {
    "16s": "/home/bayegy/Databases/amplicon/gg_13_8_otus/rep_set/99_otus.fasta",
    "18s": "/home/bayegy/Databases/amplicon/silva_18S_fungi_fixed/silva_132_99_18S.fna",
    "its": "/home/bayegy/Databases/amplicon/UNITE_release/sh_qiime_release_s/sh_refs_qiime_ver7_dynamic_s_01.12.2017.fasta"
}
# "16s_silva": "/home/bayegy/Databases/amplicon/silva_16S/silva_132_99_16S.fna"
