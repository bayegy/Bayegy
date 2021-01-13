# install qiime1
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

# still need to install  backports.functools-lru-cache
source activate qiime1
pip install backports.functools-lru-cache

# test qiime1
print_qiime_config.py -t