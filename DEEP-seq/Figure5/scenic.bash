#运行change.py
python change.py

##运行pySCENIC
# 不同物种的数据库不一样，这里是人类是human 
dir=/data01/xgt/pyscenic/cisTarget_databases/ #改成自己的目录
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
# 一定要保证上面的数据库文件完整无误哦 
input_loom=./sample.loom
ls $tfs  $feather  $tbl  

#2.1 grn
pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$tfs #转录因子文件

#2.2 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 10  \
--mask_dropouts

#2.3 AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 10
