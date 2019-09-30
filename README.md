### 覆盖度评估 ###
覆盖度评估，核心代码：
```
# 基于不同的reads，选择各种软件生成bam
samtools depth $BAM -a -a | gzip -c > mapping.depth.gz

# 全基因组coverage评估
python ~/src/depth_plot.py mapping.depth.gz mapping.depth &

# BUSCO基因区coverage评估, run_ref.fa/full*tsv是BUSCO输出的full*tsv文件
python ~/src/depth_plot_busco.py mapping.depth.gz mapping.depth.busco run_ref.fa/full*tsv & 
wait
```
