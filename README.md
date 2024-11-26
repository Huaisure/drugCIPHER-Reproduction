# drugCIPHER

本仓库为论文  [《Network-Based Relating Pharmacological and Genomic Spaces for Drug Target Identification》](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011764) 的复现项目，实现了论文中的三种计算药物和靶点相关性的方法，分别为drugCIPHER-TS、drugCIPHER-CS、drugCIPHER-MS。

### 项目介绍

该项目中使用到的数据集为 DrugBank、KEGG Compound Dataset 以及 STRING，项目中涉及到的数据以`.csv` 格式存储，一些提前计算好的相关性系数的矩阵以`.npy`格式存储。

`data_preprocessing` 目录下为对数据集内元素的清洗和映射的预处理操作。

`heatmap` 目录下为根据药物之间的三种分数，TS、CS和GR，绘制热力图的脚本。

`kegg` 目录下为根据kegg API 下载KEGG相关数据集的脚本，以及存储了药物对应的2D化学结构 `.mol` 格式文件。

`PPI` 目录下根据蛋白质关联数据构建ppi网络。

`target_drug_relevance` 内部为计算药物和靶点对应关联度的相关脚本文件。

根目录下的 `drugCIPHER_TS.py`、`drugCIPHER_CS.py` 和`drugCIPHER_MS.py` 为计算数据集中测试集对训练集的三种分数的脚本

论文原件附在根目录下。

### How to Start

```
python drugCIPHER_MS.py
```

会输出测试集中药物与其对应靶点匹配的准确率，分别依据三种评分标准。
