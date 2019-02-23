**Simplified Bioinformatics Module (SBM)**

**使用说明**

**概述**

本模块是一个对生物信息学研究中，经常遇到的Fasta、Sequence、Table等格式的文件进行处理的简单、易用、人性化的工具，旨在使相关的分析更加方便、快捷。

通过创建Fasta、Sequence、Table等对象和调用不同的方法，可以实现如下功能：

- Fasta文件：通过identifier快速得到对应的具体序列，支持不完整的identifier（但必须唯一）；

  例如：

  `fa = Fasta('example/fasta.fas')`

  