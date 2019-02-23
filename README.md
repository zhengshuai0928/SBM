**Simplified Bioinformatics Module (SBM)**

**使用说明**

**1. 概述**

本模块是一个对生物信息学研究中，经常遇到的Fasta、Sequence、Table等格式的文件进行处理的简单、易用、人性化的工具，旨在使相关的分析更加方便、快捷。

**2. 操作系统与python版本**

Linux; python 3.6

**3. 使用 方式**

```python
from sbm import *
```

通过创建Fasta、Sequence、Table等对象和调用不同的方法，可以实现如下功能：

- Fasta文件：

  - ```python
    >> fa = Fasta('example/fasta.fas')   #创建Fasta对象
    ```

  - 通过identifier快速得到对应的具体序列，支持不完整的identifier（但必须唯一）

  ```python
  >>> fa['>scaffold16_size95689_859_6011']  #获得id为>scaffold16_size95689_859_6011的序列
  'CAGCCTTATTTTAGGAGGCATAACTTTTGAAGCTTCCCCCCACATATTAGGTGGCATTACTTCCAGAGCAATCCCCCCCTTATTAGGAGGCATTGCTTTCGGAGCA'
  >>> fa['scaffold16_size95689']  #不完整但唯一的id
  'CAGCCTTATTTTAGGAGGCATAACTTTTGAAGCTTCCCCCCACATATTAGGTGGCATTACTTCCAGAGCAATCCCCCCCTTATTAGGAGGCATTGCTTTCGGAGCA'
  # 遍历fasta文件中每一个id和对应序列
  for key in fa.keys():
      print(key)
      print(fa[key])
  ```

  

  