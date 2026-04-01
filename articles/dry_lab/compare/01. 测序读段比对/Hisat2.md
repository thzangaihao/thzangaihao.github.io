<!-- vscode-markdown-toc -->
* 1. [HISAT2双端比对流程](#HISAT2)
* 2. [基本命令](#)
* 3. [Main arguments](#Mainarguments)
* 4. [Options](#Options)
	* 4.1. [Input options](#Inputoptions)
	* 4.2. [Alignment options](#Alignmentoptions)
	* 4.3. [Scoring options](#Scoringoptions)
	* 4.4. [Spliced alignment options](#Splicedalignmentoptions)
	* 4.5. [Reporting options](#Reportingoptions)
	* 4.6. [Paired-end options](#Paired-endoptions)
	* 4.7. [Output options](#Outputoptions)
	* 4.8. [SAM options](#SAMoptions)
	* 4.9. [Performance options](#Performanceoptions)
* 5. [Other options](#Otheroptions)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->
# HISAT2
##  1. <a name='HISAT2'></a>HISAT2双端比对流程
双端测序无论是做RNA-seq还是DNA-seq，都有如下情况：

假设极限读长为20bp，不同情况随着片段不同：

1. 完全一致的反向互补平行

        5'-ATAGTCATGCAGTCGTATGC-3'
        20b====================
        3'-TATCAGTACGTCAGCATACG-5'

2. 部分重叠型

        5'-TGCAGTCTAGAGTCGTATTA-3'
        30b==============================
                    TCTCAGCATAATCCGCAACCAT-5'

3. 存在未知序列

        5'-GCAGTCATGCGCGTAGCTAT-3'
        50b==================================================
                                        3'-ACTGCATGCGTACGTGACGT-5'
        此情况存在技术型gap

基本流程：\
&emsp;&emsp;分别定位：先为每个read单独寻找候选比对位置\
&emsp;&emsp;配对验证：检查这些位置是否满足双端约束\
&emsp;&emsp;得分计算：计算整体配对得分

重叠序列的处理：\
当reads有重叠时（情况1，2）：\
&emsp;&emsp;重叠区域的得分处理：\
&emsp;&emsp;1. 一致性检查\
&emsp;&emsp;&emsp;&emsp;重叠区域必须完全匹配\
&emsp;&emsp;&emsp;&emsp;如果重叠部分序列不一致，该配对会被拒绝或罚分\
&emsp;&emsp;2. 得分计算\
&emsp;&emsp;&emsp;&emsp;HISAT2使用全局最优策略：\
&emsp;&emsp;&emsp;&emsp;不重复计算重叠区域的得分\
&emsp;&emsp;&emsp;&emsp;将两个reads视为一个连续片段进行评分\
&emsp;&emsp;&emsp;&emsp;重叠部分只贡献一次得分

&emsp;&emsp;未重叠时的gap处理（情况3）\
&emsp;&emsp;&emsp;&emsp;gap延伸罚分体系：\
&emsp;&emsp;&emsp;&emsp;HISAT2使用仿射gap罚分：\
&emsp;&emsp;&emsp;&emsp;1. gap开放罚分 (--rdg / --rfg)\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;对每个新gap的一次性罚分\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;默认：5分\
&emsp;&emsp;&emsp;&emsp;2. gap延伸罚分 (--rdg / --rfg)\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;对gap每个碱基的罚分\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;默认：3分/碱基\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;罚分计算示例：\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;gap长度 = 10bp\
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;gap罚分 = 开放罚分(5) + 延伸罚分(3 × 10) = 35分

得分组成：\
&emsp;&emsp;$双端总得分 = read1得分 + read2得分 - gap罚分$\
&emsp;&emsp;其中：\
&emsp;&emsp;$read得分 = 匹配得分 - 错配罚分 - 该read内的gap罚分$\
&emsp;&emsp;$gap罚分 = 读段间gap的开放罚分 + 延伸罚分$

至于单端数据，那就很简单了。

##  2. <a name=''></a>基本命令
```bash
hisat2 \
    [options]* \

-x <hisat2-idx> \
    {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} \
    [-S <hit>]
```

##  3. <a name='Mainarguments'></a>Main arguments
`-x <hisat2-idx>`

The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final / etc. looks for the specified index first in the current directory, then in the directory specified in the environment variable..1.ht2hisat2HISAT2_INDEXES

参考基因组索引的基名。basename是任何索引文件的名称，但不包括最终的“/”。HISAT2会首先在当前目录中查找指定的索引，然后在环境变量中指定的目录中查找。

`-1 <m1>`

Comma-separated list of files containing mate 1s (filename usually includes ), e.g. . Sequences specified with this option must correspond file-for-file and read-for-read with those specified in . Reads may be a mix of different lengths. If is specified, will read the mate 1s from the “standard in” or “stdin” filehandle._1-1 flyA_1.fq,flyB_1.fq<m2>-hisat2

双端测序选项。支持多文件映射，以逗号分隔的文件列表，至少要包含一个R1端。使用此选项指定的序列必须与中指定的file-for-file和read-for-read对应。如果指定，将从另一个fastq文件中读取第二个read，即另一端。

`-2 <m2>`

Comma-separated list of files containing mate 2s (filename usually includes ), e.g. . Sequences specified with this option must correspond file-for-file and read-for-read with those specified in . Reads may be a mix of different lengths. If is specified, will read the mate 2s from the “standard in” or “stdin” filehandle._2-2 flyA_2.fq,flyB_2.fq<m1>-hisat2

双端测序选项，即指定另一端的文件，和-1 <m1>一样。

`-U <r>`

Comma-separated list of files containing unpaired reads to be aligned, e.g. . Reads may be a mix of different lengths. If is specified, gets the reads from the “standard in” or “stdin” filehandle.lane1.fq,lane2.fq,lane3.fq,lane4.fq-hisat2

单端测序选项。同样支持多文件映射，以逗号分隔的文件列表，至少要包含一个fastq文件。

`--sra-acc <SRA accession number>`

Comma-separated list of SRA accession numbers, e.g. . Information about read types is available at , where sra-acc is SRA accession number. If users run HISAT2 on a computer cluster, it is recommended to disable SRA-related caching (see the instruction at SRA-MANUAL).--sra-acc SRR353653,SRR353654 \
http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?sp=runinfo&acc=sra-acc&retmode=xml

NCBI子数据库SRA库相关。

`-S <hit>`

File to write SAM alignments to. By default, alignments are written to the “standard out” or “stdout” filehandle (i.e. the console).

用于指定输出的SAM格式比对结果文件的名称和路径。

##  4. <a name='Options'></a>Options
###  4.1. <a name='Inputoptions'></a>Input options
`-q`

Reads (specified with \<m1>, \<m2>, \<s>) are FASTQ files. FASTQ files usually have extension .fq or .fastq. FASTQ is the default format. See also: --solexa-quals and --int-quals.

指定读取的\<m1>, \<m2>, \<s>的文件为FASTQ文件。

FASTQ文件通常有扩展名.Fq或.fastq.FASTQ是默认格式。参见：`--solexa-quals`和`--int-quals`。

`--qseq`

Reads (specified with \<m1>, \<m2>, \<s>) are QSEQ files. QSEQ files usually end in _qseq.txt. See also: --solexa-quals and --int-quals.

读取指定的\<m1>, \<m2>, \<s>的文件为QSEQ文件。

QSEQ文件通常以_qseq.txt结尾。参见：`--solexa-quals`和`--int-quals`。


`-f`

Reads (specified with \<m1>, \<m2>, \<s>) are FASTA files. FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar. FASTA files do not have a way of specifying quality values, so when -f is set, the result is as if --ignore-quals is also set.

读取指定的\<m1>, \<m2>, \<s>文件为FASTA文件（注意这里指的是读段信息reads是以fasta文件记录的）。

FASTA文件通常有扩展名.fa.fasta.mfa.Fna或类似的。FASTA文件没有指定质量值的方法，所以当设置-f时，也要设置参数`--ignore-quals`。

`-r`

Reads (specified with \<m1>, \<m2>, \<s>) are files with one input sequence per line, without any other information (no read names, no qualities). When -r is set, the result is as if --ignore-quals is also set.

读取指定的\<m1>, \<m2>, \<s>文件为FASTQ文件。

注意，如果你的FASTQ中每行只有一个序列时，没有任何其他信息（没有读名称，没有质量），当设置-r时，还要设置参数`--ignore-quals`。

`-c`

The read sequences are given on command line. I.e. \<m1>, \<m2> and <singles> are comma-separated lists of reads rather than lists of read files. There is no way to specify read names or qualities, so -c also implies --ignore-quals.

说明读取序列在命令行中给出。

\<m1>, \<m2>和\<singles>是逗号分隔的读段列表，而不是读文件列表。没有办法指定读取的名称或质量时，设置有`-c`时要设置参数`--ignore-quals`。

`-s/--skip <int>`

Skip (i.e. do not align) the first \<int> reads or pairs in the input.

跳过前\<int>个读段，常与后面的`-u/--upto <int>`组合进行测试、调试或快速检查数据时。

`-u/--upto <int>`

Align the first \<int> reads or read pairs from the input (after the -s/--skip reads or pairs have been skipped), then stop. Default: no limit.

读到第\<int>个读段时停止，用于限制HISAT2处理读段数量的参数，主要用在测试、调试或快速检查数据时。

默认值：无限制。

`-5/--trim5 <int>`

Trim \<int> bases from 5’ (left) end of each read before alignment (default: 0).

从每个读段的5'端（开头）修剪掉指定数量\<int>的碱基。

测序质量通常在读段末端下降，在比对前精确地修剪掉读段两端不需要的部分，从而提高比对的质量和准确性。

`-3/--trim3 <int>`

Trim \<int> bases from 3’ (right) end of each read before alignment (default: 0).

从每个读段的3'端修剪掉指定数量\<int>的碱基。道理同上。

`--phred33`

Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called the “Phred+33” encoding, which is used by the very latest Illumina pipelines.

用于指定测序数据质量值编码格式为phred33格式。2009年以后的Illumina数据常用此格式。默认情况下，HISAT2会尝试自动检测。

`--phred64`

Input qualities are ASCII chars equal to the Phred quality plus 64. This is also called the “Phred+64” encoding.

用于指定测序数据质量值编码格式为phred64格式。2009年以前的Illumina数据常用此格式。

`--solexa-quals`

Convert input qualities from Solexa (which can be negative) to Phred (which can’t). This scheme was used in older Illumina GA Pipeline versions (prior to 1.3). Default: off.

将输入质量从Solexa（可以为负）转换为Phred（不能为负）。

该方案用于较旧的Illumina GA Pipeline版本（1.3之前）。

默认关闭。
    
`--int-quals`

Quality values are represented in the read input file as space-separated ASCII integers, e.g., 40 40 30 40…, rather than ASCII characters, e.g., II?I…. Integers are treated as being on the Phred quality scale unless --solexa-quals is also specified. Default: off.

当fastq文件质量表示为人类可读的特殊样式，像40 40 30 40…这样，用此选修。
    
默认关闭

###  4.2. <a name='Alignmentoptions'></a>Alignment options
`--n-ceil <func>`

Sets a function governing the maximum number of ambiguous characters (usually Ns and/or .s) allowed in a read as a function of read length. For instance, specifying -L,0,0.15 sets the N-ceiling function f to f(x) = 0 + 0.15 * x, where x is the read length. See also: [setting function options]. Reads exceeding this ceiling are filtered out. Default: L,0,0.15.

根据读段长度设置模糊字符数量的过滤阈值，适用于有N出现的读段。

\<func>参数为-L, num, num。其中第一个为截距，第二个为斜率，x为读段长度，最终的f(x)为在读段x中最多允许f(x)个模糊碱基N，超过这个数目的碱基将被舍弃。

默认值为 $f(x) = 0 + 0.15 * x$ ，例如一个长度为100的读段最大允许出现 $f(100) = 0 + 15 = 15$ 个模糊碱基。筛选出去低质量碱基。

`--ignore-quals`

When calculating a mismatch penalty, always consider the quality value at the mismatched position to be the highest possible, regardless of the actual value. I.e. input is treated as though all quality values are high. This is also the default behavior when the input doesn’t specify quality values (e.g. in -f, -r, or -c modes).

使用 `--ignore-quals` 时，所有碱基都被当作高质量处理，无论实际的测序质量如何。因此，在低质量碱基错配时，罚分也会很大。

当使用 `-f`（FASTA格式）、`-r`（原始序列）或 `-c（命令行直接输入序列）模式时，由于没有质量信息，这是默认行为。

`--nofw/--norc`

If --nofw is specified, hisat2 will not attempt to align unpaired reads to the forward (Watson) reference strand. If --norc is specified, hisat2 will not attempt to align unpaired reads against the reverse-complement (Crick) reference strand. In paired-end mode, --nofw and --norc pertain to the fragments; i.e. specifying --nofw causes hisat2 to explore only those paired-end configurations corresponding to fragments from the reverse-complement (Crick) strand. Default: both strands enabled.

用于控制比对时考虑的参考基因组链方向。

`--nofw`为禁止将读段比对到正向链，只考虑反向互补链的比对，为"no forward"的缩写。

`--norc`为禁止将读段比对到反向互补链，只考虑正向链的比对，为"no reverse complement"的缩写。

在双端测序中，这两个参数控制片段来源的链方向。但是使用该参数需要已知库的链的方向性。

该参数默认禁止，Hisat2会搜索两条链。

###  4.3. <a name='Scoringoptions'></a>Scoring options
`--mp MX,MN`

Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. A number less than or equal to MX and greater than or equal to MN is subtracted from the alignment score for each position where a read character aligns to a reference character, the characters do not match, and neither is an N. If --ignore-quals is specified, the number subtracted quals MX. Otherwise, the number subtracted is MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) ) where Q is the Phred quality value. Default: MX = 6, MN = 2.

错配罚分取值范围设置，当发生错配时，Hisat2会根据每一个碱基质量从你设置的罚分范围中计算出一个合适的罚分，其计算公式为：

$$罚分=MN+floor((MX-MN)×(MIN(Q, 40.0)/40.0))$$

当使用参数`--ignore-quals`模式时（即忽略碱基质量），所有罚分都以最大罚分执行。

默认值为MX = 6, MN = 2。

`--sp MX,MN`

Sets the maximum (MX) and minimum (MN) penalties for soft-clipping per base, both integers. A number less than or equal to MX and greater than or equal to MN is subtracted from the alignment score for each position. The number subtracted is MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) ) where Q is the Phred quality value. Default: MX = 2, MN = 1.

软剪辑罚分范围。

软剪辑是指读段的一部分无法比对到参考基因组，但这部分序列仍然保留在SAM/BAM记录中（用S表示）。

$$罚分 = MN + floor( (MX-MN) × (MIN(Q, 40.0)/40.0) )$$

与错配罚分类似，但专门针对软剪辑区域。

与错配罚分类似，但专门针对软剪辑区域。

默认值为MX = 2, MN = 1。

`--no-softclip`

Disallow soft-clipping.

不允许软剪辑存在，当软剪辑存在于该read时，Hisat2会将这个read直接丢弃。

`--np <int>`

Sets penalty for positions where the read, reference, or both, contain an ambiguous character such as N. Default: 1.

设置当读段或参考基因组中出现模糊字符N时的罚分，默认为1。

`--rdg <int1>,<int2>`

Sets the read gap open (\<int1>) and extend (\<int2>) penalties. A read gap of length N gets a penalty of \<int1> + N * \<int2>. Default: 5, 3.

设置读段中间隙的罚分

\<int1>：间隙打开罚分（gap open penalty）\
\<int2>：间隙扩展罚分（gap extend penalty）

$$总罚分 = <int1> + N × <int2>（N为间隙长度）$$

默认值：5,3

值得注意的是，`--rdg <int1>,<int2>`不仅仅用在生物学层面导致的gap上（如reads碱基缺失，fasta插入造成的gap），双端测序时，两个reads之间的技术型gap（看最后面）也是按照这个罚分规则进行的。

`--rfg <int1>,<int2>`

Sets the reference gap open (\<int1>) and extend (\<int2>) penalties. A reference gap of length N gets a penalty of \<int1> + N * \<int2>. Default: 5, 3.

设置参考基因组中间隙的罚分

\<int1>：间隙打开罚分\
\<int2>：间隙扩展罚分

$$总罚分 = <int1> + N × <int2>（N为间隙长度）$$

默认值：5,3

`--score-min L,<base>,<slope>`

Sets a function governing the minimum alignment score needed for an alignment to be considered “valid” (i.e. good enough to report). This is a function of read length. For instance, specifying L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, where x is the read length. See also: [setting function options]. The default is L,0,-0.2.

该参数是HISAT2中一个过滤参数，它用于根据读段长度设置最低比对得分阈值，只有得分高于这个阈值的比对才会被报告。

因为它根据读段长度设置最低比对得分阈值，所以它是一个动态阈值。

$$f(x) = base + slope × x $$

$f(x)$为最终阈值，$base$是基准阈值，$slope × x$是根据读段长度调整的动态阈值。

默认为`L,0,-0.2`。

###  4.4. <a name='Splicedalignmentoptions'></a>Spliced alignment options
`--pen-cansplice <int>`

Sets the penalty for each pair of canonical splice sites (e.g. GT/AG). Default: 0.

典型剪接位点：最常见的剪接信号（GT-AG，GC-AG，AT-AC）

遇到该剪切位点的罚分，默认不罚分。

该参数是罚分修订，更符合“内含子”的跨越。

这些参数只在对RNA-seq数据进行比对时有意义，因为DNA-seq没有剪接。

`--pen-noncansplice <int>`

Sets the penalty for each pair of non-canonical splice sites (e.g. non-GT/AG). Default: 12.

罕见剪切位点罚分，在研究特殊剪切位点时调整。

默认为12分

`--pen-canintronlen <func>`

Sets the penalty for long introns with canonical splice sites so that alignments with shorter introns are preferred to those with longer ones. Default: G,-8,1

设置典型剪接位点的内含子长度罚分。适用于对内含子长度偏好性的罚分。

例如默认Default `G,-8,1`：

$$罚分 = -8 + 1 × log(内含子长度)$$

| 内含子长度 | 计算过程                           | 罚分 |
| ---------- | ---------------------------------- | ---- |
| 100 bp     | -8 + log(100) = -8 + 4.6 = -3.4    | -3.4 |
| 1,000 bp   | -8 + log(1000) = -8 + 6.9 = -1.1   | -1.1 |
| 10,000 bp  | -8 + log(10000) = -8 + 9.2 = 1.2   | 1.2  |
| 100,000 bp | -8 + log(100000) = -8 + 11.5 = 3.5 | 3.5  |

`--pen-noncanintronlen <func>`

Sets the penalty for long introns with noncanonical splice sites so that alignments with shorter introns are preferred to those with longer ones. Default: G,-8,1

设置非典型剪接位点的内含子长度罚分，道理同上。

`--min-intronlen <int>`\
`--max-intronlen <int>`\
Sets minimum intron length. Default: 20\
Sets maximum intron length. Default: 500000

这些参数为 HISAT2 在寻找剪接比对时设置内含子长度的上下边界：\
最小内含子长度：比这更短的不被认为是内含子，默认为20\
最大内含子长度：比这更长的不被认为是单个内含子，默认为50

当这个内含子超过此阈值范围时，Hisat2会将其用gap罚分表，而不是内含子长度偏好罚分表。

不同物种的推荐设置：
| 物种类型 | 最小内含子长度 | 最大内含子长度    | 理由             |
| -------- | -------------- | ----------------- | ---------------- |
| 酵母     | 20-50          | 1,000-2,000       | 内含子通常很短   |
| 昆虫     | 20-50          | 50,000-100,000    | 中等长度内含子   |
| 植物     | 20-50          | 10,000-50,000     | 可变性较大       |
| 哺乳动物 | 20-50          | 500,000-1,000,000 | 含有超长内含子   |
| 人类     | 20-50          | 1,000,000+        | 已知有极长内含子 |

`--known-splicesite-infile <path>`

With this mode, you can provide a list of known splice sites, which HISAT2 makes use of to align reads with small anchors.
You can create such a list using python hisat2_extract_splice_sites.py genes.gtf > splicesites.txt, where hisat2_extract_splice_sites.py is included in the HISAT2 package, genes.gtf is a gene annotation file, and splicesites.txt is a list of splice sites with which you provide HISAT2 in this mode. Note that it is better to use indexes built using annotated transcripts (such as genome_tran or genome_snp_tran), which works better than using this option. It has no effect to provide splice sites that are already included in the indexes.

使用已知的剪接位点信息（.gtf文件）来指导RNA-seq数据的比对。

```bash
# 使用 HISAT2 自带的 Python 脚本提取剪接位点
python hisat2_extract_splice_sites.py genes.gtf > splicesites.txt
```
```bash
# 使用剪接位点文件进行比对
hisat2 -x genome_index -U rna_seq.fastq -S output.sam --known-splicesite-infile splicesites.txt
```

`--novel-splicesite-outfile <path>`

In this mode, HISAT2 reports a list of splice sites in the file :
    chromosome name \<tab> genomic position of the flanking base on the left side of an intron \<tab> genomic position of the flanking base on the right \<tab> strand (+, -, and .)

’.’ indicates an unknown strand for non-canonical splice sites.

输出在比对过程中新发现的剪接位点。

在 RNA-seq 比对过程中，HISAT2 可能会发现一些未被已知注释包含的剪接位点，这些就是"新"剪接位点。

指定一个输出文件路径\<path>，HISAT2 会将新发现的剪接位点写入该文件。

默认使用该选修，只不过在位表明\<path>时并不记录新位点。

`--novel-splicesite-infile <path>`

With this mode, you can provide a list of novel splice sites that were generated from the above option “–novel-splicesite-outfile”.

允许使用之前发现的新剪接位点来指导新的比对分析。

此路径应当和上述\<path>一致。

这是Hisat2动态学习的一部分。

`--no-temp-splicesite`

HISAT2, by default, makes use of splice sites found by earlier reads to align later reads in the same run, in particular, reads with small anchors (<= 15 bp).

The option disables this default alignment strategy.

禁用HISAT2的实时剪接位点发现功能。

默认不使用该参数。

使用 --no-temp-splicesite 后：\
&emsp;&emsp;静态比对：仅使用预先提供的剪接位点信息（如通过 --known-splicesite-infile 指定的）\
&emsp;&emsp;无动态学习：不会利用当前运行中已比对reads发现的剪接位点\
&emsp;&emsp;一致性降低：可能会降低某些复杂剪接事件的检测灵敏度

`--no-spliced-alignment`

Disable spliced alignment.

完全禁用剪接比对。

这个参数明显不用于RNA-seq分析。

当使用此参数时，HISAT2将：\
    不进行任何跨内含子的比对\
    将每个read作为一个连续的片段进行比对\
    忽略所有剪接位点（无论是已知的还是新发现的）\
    如果read跨越了内含子区域，它将无法被比对，或者只能比对上其中一侧的外显子，因为超过阈值被丢弃了

`--rna-strandness <string>`

Specify strand-specific information: the default is unstranded.\
For single-end reads, use F or R.\
&emsp;&emsp;‘F’ means a read corresponds to a transcript.\
&emsp;&emsp;‘R’ means a read corresponds to the reverse complemented counterpart of a transcript.

For paired-end reads, use either FR or RF.\
With this option being used, every read alignment will have an XS attribute tag:\
&emsp;&emsp;’+’ means a read belongs to a transcript on ‘+’ strand of genome.\
&emsp;&emsp;‘-‘ means a read belongs to a transcript on ‘-‘ strand of genome.

(TopHat has a similar option, –library-type option, where fr-firststrand corresponds to R and RF; fr-secondstrand corresponds to F and FR.)

告诉HISAT2你的RNA-seq实验是否是链特异性建库，以及具体的链方向。

参数取值\
&emsp;&emsp;单端测序：\
&emsp;&emsp;&emsp;&emsp;F：read对应于转录本的正向序列\
&emsp;&emsp;&emsp;&emsp;R：read对应于转录本的反向互补序列\
&emsp;&emsp;双端测序：\
&emsp;&emsp;&emsp;&emsp;FR：第一个read是转录本方向，第二个read是反向互补\
&emsp;&emsp;&emsp;&emsp;RF：第一个read是反向互补，第二个read是转录本方向\

启用此参数后：\
&emsp;&emsp;每个比对的read都会添加 XS标签：\
&emsp;&emsp;&emsp;&emsp;+：read属于基因组正链上的转录本\
&emsp;&emsp;&emsp;&emsp;-：read属于基因组负链上的转录本\
&emsp;&emsp;比对时会优先考虑正确的链方向

如果数据是非链特异性：不要使用此参数或使用默认值。该参数也是默认不使用。

`--tmo/--transcriptome-mapping-only`

Report only those alignments within known transcripts.

用于限制只输出已知转录本内的比对结果。

使用此参数必须提供转录本注释信息，在构建索引时使用`--known-splicesite-infile` 预加载转录本信息。

优点：\
&emsp;&emsp;提高转录本定量的准确性\
&emsp;&emsp;减少数据文件大小\
&emsp;&emsp;简化下游分析\

缺点：\
&emsp;&emsp;会丢失新转录本的信息\
&emsp;&emsp;无法发现新的基因或剪接变体\
&emsp;&emsp;依赖注释文件的质量和完整性

`--dta/--downstream-transcriptome-assembly`

Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.

让HISAT2生成适合下游转录本组装工具的比对结果。当计划使用StringTie、Cufflinks等工具进行转录本组装时，启用此选项。

`--dta-cufflinks`

Report alignments tailored specifically for Cufflinks. In addition to what HISAT2 does with the above option (–dta), With this option, HISAT2 looks for novel splice sites with three signals (GT/AG, GC/AG, AT/AC), but all user-provided splice sites are used irrespective of their signals. HISAT2 produces an optional field, XS:A:[+-], for every spliced alignment.

专门为Cufflinks工具优化的比对报告模式。

`--dta-cufflinks` 是 `--dta` 的增强版本，专门针对Cufflinks转录本组装工具进行额外优化：
1. 包含 `--dta 的所有功能：\
&emsp;更长的锚定序列要求\
&emsp;减少短锚定比对\
&emsp;优化转录本组装

2. 额外的剪接位点检测\
&emsp;主动寻找三种经典剪接信号：\
&emsp;GT/AG（最常见的）\
&emsp;GC/AG（次常见的）\
&emsp;AT/AC（U12型内含子）\
&emsp;但用户提供的剪接位点（如通过GTF文件）会无条件使用，不受信号类型限制

3. 强制添加XS标签\
&emsp;为每个剪接比对添加 XS:A:[+-] 字段\
&emsp;明确指示该比对属于哪条链

`--avoid-pseudogene`

Try to avoid aligning reads to pseudogenes. Note this option is experimental and needs further investigation.

这是一个实验性参数，用于帮助减少reads在假基因上的错误比对。

假基因和真基因具有着相似的结构，但是它的核心启动元件等遭到了破坏，无法行使转录翻译等功能。

`--no-templatelen-adjustment`

Disables template length adjustment for RNA-seq reads.

双端测序相关，禁用RNA-seq reads的模板长度调整。

    基因组位置：       exon1=======[内含子]=======exon2
                        ^ Read1位置                ^ Read2位置
    基因组距离：        1000 bp
    实际CDS长度：       100 bp（因为中间有900bp内含子被剪接）

在SAM/BAM格式中对应 TLEN 字段\
默认情况下，HISAT2会报告TLEN ≈ 100 bp，即不使用此参数。

使用 `--no-templatelen-adjustment` 后：\
&emsp;&emsp;计算模板长度时不考虑剪接\
&emsp;&emsp;使用基因组上的原始物理距离\
&emsp;&emsp;报告未调整的模板长度\
&emsp;&emsp;同样的例子会报告TLEN ≈ 1000 bp

###  4.5. <a name='Reportingoptions'></a>Reporting options
`-k <int>`

It searches for at most distinct, primary alignments for each read. Primary alignments mean alignments whose alignment score is equal to or higher than any other alignments. The search terminates when it cannot find more distinct valid alignments, or when it finds , whichever happens first. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM ‘secondary’ bit (which equals 256) set in its FLAGS field. For reads that have more than distinct, valid alignments, hisat2 does not guarantee that the alignments reported are the best possible in terms of alignment score. Default: 5 (linear index) or 10 (graph index).

Note: HISAT2 is not designed with large values for in mind, and when aligning reads to long, repetitive genomes, large could make alignment much slower.-k-k

`-k` 参数控制HISAT2为每条read搜索的最大不同比对位置数量。\
`-k <int>` 指定每条read最多报告的不同的和主要比对的数量。

假设read有n个高于--score-min的比对位点应当被记录，但是-k参数表明会只记录前k个比对位点。实际上，当一个read被记录超过k个比对位点时，将不再执行--score-min比对得分计算以节省算力。注意！Hisat2也不能保证在k个记录后面还有更大的得分位点作为主要比对。

1. 主要比对\
    指比对得分等于或高于任何其他比对的比对结果\
    对于双端reads，得分是两个mate比对得分的总和

2. 不同的比对（次要比对）\
    指比对到基因组上不同位置的比对\
    即使得分相同，但位置不同就算不同的比对

性能考虑：\
&emsp;&emsp;HISAT2不是为大的-k值设计的\
&emsp;&emsp;比对到长而重复的基因组时，大的-k值会显著降低比对速度

质量保证：\
&emsp;&emsp;当reads有超过\<int>个有效比对时\
&emsp;&emsp;HISAT2不保证报告的比对是绝对最优的（按比对得分）\
&emsp;&emsp;大多数表达量定量工具（像RSEM）会默认忽略次要比对

默认值:\
&emsp;&emsp;线性索引：-k 5\
&emsp;&emsp;图索引：-k 10

| 特性       | 线性索引       | 图索引               |
| ---------- | -------------- | -------------------- |
| 参考基础   | 单一参考序列   | 包含变异的图结构     |
| 变异处理   | 忽略或作为错误 | 整合为合法路径       |
| 比对准确性 | 在变异区域较低 | 在变异区域更高       |
| 索引大小   | 较小           | 较大（包含更多信息） |
| 构建复杂度 | 简单           | 复杂                 |
| 默认-k值   | 5              | 10                   |

    线性索引是这样记录参考基因组的：
        5'-ACCTAACTTTCAACT-3'
    图索引是这样的：
        5'-ACCTAACTTTCAACT-3'
                TG   TCC
                C

图索引包含变异信息：将已知的SNP、indel等变异整合到索引中\
图结构：使用图基因组表示，包含多个可能的路径\
图索引更生物学真实：反映群体遗传多样性

线性索引和图索引都是以.fasta记录基因组的，但是可能会有一个.vcf文件记录变异信息。\
只依赖.fasta文件构建的索引为线性索引，结合.fasta和.vcf两个文件的索引是图索引。

`--max-seeds <int>`

HISAT2, like other aligners, uses seed-and-extend approaches. HISAT2 tries to extend seeds to full-length alignments. In HISAT2, --max-seeds is used to control the maximum number of seeds that will be extended. For DNA-read alignment (--no-spliced-alignment), HISAT2 extends up to these many seeds and skips the rest of the seeds. For RNA-read alignment, HISAT2 skips extending seeds and reports no alignments if the number of seeds is larger than the number specified with the option, to be compatible with previous versions of HISAT2. Large values for --max-seeds may improve alignment sensitivity, but HISAT2 is not designed with large values for --max-seeds in mind, and when aligning reads to long, repetitive genomes, large --max-seeds could make alignment much slower. The default value is the maximum of 5 and the value that comes with -k times 2.

该参数控制HISAT2在种子扩展阶段处理的最大种子数量，是影响比对速度和敏感性平衡的重要参数。

种子是进行比对的一个算法概念，不同于NW算法，Hisat2使用BWT转化和类似于BLAST算法的种子-启发式算法，在损失一点精度的同时大幅度提高运算速度。

$$默认值 = max(5, 2 × -k 值),<int>=5, 2$$

&emsp;&emsp;线性索引：$-k 5 → --max-seeds max(5, 10) = 10$\
&emsp;&emsp;图索引：$-k 10 → --max-seeds max(5, 20) = 20$

这个参数在速度和敏感性之间提供重要权衡：\
&emsp;&emsp;小值：快速但可能漏掉真实比对\
&emsp;&emsp;大值：彻底但速度慢，可能引入噪音\
&emsp;&emsp;默认值：在大多数情况下提供良好平衡

`-a/--all`

HISAT2 reports all alignments it can find. Using the option is equivalent to using both --max-seeds and -k with the maximum value that a 64-bit signed integer can represent (9,223,372,036,854,775,807).

极端参数，让HISAT2报告所有能找到的比对，不设任何限制。

使用这个参数相当于同时设置：

    --max-seeds 9223372036854775807   # 64位有符号整数最大值
    -k 9223372036854775807            # 同样的最大值

该参数是一个学术研究工具而非生产分析工具：\
&emsp;&emsp;理论上：确保不遗漏任何可能的比对\
&emsp;&emsp;实践中：计算上不可行，输出难以处理\
&emsp;&emsp;建议：除非有极其特殊的理由，否则不要使用

--secondary

Report secondary alignments.

明确启用次要比对的报告，并报道，优先级要高于k参数。

与 `-k` 参数配合使用时有重叠效果，不使用此参数时的默认行为取决于其他参数设置。

###  4.6. <a name='Paired-endoptions'></a>Paired-end options
`-I/--minins <int>`

The minimum fragment length for valid paired-end alignments.This option is valid only with --no-spliced-alignment. E.g. if is specified and a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied). A 19-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the -I constraint is applied with respect to the untrimmed mates. -I 60

The larger the difference between -I and -X, the slower HISAT2 will run. This is because larger differences between -I and -X require that HISAT2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), HISAT2 is very efficient.

Default: 0 (essentially imposing no minimum)

设置有效的双端比对的最小插入片段长度为\<int>。

仅在使用 `--no-spliced-alignment` 时有效

这个参数主要针对DNA-seq数据或非剪接比对

这个DNA-seq就是将DNA片段黏贴到基因组DNA上，例如Chip-seq染色质免疫共沉淀的DNA映射。

对于RNA-seq数据（默认允许剪接），插入片段长度计算方式不同，概念不同，在这里需要注意理解。

注意，二代测序时，对于一个长一点的测序片段，我们只能测得它两端的一部分，而中间未知的序列就是DNA-seq中所谓的“插入序列”

    read1      -----
    基因组      =====|==============|=====
    read2                           -----
    插入序列        I|==============|X

你看，双端测序的序列并不一定是反向互补平行的。

同样，在测序片段不是很长的情况下，双端reads可能会重叠或反向互补平行。

通常与 `-X` 一起使用。

默认设置为0。

`-X/--maxins <int>`

The maximum fragment length for valid paired-end alignments. This option is valid only with --no-spliced-alignment. E.g. if is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). A 61-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the constraint is applied with respect to the untrimmed mates, not the trimmed mates. -X 100-X

The larger the difference between -I and -X, the slower HISAT2 will run. This is because larger differences between -I and -X require that HISAT2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), HISAT2 is very efficient.

Default: 500.

设置有效的双端比对的最大插入片段长度为\<int>。

和上述参数一起构成区间范围。

    I|==============|X

默认设置为500。

`--fr/--rf/--ff`

The upstream/downstream mate orientations for a valid paired-end alignment against the forward reference strand. E.g., if --fr is specified and there is a candidate paired-end alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears upstream of the reverse complement of mate 1 and all other constraints are met, that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a downstream mate 2 to be forward-oriented. Default: --fr (appropriate for Illumina’s Paired-end Sequencing Assay).

这个参数控制双端测序reads的预期相对方向。

&emsp;&emsp;`--fr` (Forward-Reverse) 默认值

        参考基因组正向链：  5' -------------> 3'
        Read1 (mate1):    --------->      (正向)
        Read2 (mate2):            <---------  (反向互补)
        特点：
            Read1比对到正向链
            Read2比对到反向互补链
            Read1在Read2的上游
            Illumina标准建库方式

&emsp;&emsp;`--rf` (Reverse-Forward) 反向方向：

        参考基因组正向链：  5' -------------> 3'
        Read1 (mate1):            <--------- (反向互补)
        Read2 (mate2):    --------->      (正向)
        特点：
            Read1比对到反向互补链
            Read2比对到正向链
            Read1在Read2的上游

&emsp;&emsp;`--ff` (Forward-Forward) 两个reads都在正向链：

        参考基因组正向链：  5' -------------> 3'
        Read1 (mate1):    -------------> (正向)
        Read2 (mate2):         -------------> (正向)
        特点：
            两个reads都比对到正向链
            Read1在Read2的上游

`--no-mixed`

By default, when hisat2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior.

默认行为（混合模式）

当HISAT2找不到一致的双端比对时：\
&emsp;&emsp;尝试找到两个reads各自的最佳单端比对\
&emsp;&emsp;分别报告这两个单端比对

`--no-discordant`

By default, hisat2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior.

禁用混合报告

如果找不到配对成功的双端比对，就完全不报告该read pair

只输出"完美配对"的结果

###  4.7. <a name='Outputoptions'></a>Output options
`-t/--time`

Print the wall-clock time required to load the index files and align the reads. This is printed to the “standard error” (“stderr”) filehandle. Default: off.

打印加载索引文件和比对读取所需的实际时间。此信息会打印到“标准错误”（“stderr”）文件。

默认：关闭。
    
`--un <path>`, `--un-gz <path>`, `--un-bz2 <path>`

Write unpaired reads that fail to align to file at \<path>. These reads correspond to the SAM records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set. If --un-gz is specified, output will be gzip compressed. If --un-bz2 is specified, output will be bzip2 compressed. Reads written in this way will appear exactly as they did in the input file, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the input.

单端测序数据：保存完全无法比对的reads，-gz，-bz2都是压缩为该文件的压缩包。

`--al <path>`, `--al-gz <path>`, `--al-bz2 <path>`

Write unpaired reads that align at least once to file at \<path>. These reads correspond to the SAM records with the FLAGS 0x4, 0x40, and 0x80 bits unset. If --al-gz is specified, output will be gzip compressed. If --al-bz2 is specified, output will be bzip2 compressed. Reads written in this way will appear exactly as they did in the input file, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the input.

单端测序数据：保存成功比对至少一次的reads，压缩同理。

`--un-conc <path>`, `--un-conc-gz <path>`, `--un-conc-bz2 <path>`

Write paired-end reads that fail to align concordantly to file(s) at \<path>. These reads correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set (depending on whether it’s mate #1 or #2). .1 and .2 strings are added to the filename to distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in \<path>, the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, .1 or .2 are added before the final dot in \<path> to make the per-mate filenames. Reads written in this way will appear exactly as they did in the input files, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the inputs.

双端测序数据：未能一致比对的双端reads
&emsp;&emsp;保存未能找到一致比对的paired-end reads
&emsp;&emsp;会生成两个文件（mate1和mate2）
&emsp;&emsp;压缩同理
&emsp;&emsp;SAM标志：0x4位设置，且0x40或0x80位设置

`--al-conc <path>`, `--al-conc-gz <path>`, `--al-conc-bz2 <path>`

Write paired-end reads that align concordantly at least once to file(s) at \<path>. These reads correspond to the SAM records with the FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set (depending on whether it’s mate #1 or #2). .1 and .2 strings are added to the filename to distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in \<path>, the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, .1 or .2 are added before the final dot in \<path> to make the per-mate filenames. Reads written in this way will appear exactly as they did in the input files, without any modification (same sequence, same name, same quality string, same quality encoding). Reads will not necessarily appear in the same order as they did in the inputs.

双端测序数据：成功一致比对的双端reads
&emsp;&emsp;保存成功找到一致比对的paired-end reads
&emsp;&emsp;会生成两个文件（mate1和mate2）
&emsp;&emsp;压缩同理
&emsp;&emsp;SAM标志：0x4位未设置，且0x40或0x80位设置

`--quiet`

Print nothing besides alignments and serious errors.

输出静默控制，完全静默模式
&emsp;&emsp;只输出：比对结果 + 严重错误
&emsp;&emsp;不输出：进度信息、警告、统计信息等

`--summary-file <path>`

Print alignment summary to this file.

摘要报告控制
&emsp;&emsp;将比对摘要统计写入指定文件路径中
&emsp;&emsp;包含：总reads数、比对率、唯一比对数等

`--new-summary`

Print alignment summary in a new style, which is more machine-friendly.

使用新的摘要格式（更机器友好）

便于后续程序解析处理

`--met-file <path>`

Write hisat2 metrics to file <path>. Having alignment metric can be useful for debugging certain problems, especially performance issues. See also: --met. Default: metrics disabled.

将详细的性能指标写入指定文件

用于调试性能问题

默认关闭
    
`--met-stderr`

Write hisat2 metrics to the “standard error” (“stderr”) filehandle. This is not mutually exclusive with --met-file. Having alignment metric can be useful for debugging certain problems, especially performance issues. See also: --met. Default: metrics disabled.

将性能指标输出到标准错误流

可与 `--met-file` 同时使用

默认关闭

`--met <int>`

Write a new hisat2 metrics record every \<int> seconds. Only matters if either --met-stderr or --met-file are specified. Default: 1.

指定性能指标记录的时间间隔（秒）

默认：1秒记录一次

###  4.8. <a name='SAMoptions'></a>SAM options
`--no-unal`

Suppress SAM records for reads that failed to align.

抑制未比对的SAM记录

不输出任何比对失败的reads的SAM记录

显著减少输出文件大小

`--no-hd`

Suppress SAM header lines (starting with @).

完全抑制SAM头信息

不输出任何以@开头的头行

`--no-sq`

Suppress @SQ SAM header lines.

抑制参考序列信息

不输出@SQ头行（包含参考序列名称和长度）

`--rg-id <text>`

Set the read group ID to <text>. This causes the SAM @RG header line to be printed, with <text> as the value associated with the ID: tag. It also causes the RG:Z: extra field to be attached to each SAM output record, with value set to <text>.

设置读组ID

自动添加@RG头行和RG:Z:标签到每个记录

`--rg <text>`

Add \<text> (usually of the form TAG:VAL, e.g. SM:Pool1) as a field on the @RG header line. Note: in order for the @RG line to appear, --rg-id must also be specified. This is because the ID tag is required by the SAM Spec. Specify --rg multiple times to set multiple fields. See the SAM Spec for details about what fields are legal.

添加读组字段（格式：TAG:VAL）

必须与 `--rg-id` 一起使用

可多次指定添加多个字段

`--remove-chrname`

Remove ‘chr’ from reference names in alignment (e.g., chr18 to 18)

从参考序列名称中移除'chr'前缀

例如：chr1 → 1

`--add-chrname`

Add ‘chr’ to reference names in alignment (e.g., 18 to chr18)

向参考序列名称中添加'chr'前缀

例如：1 → chr1

`--omit-sec-seq`

When printing secondary alignments, HISAT2 by default will write out the SEQ and QUAL strings. Specifying this option causes HISAT2 to print an asterisk in those fields instead.

在次要比对中省略序列和质量信息

用*代替实际的SEQ和QUAL字段

###  4.9. <a name='Performanceoptions'></a>Performance options
`-o/--offrate <int>`

Override the offrate of the index with <int>. If <int> is greater than the offrate used to build the index, then some row markings are discarded when the index is read into memory. This reduces the memory footprint of the aligner but requires more time to calculate text offsets. <int> must be greater than the value used to build the index.

覆盖索引的offrate值

offrate控制索引中标记点的密度\
较高的offrate = 较少标记点 = 较小内存占用，但较慢的定位速度\
较低的offrate = 较多标记点 = 较大内存占用，但较快的定位速度

| offrate | 内存占用 | 搜索速度 |
| ------- | -------- | -------- |
| 1（低） | 高       | 非常快   |
| 2       | 中       | 快       |
| 4       | 低       | 中等     |
| 8（高） | 非常低   | 慢       |

`-p/--threads NTHREADS`

Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases HISAT2’s memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if HISAT2 is linked with the pthreads library.

设置并行搜索线程数

默认：1个线程

比对搜索高度并行，加速接近线性

会增加内存占用

| 线程数 | 相对速度 | 内存增长 |
| ------ | -------- | -------- |
| 1      | 1.0x     | 基准     |
| 4      | ~3.8x    | +200MB   |
| 8      | ~7.5x    | +400MB   |
| 16     | ~14x     | +800MB   |

`--reorder`

Guarantees that output SAM records are printed in an order corresponding to the order of the reads in the original input file, even when -p is set greater than 1. Specifying --reorder and setting -p greater than 1 causes HISAT2 to run somewhat slower and use somewhat more memory then if --reorder were not specified. Has no effect if -p is set to 1, since output order will naturally correspond to input order in that case.

保证输出顺序与输入文件顺序一致

在多线程模式下（-p > 1）特别有用

会稍微降低速度并增加内存使用

`--mm`

Use memory-mapped I/O to load the index, rather than typical file I/O. Memory-mapping allows many concurrent bowtie processes on the same computer to share the same memory image of the index (i.e. you pay the memory overhead just once). This facilitates memory-efficient parallelization of bowtie in situations where using -p is not possible or not preferable.

使用内存映射I/O加载索引

允许多个HISAT2进程共享同一索引的内存映像

节省内存，特别适合并行运行多个实例

##  5. <a name='Otheroptions'></a>Other options
`--qc-filter`

Filter out reads for which the QSEQ filter field is non-zero. Only has an effect when read format is --qseq. Default: off.

过滤QSEQ格式中质量失败的reads

仅当使用--qseq输入格式时有效

过滤掉QSEQ格式中filter字段非零的reads

默认关闭

`--seed <int>`

Use \<int> as the seed for pseudo-random number generator. Default: 0.

设置伪随机数生成器的种子值

默认种子：0

同样的种子确保结果可重现

`--non-deterministic`

Normally, HISAT2 re-initializes its pseudo-random generator for each read. It seeds the generator with a number derived from (a) the read name, (b) the nucleotide sequence, (c) the quality sequence, (d) the value of the --seed option. This means that if two reads are identical (same name, same nucleotides, same qualities) HISAT2 will find and report the same alignment(s) for both, even if there was ambiguity. When --non-deterministic is specified, HISAT2 re-initializes its pseudo-random generator for each read using the current time. This means that HISAT2 will not necessarily report the same alignment for two identical reads. This is counter-intuitive for some users, but might be more appropriate in situations where the input consists of many identical reads.

启用非确定性模式

使用当前时间作为随机种子

`--version`

Print version information and quit.

版本信息

`-h` /` --help`

Print usage information and quit.

帮助信息