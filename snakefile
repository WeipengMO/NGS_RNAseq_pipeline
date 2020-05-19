'''
@Author       : windz
@Date         : 2020-05-18 20:38:48
@LastEditTime : 2020-05-19 15:51:44
@Description  : 

步骤顺序： Header_Hisat2， MarkDump, StatDuplication, , RemoveDump_StringtieForBallgown,  RemoveDump_ExtractRpkm, RemoveDump_MergeGeneRpkm, RemoveDump_Extract_Gene_Count, Cal_IR_By_Splicing_Read_Type_From_Bam
'''


import os
import glob


configfile: 'config.yml'

'''
📌 把数据放在origin_data/目录里面，文件目录如下：
.
├── config.yml
├── origin_data
│   ├── mutant-CB-1-1
│   │   ├── mutant-CB-1-1_1.fastq.gz
│   │   ├── mutant-CB-1-1_2.fastq.gz
│   │   └── SRR8757099
│   ├── mutant-CB-1-2
│   │   ├── mutant-CB-1-2_1.fastq.gz
│   │   ├── mutant-CB-1-2_2.fastq.gz
│   │   └── SRR8757100
│   ├── mutant-CB-2-1
│   │   ├── mutant-CB-2-1_1.fastq.gz
│   │   ├── mutant-CB-2-1_2.fastq.gz
│   │   └── SRR8757101
'''

SAMPLE_NAME = os.listdir('./origin_data')


rule all:
    input:
        expand('stringtie/gene_count_matrix.csv', sample_name=SAMPLE_NAME),
        expand('irratio_all/{sample_name}.irratio.txt', sample_name=SAMPLE_NAME)

# 👉比对
rule run_hisat2:
    input:
        fq1='origin_data/{sample_name}/{sample_name}_1.fastq.gz',
        fq2='origin_data/{sample_name}/{sample_name}_2.fastq.gz'
    output:
        'align_data/{sample_name}/{sample_name}.sorted.bam'
    params:
        genome=config['genome']
    threads: 64
    shell:
        '''
hisat2 -x {params.genome} -p {threads} --min-intronlen 20 --max-intronlen 10000 --dta --time -1 {input.fq1} -2 {input.fq2} | samtools sort -@ {threads} -O bam -o {output} - && samtools index {output}
        '''

# 👉去重
rule MarkDuplicates:
    input:
        'align_data/{sample_name}/{sample_name}.sorted.bam'
    output:
        'markdup_data/{sample_name}/{sample_name}.sorted.rmdup.bam'
    threads: 1
    shell:
        '''
java -jar /public/apps/picard_2.20.2/picard.jar MarkDuplicates REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 I={input} O={output} M={output}.markdump.txt && samtools index {output} 
        '''

# 👉转录本拼接
rule run_stringtie:
    input:
        'markdup_data/{sample_name}/{sample_name}.sorted.rmdup.bam'
    output:
        'stringtie/{sample_name}/{sample_name}.gene.gtf'
    threads: 16
    params:
        gff=config['ann_gff'],
        gene_abund='stringtie/{sample_name}/{sample_name}.gene_abund.tab'
    shell:
        '''
stringtie -A {params.gene_abund} -e --rf -B -p {threads}  -G {params.gff} -o {output} {input}
        '''

# 👉提取表达量
rule extract_rpkm:
    input:
        'stringtie/{sample_name}/{sample_name}.gene.gtf'
    output:
        mrna='stringtie/{sample_name}/{sample_name}.RNA.rpkm.txt',
        gene='stringtie/{sample_name}/{sample_name}.gene.rpkm.txt'
    params:
        label='{sample_name}',
        path='stringtie/{sample_name}/'
    threads: 1
    shell:
        '''
~/miniconda3/envs/r/bin/Rscript script/extract_rpkm_from_ballgown.R {params.label} {params.path} {output.mrna} {output.gene}
        '''


# 👉整合stringtie结果，提取基因/转录本count
# prepDE.py 为stringtie的一个脚本
# /public/home/mowp/softwares/bio/bin/prepDE.py
rule extract_gene_count:
    input:
        [f'stringtie/{sample_name}/{sample_name}.RNA.rpkm.txt' for sample_name in SAMPLE_NAME]
    output:
        gene_count_matrix='stringtie/gene_count_matrix.csv',
        transcript_count_matrix='stringtie/transcript_count_matrix.csv'
    threads: 1
    params:
        in_dir='stringtie'
    shell:
        '''
 prepDE.py -g {output.gene_count_matrix} -t {output.transcript_count_matrix} -i {params.in_dir}
        '''


# 👉计算IR比例
# 🤣太难了，有空慢慢学习
# 把阿布的buseq里面用到的函数整在了脚本之中
rule Cal_IR_By_Splicing_Read_Type_From_Bam:
    input:
        'markdup_data/{sample_name}/{sample_name}.sorted.rmdup.bam'
    output:
        'irratio_all/{sample_name}.irratio.txt'
    threads: 1
    params:
        intron_pos='supplementary_data/intron_pos.repr.txt'
    shell:
        '''
python script/cal_ir_ratio.from_bam.by_splicing_Read_type.version2.py {input} {params.intron_pos} {output} 1 4 1
        '''