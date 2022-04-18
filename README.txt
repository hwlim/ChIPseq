ChIP-seeq data analysis pipelines

Snakemake.ChIP_SE:
- Default pipeline that can process both single-end and paired-end data
- Eventually in peak calling and bigwig file generation, paired-end data is processed like single-end by treating each read as a separate read
Snakemake.PE
- Developed for fragment-level analysis of paired-end sequencing data
- Thus, it include fragment level bigwig file generation and domain calling
- Yet, it can do peak calling in single-end mode like Snakemake.ChIP_SE
