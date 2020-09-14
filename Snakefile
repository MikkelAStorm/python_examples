# modify snp data with plink
# filter the data to remove SNPs that are missing in more than $maf of the subjects and where minor 
# allele frequenc is less than $gen, To keep power of analysis high. These limits might be made tighter 
# if needed for the analysis. 
rule plink:
    input:
        bed=config["gen_data"] + "omniexp24_merge.bed",
        bim=config["gen_data"] + "omniexp24_merge.bim",
        fam=config["gen_data"] + "omniexp24_merge.fam", 
        ids=config["ids"]
    params:
        maf=config["maf"],
        geno=config["geno"]
    output:
        snp=config["data_dir"] + 'snp.bed'
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 3600
    shell:
        'export file={output.snp}; export snp=$(echo "${{file%%.*}}"); plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --maf {params.maf} --geno {params.geno} --make-bed --keep {input.ids} --snps-only --out $snp'

# create .raw and .map file
rule plink_raw:
    input:
        bed=config["data_dir"] + 'snp.bed',
        bim=config["data_dir"] + 'snp.bim',
        fam=config["data_dir"] + 'snp.fam',
    output:
        raw=config["data_dir"] + 'snp.raw',
        map=config["data_dir"] + 'snp.map'
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 3600
    shell:
        'export file={output.raw}; export snp=$(echo "${{file%%.*}}"); plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --recodeAD --out $snp; plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --recode12 --out $snp'

# make vcf file 
rule plink_vcf:
    input:
        bed=config["data_dir"] + 'snp.bed',
        bim=config["data_dir"] + 'snp.bim',
        fam=config["data_dir"] + 'snp.fam'
    output:
        vcf=config["data_dir"] + 'snp.vcf'
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 3600
    shell:
        'export file={output.vcf}; export snp=$(echo "${{file%%.*}}"); plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --recode vcf --out $snp'

# compress snp.vcf
rule bgizp_vcf:
    input:
        vcf=config["data_dir"] + 'snp.vcf'
    output:
        vcf=config["data_dir"] + 'snp.vcf.gz',
        tbi=config["data_dir"] + 'snp.vcf.gz.tbi'
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 3600
    shell:
        'bgzip -c {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf}'


# transpose the raw file
rule t_raw:
    input:
        raw=config["data_dir"] + 'snp.raw',
    output:
        raw=config["data_dir"] + 'snp_transpose.raw'
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 3600
    shell:
        '/home/projects/ip_10000-1/people/mamsto/scripts/transpose.sh {input.raw} > {output.raw}'


# make pheontype files
rule phenotype:
    input:
        ids=config["ids"],
        mrna=config["mrna"]
    output:
        subt=config["data_dir"] + "txi_norm_counts_subt.txt",
        counts=config["data_dir"] + "txi_norm_counts.txt",
        covariates=config["data_dir"] + "covariates_subt_numeric.txt",
        covariates2=config["data_dir"] + "covariates_subt_numeric.cov"
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 7200
    shell:
        'export file={output.subt}; export data_dir=$(dirname $file);/home/projects/ip_10000-1/people/mamsto/scripts/01_tximport_phenotype.R {input.ids} $data_dir"/" {input.mrna}'

# make count bed file
rule counts_bed:
    input:
        count=config["data_dir"] + "txi_norm_counts_subt.txt",
        pos=config["data"] + "gene_pos.txt"
    output:
        bed=config["data_dir"] + "txi_norm_counts_subt.bed.gz",
        index=config["data_dir"] + "txi_norm_counts_subt.bed.gz.tbi"
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 3600
    shell:
        '/home/projects/ip_10000-1/people/mamsto/scripts/01_phenotype_bed.sh {input.count} {input.pos} {output.bed}'


# matrix eqtl
rule matrixEQTL_lin:
    input:
        snp=config["data_dir"] + 'snp_transpose.raw',
        counts=config["data_dir"] + "txi_norm_counts_subt.txt",
        covariates=config["data_dir"] + "covariates_subt_numeric.txt"
    output:
        out=config["result_dir"] + "matrixEQTL_cis" + config["cis_dist"] + "_modelLINEAR.txt"        
    params:
        cis=config["cis_dist"]
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 7200
    shell:
        '/home/projects/ip_10000-1/people/mamsto/scripts/02_matrixEQTL.R {input.snp} {input.counts} {input.covariates} {output.out} {params.cis} modelLINEAR '

# matrix eqtl anova
rule matrixEQTL_anova:
    input:
        snp=config["data_dir"] + 'snp_transpose.raw',
        counts=config["data_dir"] + "txi_norm_counts_subt.txt",
        covariates=config["data_dir"] + "covariates_subt_numeric.txt"
    output:
        out=config["result_dir"] + "matrixEQTL_cis" + config["cis_dist"] + "_modelANOVA.txt"
    params:
        cis=config["cis_dist"]
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 7200
    shell:
        '/home/projects/ip_10000-1/people/mamsto/scripts/02_matrixEQTL.R {input.snp} {input.counts} {input.covariates} {output.out} {params.cis} modelANOVA'



# run fastQTL
rule fastQTL:
    input:
        snp=config["data_dir"] + 'snp.vcf.gz',
        counts=config["data_dir"] + 'txi_norm_counts_subt.bed.gz'
    output:
        out=config["result_dir"] + "fastQTL_cis" + config["cis_dist"] + ".tab" 
    params:
        cis=config["cis_dist"]
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 7200
    shell:
        '/home/projects/ip_10000-1/people/mamsto/scripts/03_fastQLT.py {input.snp} {input.counts} {params.cis} {output.out}'


# add p value to he fastQLT
rule fastQTL_p:
    input:
        config["result_dir"] + "fastQTL_cis" + config["cis_dist"] + ".tab"
    output:
        out=config["result_dir"] + "fastQTL_cis" + config["cis_dist"] + "_pval.tab"
    threads: 1
    resources:
        mem_gb = 50,
        runtime = 7200
    shell:
        '/home/projects/ip_10000-1/people/mamsto/scripts/04_corrected_p.R {input} {output.out}'

