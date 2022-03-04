#!/bin/env nextflow

/*
 * Copyright (c) 2021, Diego Garrido-Martín
 *
 * Simulation setting to benchmark the running time of multivariate 
 * methods (MANTA, MANOVA and GEMMA) in the context of multi-trait 
 * GWAS studies.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


/*
 *  Define parameters
 */

params.genotype = null
params.metadata = null
params.dir = 'result'
params.out = 'simulationRT.tsv'
params.n = 1000
params.q = 3
params.A = 10
params.p = 10000
params.b = 1000 
params.k = 20
params.hg2 = 0.2
params.seed = 123
params.r = 1
params.t = 1
params.gemma = true
params.manova = true
params.fx = "$baseDir/supp"
params.help = false


/*
 *  Print usage and help
 */

if (params.help) {
    log.info ''
    log.info 'S I M U L A T E R T - N F'
    log.info '======================================================================='
    log.info 'Benchmark running time and RAM usage of MANTA, MANOVA and GEMMA'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run simulateRT.nf [options]'
    log.info ''
    log.info 'Parameters:'
    log.info ' --genotype GENOTYPES        genotype VCF file from 1000G Phase 3 no duplicates'
    log.info ' --metadata METADATA         metadata from 1000G Phase 3'
    log.info ' --n INDIVIDUALS             number of individuals (default: 1000)'
    log.info ' --q RESPONSES               number of response variables (default: 3)'
    log.info ' --p VARIANTS                number of variants to test (default: 10000)'
    log.info ' --b BLOCKSIZE               variants per block (default: 1000)'
    log.info ' --A ANCESTORS               number of ancestors (default: 10)'
    log.info ' --k NUMBER PC               number of PCs used to correct population stratification (default: 20)'
    log.info ' --hg2 REL HERITABILITY      average fraction of variance explained by relatedness across traits (default: 0)'
    log.info ' --s SEED                    seed (default: 123)'
    log.info ' --r REPLICATE NUMBER        replicate number (default: 1)'
    log.info ' --t OPENBLAS THREADS        OpenBLAS number of threads (default: 1)'
    log.info ' --gemma GEMMA               run GEMMA in addition to MANTA (default: true)'
    log.info ' --manova MANOVA             run MANOVA in addition to MANTA (default: true)'
    log.info ' --fx FUNCTIONS              path to helper functions and precomputed datasets (default: ./supp)'
    log.info ' --dir DIRECTORY             output directory (default: result)'
    log.info ' --out OUTPUT                output file prefix (default: simulationRT.tsv)'
    log.info ''
    exit(1)
}


/*
 *  Print parameter selection
 */

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Genotype data                : ${params.genotype}"
log.info "Metadata                     : ${params.metadata}"
log.info "No. of individuals           : ${params.n}"
log.info "Number of responses          : ${params.q}"
log.info "No. of variants              : ${params.p}"
log.info "No. of variants per block    : ${params.b}"
log.info "No. of ancestors             : ${params.A}"
log.info "No. of PCs                   : ${params.k}"
log.info "Relatedness heritability     : ${params.hg2}"
log.info "Replicate number             : ${params.r}"
log.info "OpenBLAS number of threads   : ${params.t}"
log.info "Seed                         : ${params.seed}"
log.info "Run GEMMA                    : ${params.gemma}"
log.info "Run MANOVA                   : ${params.manova}"
log.info "Helper functions             : ${params.fx}"
log.info "Output directory             : ${params.dir}"
log.info "Output file prefix           : ${params.out}"
log.info ''


/*
 * Checks 
 */

// Mandatory options

if (!params.genotype) {
    exit 1, "Genotype file not specified."
} else if (!params.metadata){
    exit 1, "Metadata file not specified."
}

// Total no. of variants and block size should be compatible

if( params.p % params.b != 0) {
    exit 1, sprintf('Error: %s %% %s != 0', params.p, params.b)
} 

// Check | we iterate either over n or over q
//if (params['n'] =~ /,/ && params['q'] =~ /,/){
//   exit 1, "Provide multiple values EITHER for n or q"
//} 

// Expand n (allow k for thousands)
if (params['n'] =~ /,/){
    Channel.fromList(params['n'].replaceAll("k", "000").tokenize(',')).set{n_ch}
} else {
    Channel.of(params['n'].toString().replaceAll("k", "000")).set{n_ch}
}

// Expand q
if (params['q'] =~ /,/){
    Channel.fromList(params['q'].tokenize(',')).set{q_ch}
} else {
    Channel.of(params['q']).set{q_ch}
}


/* 
 *  Subset VCF
 */

process subset {

    input:
    file(vcf) from file(params.genotype)
    val(p) from params.p

    output:
    file("ss.vcf") into ssvcf_ch

    script:
    """
    # Thin (subset p)
    plink2 --vcf $vcf --make-bed --out geno --threads 1
    plink2 --bfile geno --seed ${params.seed} --thin-count $params.p --make-bed --out geno.thin --threads 1
    plink2 --bfile geno.thin --export vcf --out ss --threads 1
    """
}


/*
 *  Generate ancestors
 */

process generate {

    tag { "n:$n" }

    input:
    file(meta) from file(params.metadata)
    each n from n_ch

    output:
    tuple val(n), file('ancestors.pickle') into pickle_ch

    script:
    """
    generate_ancestors.py -m $meta -A ${params.A} -n $n -s ${params.seed} -M p -p ancestors.pickle
    """
}


/*
 *  Simulate genotypes and compute PCA
 */

process simulateGT {

    tag { "n:$n,r:$r" }

    input:
    each r from Channel.from(1..params.r)
    file(vcf) from ssvcf_ch
    tuple val(n), file(pickle) from pickle_ch

    output:
    tuple val(n), val(r), val("geno"), file ("geno.bed"), file("geno.bim"), file("geno.fam"), file("geno.eigenvec") into gt_ch
    file("runtime.pca.txt") into runtime_pca_ch

    script:
    """
    simulate.py -g $vcf -A ${params.A} -p $pickle -n $n -b ${params.b} -s $r -o sim.vcf --named

    # convert to PLINK format
    plink2 --vcf sim.vcf --make-bed --out geno --threads 1

    # Prune and compute PCA (assumed biallelic variants, MAF >0.05, no missing genotypes)
    start=\$(date +%s)
    plink2 --bfile geno --indep-pairwise 50 5 0.1 --out geno --threads 1
    plink2 --bfile geno --extract geno.prune.in --out geno.pruned --make-bed --threads 1
    if [[ $n -ge 5000 ]]; then approx="approx"; else approx=""; fi
    plink2 --bfile geno.pruned --pca ${params.k} \$approx --out geno --threads 1
    end=\$(date +%s)
    touch runtime.pca.txt
    for q in {${params.q},}; do
        for method in {MANTA,MANOVA}; do
            echo -e "$n\t\$q\t$r\t\$method\tpca\t\$((end-start))" >> runtime.pca.txt
        done
    done
    """
}


/*  
 *  Compute kinship
 */

process kinship {

    tag { "n:$n,r:$r" }

    label "high_mem"

    input:
    tuple val(n), val(r), val(prefix), file(bed),file(bim),file(fam),file(pcs) from gt_ch
 
    output:
    tuple val(n), val(r), val(prefix), file(bed),file(bim),file(fam),file(pcs),file("kinship.sXX.txt.gz") into gt2pt_ch
    file('runtime.kinship.txt') into runtime_kinship_ch

    script:
    """
    # Compute kinship
    sed -i 's/-9/1/' $fam    
    export OPENBLAS_NUM_THREADS=${params.t}
    start=\$(date +%s)
    gemma -gk 2 -bfile $prefix -outdir . -o kinship
    gzip kinship.sXX.txt
    end=\$(date +%s)
    for q in {${params.q},}; do
        echo -e "$n\t\$q\t$r\tGEMMA\tkinship\t\$((end-start))" >> runtime.kinship.txt
    done 
    """
}


/*
 *  Simulate phenotype
 */ 

process simulatePT {

    tag { "n:$n,q:$q,r:$r" }

    label 'high_mem'

    input:
    each q from q_ch
    tuple val(n), val(r), val(prefix), file(bed),file(bim),file(fam),file(pcs),file(kinship) from gt2pt_ch
   
    output:
    tuple val(n), val(q), val(r), val(prefix), file(bed),file(bim),file(fam),file(pcs),file(kinship),file('pheno.txt') into totime_ch

    script:
    """ 
    export R_DATATABLE_NUM_THREADS=1
    # Simulate phenotype
    simulatePT.R -s $r -n $n -q $q --geno $prefix --kinship $kinship --hg2 ${params.hg2} --fx ${params.fx} -o pheno.txt
    """
}


/*
 *  Timing
 */

if (params.gemma & params.manova){
    methods_ch = Channel.fromList(["MANTA", "GEMMA", "MANOVA"]) 
} else if (!params.gemma & params.manova) {
    methods_ch = Channel.fromList(["MANTA", "MANOVA"])
} else if (params.gemma & !params.manova) {
    methods_ch = Channel.fromList(["MANTA", "GEMMA"])
} else if (!params.gemma & !params.manova) {
    methods_ch = Channel.fromList(["MANTA"])
} 

process time {

    tag { "$method n:$n,q:$q,r:$r" }
  
    label 'high_mem'

    input:
    tuple val(n), val(q), val(r), val(prefix), file(bed),file(bim),file(fam),file(pcs),file(kinship),file(pheno) from totime_ch
    each method from methods_ch
   
    output:
    file("runtime.txt") into out_ch

    script:
    pids = (1..q.toInteger()).join(' ')
    if(method == "GEMMA"){
    """ 
    paste <(cut -f1-5 $fam) $pheno > tmpfile; mv tmpfile $fam
    export OPENBLAS_NUM_THREADS=${params.t}
    (timeout 600 gemma -lmm -b $prefix -k $kinship -n $pids -outdir . -o gemma &> STATUS || exit 0)
    if [[ \$(grep ERROR STATUS) ]]; then
        echo -e "$n\t$q\t$r\t$method\tgemma\tNA)" > runtime.txt
    else
        start=\$(date +%s)
        gemma -lmm -b $prefix -k $kinship -n $pids -outdir . -o gemma
        end=\$(date +%s)
        echo -e "$n\t$q\t$r\t$method\tgemma\t\$((end-start))" > runtime.txt
    fi
    """
    } else if (method == "MANTA"){
    """
    export R_DATATABLE_NUM_THREADS=1
    manta.R -p $pheno -g $prefix -c $pcs -k ${params.k} --manta manta.assoc.txt --runtime -i $r > runtime.txt
    """
    } else if (method == "MANOVA"){
    """
    export R_DATATABLE_NUM_THREADS=1
    manta.R -p $pheno -g $prefix -c $pcs -k ${params.k} --manova manova.assoc.txt --runtime -i $r > runtime.txt
    """
    }
}

runtime_kinship_ch.concat(runtime_pca_ch).concat(out_ch).collectFile(name: "${params.out}", sort: { it.text }).set{pub_ch}

process end {

   publishDir "${params.dir}", mode: 'copy'

   input:
   file(runtimes) from pub_ch

   output:
   file(runtimes) into end_ch

   script:
   """
   sed -i "1 s/^/n\tq\tr\tmethod\tstep\truntime\\n/" $runtimes
   """
}

