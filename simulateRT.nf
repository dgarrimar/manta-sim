/*
 * Simulation setting to benchmark multivariate methods using real genotype data (running time)
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.genotype = 'data/genotypes.vcf.gz'
params.metadata = 'data/metadata.tsv'
params.dir = 'result'
params.out = 'simulateRT'
params.n = 1000
params.q = 3
params.A = 10
params.p = 10000
params.b = 1000 
params.k = 20
params.hs2 = 0 
params.hg2 = 0.2
params.seed = 123
params.help = false

/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'SIMULATE RT'
  log.info '======================================================================='
  log.info 'Benchmark running time of GEMMA, GAMMA and MLM in a simulation using real GT data'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run simulateRT.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info ' --genotype GENOTYPES        genotype VCF file from 1000G Phase 3 no duplicates (default: data/genotypes.vcf.gz)'
  log.info ' --metadata METADATA         metadata from 1000G Phase 3 (default: data/metadata.tsv)'
  log.info ' --n INDIVIDUALS             number of individuals (default: 1000)'
  log.info ' --q RESPONSES               number of response variables (default: 3)'
  log.info ' --p VARIANTS                number of variants to test (default: 10000)'
  log.info ' --b BLOCKSIZE               variants per block (default: 1000)'
  log.info ' --A ANCESTORS               number of ancestors (default: 10)'
  log.info ' --k NUMBER PC               Number of PCs used to correct population stratification (default: 20)'
  log.info ' --hs2 SNP HERITABILITY      average fraction of variance explained by causal variants across traits (default: 0)'
  log.info ' --hg2 REL HERITABILITY      average fraction of variance explained by relatedness across traits (default: 0)'

  log.info ' --s SEED                    seed (default: 123)'
  log.info ' --dir DIRECTORY             output directory (default: result)'
  log.info ' --out OUTPUT                output file prefix (default: simulated)'
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
log.info "Causal variant heritability  : ${params.hs2}"
log.info "Relatedness heritability     : ${params.hg2}"
log.info "Seed                         : ${params.seed}"
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

// Expand n (allow k for thousands)
if (params['n'] =~ /,/){
    Channel.fromList(params['n'].replaceAll("k", "000").tokenize(',')).set{n_ch}
} else {
    Channel.of(params['n'].toString().replaceAll("k", "000")).set{n_ch}
}


/* 
 *  Subset
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
    plink2 --vcf $vcf --make-bed --out geno
    plink2 --bfile geno --seed ${params.seed} --thin-count $params.p --make-bed --out geno.thin
    plink2 --bfile geno.thin --export vcf --out ss     
    """
}

/*
 *  Generate ancestors
 */

process generate {

    tag { n }

    label 'high_mem'

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

    tag { n }

    input:
    file(vcf) from ssvcf_ch
    tuple val(n),file(pickle) from pickle_ch

    output:
    tuple val(n), val("geno"), file ("geno.bed"), file("geno.bim"), file("geno.fam"), file("geno.eigenvec") into gt_ch
    file("runtime.pca.txt") into runtime_pca_ch

    script:
    """
    simulate.py -g $vcf -A ${params.A} -p $pickle -n $n -b ${params.b} -s ${params.seed} -o sim.vcf --named 

    # convert to PLINK format
    plink2 --vcf sim.vcf --make-bed --out geno

    # Prune and compute PCA (assumed biallelic variants, MAF >0.05, no missing genotypes)
    start=\$(date +%s)
    plink2 --bfile geno --indep-pairwise 50 5 0.1 --out geno
    plink2 --bfile geno --extract geno.prune.in --out geno.pruned --make-bed
    if [[ $n -ge 5000 ]]; then approx="approx"; else approx=""; fi
    plink2 --bfile geno.pruned --pca ${params.k} \$approx --out geno   
    end=\$(date +%s)
    echo -e "$n\tPCA\tpca\t\$((end-start))" > runtime.pca.txt
    """
}


/*  
 *  Compute kinship
 */

process kinship {

    tag { n }

    input:
    tuple n, prefix, file(bed),file(bim),file(fam),file(pcs) from gt_ch
 
    output:
    tuple n, prefix, file(bed),file(bim),file(fam),file(pcs),file("kinship.sXX.txt") into gt2pt_ch
    file('runtime.kinship.txt') into runtime_kinship_ch

    script:
    """
    # Compute kinship
    sed -i 's/-9/1/' $fam    
    start=\$(date +%s)
    gemma -gk 2 -bfile $prefix -outdir . -o kinship
    end=\$(date +%s)    
    echo -e "$n\tKINSHIP\tkinship\t\$((end-start))" > runtime.kinship.txt
    """
}

/*
 *  Simulate phenotype
 */ 

process simulatePT {

    tag { n }

    input:
    tuple n, prefix, file(bed),file(bim),file(fam),file(pcs),file(kinship) from gt2pt_ch
   
    output:
    tuple val(n), val(prefix), file(bed),file(bim),file(fam),file(pcs),file(kinship),file('pheno.txt') into totime_ch

    script:
    """ 
    # Simulate phenotype
    simulatePT.R -s ${params.seed} -n $n -q ${params.q} --geno $prefix --kinship $kinship --hs2 ${params.hs2} --hg2 ${params.hg2} -o pheno.txt 
    """
}

totime_ch.set{totime_ch2}

/*
 *  Timing
 */

process time {

    tag { "$t $n" }

    input:
    tuple n, prefix, file(bed),file(bim),file(fam),file(pcs),file(kinship),file(pheno) from totime_ch2
    each t from Channel.fromList(["GEMMA","GAMMA","PCA"])
   
    output:
    file("runtime.txt") into out_ch

    script:
    pids = (1..params.q.toInteger()).join(' ')
    if(t == "GEMMA"){
    """ 
    paste <(cut -f1-5 $fam) $pheno > tmpfile; mv tmpfile $fam
    start=\$(date +%s)
    gemma -lmm -b $prefix -k $kinship -n $pids -outdir . -o gemma.assoc.txt
    end=\$(date +%s)
    echo -e "$n\t$t\tgemma\t\$((end-start))" > runtime.txt
    """
    } else if (t == "GAMMA" && n.toInteger() < 5000) {
    """
    paste <(cut -f1-2 $fam) $pheno > pheno2.txt
    start=\$(date +%s`)
    vc.py -b $prefix -p pheno2.txt -k $kinship -o VC.txt -v
    end=\$(date +%s)
    echo -e "$n\t$t\tvc\t\$((end-start))" > runtime.txt

    mlm.R -p $pheno -g $prefix -t $t -c $pcs -n ${params.k} -k $kinship -v VC.txt --mlm mlm.assoc.txt --runtime >> runtime.txt
    """
    } else if (t == "GAMMA" && n.toInteger() >= 5000) {
    """
    PIDS=\$(echo $pids | sed 's/ /\\t/g')
    paste <(cut -f1-2 $fam) $pheno | sed "1 s/^/FID\\tIID\\t\$PIDS\\n/" > pheno2.txt
    start=\$(date +%s)
    for pid in {1..${params.q}}; do
        bolt --bfile $prefix --phenoFile pheno2.txt --phenoCol \$pid --reml > bolt.txt
        sigma2=\$(grep -F 'Phenotype 1 variance sigma2:' bolt.txt | sed -r 's/.*: (.+) .*/\\1/')
        h2g=\$(grep -F 'h2g (1,1):' bolt.txt | sed -r 's/.*: (.+) .*/\\1/')
        echo -e "\$(echo "\$h2g*\$sigma2" | bc -l)\\t\$(echo "(1-\$h2g)*\$sigma2" | bc -l)"
    done > VC.txt
    end=\$(date +%s)
    echo -e "$n\t$t\tvc\t\$((end-start))" > runtime.txt
    mlm.R -p $pheno -g $prefix -t $t -c $pcs -n ${params.k} -k $kinship -v VC.txt --mlm mlm.assoc.txt --runtime >> runtime.txt
    """
    } else if (t == "PCA") {
    """
    touch VC.txt
    mlm.R -p $pheno -g $prefix -t $t -c $pcs -n ${params.k} -k $kinship -v VC.txt --mlm mlm.assoc.txt --runtime > runtime.txt
    """
    }
}

runtime_kinship_ch.concat(runtime_pca_ch).concat(out_ch).collectFile(name: "${params.out}.txt", sort: { it.text }).set{pub_ch}

process end {

   publishDir "${params.dir}", mode: 'copy'

   input:
   file(runtimes) from pub_ch

   output:
   file(runtimes) into end_ch

   script:
   """
   sed -i "1 s/^/n\tmethod\tstep\truntime\\n/" $runtimes
   """
}

