/*
 * Simulate genotypes as in P. Casale (thesis)
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.genotype = 'data/genotypes.vcf.gz'
params.metadata = 'data/metadata.tsv'
params.dir = 'result'
params.out = 'simulateGT'
params.n = 1000
params.A = 10
params.l = 10000
params.b = 1000 
params.seed = 123
params.pca = false
params.mode = 'e' // e, p, ep
params.help = false

/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'SIMULATE GT'
  log.info '======================================================================='
  log.info 'Simulate genotypes with different population structure and relatedness'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run simulateGT.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info ' --genotype GENOTYPES        genotype VCF file from 1000G Phase 3 no duplicates (default: data/genotypes.vcf.gz)'
  log.info ' --metadata METADATA         metadata from 1000G Phase 3 (default: data/metadata.tsv)'
  log.info ' --n INDIVIDUALS             number of individuals (default: 10,000)'
  log.info ' --l VARIANTS/CHUNK          variants per chunk (default: 10,000)'
  log.info ' --b BLOCKSIZE               variants per block (default: 1000)'
  log.info ' --A ANCESTORS               number of ancestors (default: 10)'
  log.info ' --M MODE                    simulation mode (default: e)'
  log.info ' --pca PCA                   perform PCA (default: false)'
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
log.info "No. of variants per chunk    : ${params.l}"
log.info "No. of variants per block    : ${params.b}"
log.info "No. of ancestors             : ${params.A}"
log.info "Perform PCA                  : ${params.pca}"
log.info "Simulation mode              : ${params.mode}"
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

// Chunk and block sizes should be compatible

if( params.l % params.b != 0) {
    exit 1, sprintf('Error: %s %% %s != 0', params.l, params.b)
} 


/* 
 *  Split
 */

process split {

    input:
    file(vcf) from file(params.genotype)
    file(index) from file("${params.genotype}.tbi")

    output:
    file ("chunk*") into chunks_ch

    script:
    """
    bcftools query -f '%CHROM\t%POS\n' $vcf > positions
    split -d -a 10 -l ${params.l} positions chunk
    """
}


/*
 *  Generate ancestors
 */

process generate {

   input:
   file(meta) from file(params.metadata)
   
   output:
   file('ancestors.pickle') into pickle_ch

   script:
   """
   generate_ancestors.py -m $meta -A ${params.A} -n ${params.n} -s ${params.seed} -M ${params.mode} -p ancestors.pickle
   """
}

/*
 *  Simulate individuals
 */

process simulate {

    input:
    file(vcf) from file(params.genotype)
    file(index) from file("${params.genotype}.tbi")
    file(pickle) from pickle_ch
    each file(chunk) from chunks_ch

    output:
    file ("${chunk}.sim") into sim_ch
    file ("${chunk}.sim.vcf") optional true into simvcf_ch

    script:
    if (params.pca == true)
    """
    s=\$(echo $chunk | sed -r 's,chunk0*(.+),\\1,')  # seed is chunk id
    bcftools view -R $chunk -T $chunk -Ob $vcf | bcftools norm -d all -Ov -o ${chunk}.vcf 
    simulate.py -g ${chunk}.vcf -A ${params.A} -p $pickle -n ${params.n} -b ${params.b} -s \$s -o ${chunk}.sim -v ${chunk}.sim.vcf
    """
    else
    """
    s=\$(echo $chunk | sed -r 's,chunk0*(.+),\\1,')  # seed is chunk id
    bcftools view -R $chunk -T $chunk -Ob $vcf | bcftools norm -d all -Ov -o ${chunk}.vcf
    simulate.py -g ${chunk}.vcf -A ${params.A} -p $pickle -n ${params.n} -b ${params.b} -s \$s -o ${chunk}.sim
    """
}


/*
 *  Reduce and generate output
 */

sim_ch.collectFile(name: "${params.out}.gemma", sort: { it.name }).set{out_ch}
if ( params.pca ){ simvcf_ch.collectFile(name: "${params.out}.vcf", sort: { it.name }).set{outvcf_ch} }

process out {

    publishDir "${params.dir}"

    input:
    file(sim) from out_ch
 
    output:
    file(sim) into pub_ch
 
    script: 
    """
    echo "Done!"
    """
}

/*
 * Plot PCA
 *
 *  - We assume there are not missing genotypes
 *  - We assume biallelic variants
 *  - We assume MAF > 0.05
 *
 */

if (params.pca) {

    process pca {

        publishDir "${params.dir}"
    
        input:
        file(simvcf) from outvcf_ch
    
        output:
        file("${params.out}.eigen*") into pca_ch

        script:
        """
        ids=\$(for (( i = 1; i <= $params.n; i++ )); do echo -ne "S\$i\t" ; done | sed 's,\t\$,,')
        sed -i "1 s,^,#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\$ids\\n," $simvcf        
        
        # 1) convert to PLINK2 format
        plink2 --vcf $simvcf --out ${params.out}
        
        # 2) prune with --indep-pairwise
        plink2 --pfile ${params.out} --indep-pairwise 50 5 0.1 --out ${params.out}
        plink2 --pfile ${params.out} --extract ${params.out}.prune.in --out ${params.out}.pruned --make-pgen
        
        # 4) Compute PCA
        if [[ ${params.n} -ge 5000 ]]; then approx="approx"; else approx=""; fi
        plink2 --pfile ${params.out}.pruned --pca \$approx --out ${params.out}
        """
    }

}
