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
params.n = 10000
params.A = 10
params.l = 10000
params.b = 1000 
params.v = 100000
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
  log.info '    nextflow run simulate-gt.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info ' --genotype GENOTYPES        genotype VCF file from 1000G Phase 3 no duplicates (default: genotypes.vcf.gz)'
  log.info ' --metadata METADATA         metadata from 1000G Phase 3 (default: metadata.tsv)'
  log.info ' --n INDIVIDUALS             number of individuals (default: 10,000)'
  log.info ' --v VARIANTS                number of variants (default: 100,000)'
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
log.info "No. of variants (total)      : ${params.v}"
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
    bcftools query -f '%CHROM\t%POS\n' $vcf | shuf --random-source <(seed ${params.seed}) -n ${params.v} | sort -V > positions
    split -d -a 10 -l ${params.l} positions chunk
    """
}


/*
 *  Simulate individuals
 */

process simulate {

    echo true
    conda '/nfs/users2/rg/dgarrido/.conda/envs/ml'

    input:
    file(vcf) from file(params.genotype)
    file(index) from file("${params.genotype}.tbi")
    file(meta) from file(params.metadata) 
    each file(chunk) from chunks_ch

    output:
    file ("${chunk}.sim") into sim_ch
    file ("${chunk}.sim.vcf") optional true into simvcf_ch

    script:
    if (params.pca == true)
    """
    s=\$(echo $chunk | sed -r 's,chunk0*(.+),\\1,')  # new seed is chunk id
    bcftools view -R $chunk -T $chunk -Ob $vcf | bcftools norm -d all -Ov -o ${chunk}.vcf 
    simulate.py -g ${chunk}.vcf -m $meta -A ${params.A} -n ${params.n} -b ${params.b} -s \$s -M ${params.mode} -o ${chunk}.sim -v ${chunk}.sim.vcf
    """
    else
    """
    s=\$(echo $chunk | sed -r 's,chunk0*(.+),\\1,')  # new seed is chunk id
    bcftools view -R $chunk -T $chunk -Ob $vcf | bcftools norm -d all -Ov -o ${chunk}.vcf
    simulate.py -g ${chunk}.vcf -m $meta -A ${params.A} -n ${params.n} -b ${params.b} -s \$s -M ${params.mode} -o ${chunk}.sim
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
        plink --vcf $simvcf --pca --out ${params.out}
        """
    }

}
