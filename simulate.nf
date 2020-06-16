/*
 * Simulation setting to benchmark multivariate methods using real genotype data
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.GTdir = "$baseDir/simGT"
params.dir = 'result'
params.out = 'simulation.tsv'
params.help = false

// Simulation params
params.n = 1000 // < simulateGT.nf
params.q = 3
params.PTgen = 'matrixNorm'    
params.GTgen = 'simPopStructure' 
params.level = 0.05
params.s = 0
params.hs2 = 0.01 
params.hg2 = 0.3
params.alphaG = 0.5
params.lambda = 0.5
params.alphaH = 0.5
params.t = 'none'
params.r = false

/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'Benchmark of GEMMA, GAMMA and MLM in a simulation using real GT data'
  log.info '======================================================================='
  log.info 'Simulations given a set of parameters'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run simulate.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info ' --n SAMPLE SIZE             total number of samples (default: 100)'
  log.info ' --q RESPONSES               number of response variables (default: 3)'
  log.info ' --PTgen PHENO GENERATION    phenotype data generation: matrixNorm or copula (default: matrixNorm)'
  log.info ' --GTgen GENO GENERATION     genotype data generation: simPopStructure, simUnrelated, simRelated, simEmpirical (default: simPopStructure)'
  log.info ' --GTdir GT DIR              input genotype directory with files generated by simulateGT.nf (default: data)' 
  log.info ' --level SIG LEVEL           significance level (default: 0.05)'
  log.info ' --s CAUSAL VARIANTS         number of causal variants (default: 0)'
  log.info ' --hs2 SNP HERITABILITY      average fraction of variance explained by causal variants across traits (default: 0.01)'
  log.info ' --hg2 REL HERITABILITY      average fraction of variance explained by relatedness across traits (default: 0.3)'
  log.info ' --alphaG REL SHARED         fraction of signal from the relatedness contribution shared across traits (default: 0.5)'
  log.info ' --lambda STRUCT NOISE       fraction of structured noise (default: 0.5)'
  log.info ' --alphaH S. NOISE SHARED    fraction of structured noise that is shared across traits (default: 0.5)'
  log.info ' --t TRANSFORM               Transformation of response variables: none, PCA, GAMMA (default: none)'
  log.info ' --r RUNNING TIME            report running time rather than tie/power (default: false)'
  log.info ' --dir OUTPUT DIR            output directory (default: result)'
  log.info ' --out OUTPUT                output file (default: simulation.tsv)'
  log.info ''
  exit(1)
}

/*
 *  Print parameter selection
 */

log.info ''
log.info 'SIMULATION SCHEMA'
log.info ''
log.info 'General parameters'
log.info '------------------'
log.info "Total sample size            : ${params.n}"
log.info "Number of responses          : ${params.q}"
log.info "Phenotype data generation    : ${params.PTgen}"
log.info "Genotype data generation     : ${params.GTgen}"
log.info "Genotype input directory     : ${params.GTdir}"
log.info "Significance level           : ${params.level}"
log.info "Number of causal variants    : ${params.s}"
log.info "Causal variant heritability  : ${params.hs2}"
log.info "Relatedness heritability     : ${params.hg2}"
log.info "Fraction of hg2 shared       : ${params.alphaG}"
log.info "Fraction of st. noise        : ${params.lambda}"
log.info "Fraction of st. noise shared : ${params.alphaH}"
log.info "Transformation               : ${params.t}"
log.info "Running time                 : ${params.r}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''

/*
 *  Expand parameters
 */

def grid = [:]
params.keySet().each{
  if(it in ['n','q','PTgen','GTgen','s','hs2','hg2','alphaG','lambda','alphaH','t','pca']){
    grid[it] = params[it]
  }
}

grid.keySet().each {
    if (grid[it] =~ /,/){
        grid[it] = grid[it].tokenize(',')
    } else if (grid[it] =~ /:/) {
        def (start, end, step) = grid[it].tokenize(':')
        def seq = []
        def val = start.toFloat()
        while(val.round(4) <= end.toFloat().round(4)){
            seq += val.round(4)
            val += step.toFloat().round(4)
        }
        grid[it] = seq
    }
}

/*
 *  Simulate phenotype
 */

process simulatePT {
 
    input:
    each n from grid.n
    each q from grid.q
    each PTgen from grid.PTgen
    each GT from Channel.from(grid.GTgen).map { ["${it}", file("${params.GTdir}/${it}.gemma"), file("${params.GTdir}/${it}.sXX.txt"), file("${params.GTdir}/${it}.eigenvec")] } 
    each s from grid.s
    each hs2 from grid.hs2
    each hg2 from grid.hg2
    each alphaG from grid.alphaG
    each lambda from grid.lambda
    each alphaH from grid.alphaH

    output:    
    set val(q), GT, file('pheno.txt'), file('params.txt') into pheno_ch

    script:
    def (GTgen, geno, kinship, eigenvec) = GT
    """
    echo -e "$n\t$q\t$PTgen\t$GTgen\t$s\t$hs2\t$hg2\t$alphaG\t$lambda\t$alphaH" > params.txt
    # Generate new GT and Kinship given n if n < n_geno, use it downstream
    simulatePT.R -n $n -q $q --PTgen $PTgen --geno $geno --kinship $kinship -s $s --hs2 $hs2 --hg2 $hg2 --alphaG $alphaG --lambda $lambda --alphaH $alphaH -o pheno.txt
    """
}

pheno_ch.flatten().toList().into{pheno1_ch; pheno2_ch}

process GEMMA {

    input:
    set val(q), val(GTgen), file(geno), file(kinship), file(eigenval), file(pheno), file(par) from pheno1_ch
    
    output:
    file("gemma.txt") into gemma_ch

    script:
    def pids = (1..3).join(' ')
    """
    head -n 100 $geno > fastgeno
    gemma -lmm -g fastgeno -k $kinship -p $pheno -n $pids -outdir . -o $GTgen
    paste $par <(echo "GEMMA") <(Rscript -e 'tb <- data.table::fread("${GTgen}.assoc.txt", data.table = F); cat(mean(tb[, ncol(tb)] < ${params.level}))') > gemma.txt
    """
}

process MLM {

    input:
    set val(q), val(GTgen), file(geno), file(kinship), file(eigenval), file(pheno), file(par) from pheno2_ch
    each t from grid.t
  
    output:
    file("mlm.txt") into mlm_ch

    script:
    if (t == 'none')
    """
    mlm=\$(mlm.R -p $pheno -g $geno -t $t -o ${GTgen}.mlm)
    paste $par <(echo "MLM") <(echo \$mlm) > mlm.txt     
    """
    else if (t == 'PCA')
    """
    mlm=\$(mlm.R -p $pheno -g $geno -t $t -o ${GTgen}.mlm -c $eigenval   
    paste $par <(echo "MLM") <(echo \$mlm) > mlm.txt
    """
    else if (t == 'GAMMA')
    """
    for in {1..$q}; do 
        gemma -vc 2 -p $pheno -k $kinship -n \$i -outdir . -o VC &> /dev/null 
        grep -F "sigma2 estimates =" VC.log.txt | cut -d ' ' -f 7,9
    done > VC.txt 
    mlm=\$(mlm.R -p $pheno -g $geno -t $t -k $kinship -v VC.txt -o ${GTgen}.mlm)
    paste $par <(echo "MLM_GAMMA") <(echo \$mlm) > mlm.txt
    """
}
