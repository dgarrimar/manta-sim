/*
 * Simulation setting to benchmark multivariate methods using real genotype data
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.GTdir = "GT/1k_chr22"
params.dir = 'result'
params.out = 'simulation.tsv'
params.help = false

// Simulation params
params.n = 100 
params.q = 3
params.p = 1000
params.r = 1
params.PTgen = 'matrixNorm'    
params.GTgen = 'simPopStructure_chr22' 
params.level = 0.05
params.m = 'none'
params.s = 0
params.hs2 = 0.01 
params.hg2 = 0.3
params.alphaG = 0.5
params.lambda = 0.5
params.alphaH = 0.5
params.t = 'none'

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
  log.info ' --p VARIANTS                number of variants to test (default: 1000)'
  log.info ' --r REPEAT                  times to repeat each simulation (default: 1)'
  log.info ' --PTgen PHENO GENERATION    phenotype data generation: matrixNorm or copula (default: matrixNorm)'
  log.info ' --GTgen GENO GENERATION     genotype data generation: simPopStructure, simUnrelated, simRelated, simEmpirical (default: simPopStructure)'
  log.info ' --GTdir GT DIR              input genotype directory with files generated by simulateGT.nf (default: data)' 
  log.info ' --level SIG LEVEL           significance level (default: 0.05)'
  log.info ' --m MULTIPLE TESTING        perform multiple testing correction: none or any method in R:p.adjust (default: none)'
  log.info ' --s CAUSAL VARIANTS         number of causal variants (default: 0)'
  log.info ' --hs2 SNP HERITABILITY      average fraction of variance explained by causal variants across traits (default: 0.01)'
  log.info ' --hg2 REL HERITABILITY      average fraction of variance explained by relatedness across traits (default: 0.3)'
  log.info ' --alphaG REL SHARED         fraction of signal from the relatedness contribution shared across traits (default: 0.5)'
  log.info ' --lambda STRUCT NOISE       fraction of structured noise (default: 0.5)'
  log.info ' --alphaH S. NOISE SHARED    fraction of structured noise that is shared across traits (default: 0.5)'
  log.info ' --t TRANSFORM               transformation of response variables: none, PCA, GAMMA (default: none)'
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
log.info "Variants to test             : ${params.p}"
log.info "Repeats                      : ${params.r}"
log.info "Phenotype data generation    : ${params.PTgen}"
log.info "Genotype data generation     : ${params.GTgen}"
log.info "Genotype input directory     : ${params.GTdir}"
log.info "Significance level           : ${params.level}"
log.info "Multiple testing correction  : ${params.m}"
log.info "Number of causal variants    : ${params.s}"
log.info "Causal variant heritability  : ${params.hs2}"
log.info "Relatedness heritability     : ${params.hg2}"
log.info "Fraction of hg2 shared       : ${params.alphaG}"
log.info "Fraction of st. noise        : ${params.lambda}"
log.info "Fraction of st. noise shared : ${params.alphaH}"
log.info "Transformation               : ${params.t}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''

/*
 *  Expand parameters
 */

def grid = [:]
params.keySet().each{
  if(it in ['n','q','PTgen','GTgen','s','hs2','hg2','alphaG','lambda','alphaH','t','m']){
    grid[it] = params[it]
  }
}

def grid2ch = grid.clone() 
grid.keySet().each {
    if (grid[it] =~ /,/){
        grid[it] = grid[it].tokenize(',')
        grid2ch[it] = Channel.fromList(grid[it].clone())
    } else if (grid[it] =~ /:/) {
        def (start, end, step) = grid[it].tokenize(':')
        def seq = []
        def val = start.toFloat()
        while(val.round(4) <= end.toFloat().round(4)){
            seq += val.round(4)
            val += step.toFloat().round(4)
        }
        grid[it] = seq
        grid2ch[it] = Channel.fromList(seq)
    } else {
        grid2ch[it] = Channel.of(grid[it])
    }
}

grid2ch.n
.combine(grid2ch.q)
.combine(grid2ch.PTgen)
.combine(grid2ch.GTgen)
.combine(grid2ch.s)
.combine(grid2ch.hs2)
.combine(grid2ch.hg2)
.combine(grid2ch.alphaG)
.combine(grid2ch.lambda)
.combine(grid2ch.alphaH)
.combine(Channel.fromList(1..params.r))
.map{

    id = [it[3] + "|n=" + it[0]]
    it = id + it

}.set{grid_ch}

/*
 *  Compute kinship
 */

process kinship {

    tag { id }

    input:
    each n from grid.n
    each GT from Channel.from(grid.GTgen).map { ["${it}", file("${params.GTdir}/${it}.gemma")] }
    
    output:
    tuple val(id), file("${GT[0]}.n${n}.p${params.p}.gemma"), file("${GT[0]}.sXX.txt") into kinship_ch

    script:
    (GTgen, gemma) = GT
    id = "$GTgen|n=$n"
    """
    # Subset genotype by n
    if [[ \$(( $n + 3 )) -lt \$(head -1 $gemma | awk 'NR==1{print NF}') ]]; then
        cut -d ',' -f 1-\$(( $n + 3 )) $gemma > ${GTgen}.n${n}.gemma
    else 
        ln -s ${gemma} ${GTgen}.n${n}.gemma
    fi
    
    # Compute kinship using all variants
    Rscript -e 'write.table(file = "dummypheno", rnorm($n), col.names = F, row.names = F)'
    gemma -gk 2 -g ${GTgen}.n${n}.gemma -p dummypheno -outdir . -o $GTgen
    
    # Subset genotype by p
    head -n ${params.p} ${GTgen}.n${n}.gemma > ${GTgen}.n${n}.p${params.p}.gemma
    """
}

grid_ch.combine(kinship_ch, by: 0).set{gt2pt_ch}

/*
 *  Simulate phenotype
 */ 

process simulatePT {

    tag { par }

    input:
    tuple id,n,q,PTgen,GTgen,s,hs2,hg2,alphaG,lambda,alphaH,r,file(geno),file(kinship) from gt2pt_ch

    output:
    tuple id, par, file(geno), file(kinship), file('pheno.txt'), file('ids.txt') into pheno_ch

    script:
    par = "$n|$q|$PTgen|$GTgen|$s|$hs2|$hg2|$alphaG|$lambda|$alphaH|$r"
    """ 
    simulatePT.R -r $r -n $n -q $q --PTgen $PTgen --geno $geno --kinship $kinship -s $s --hs2 $hs2 --hg2 $hg2 --alphaG $alphaG --lambda $lambda --alphaH $alphaH -o pheno.txt -i ids.txt
    """
}

if("PCA" in grid.t) {

 /*
  *  Genotype PCA
  */
    
    process pca {
 
        tag { id }    
 
        input:
        each n from grid.n
        each GT from Channel.from(grid.GTgen).map { ["${it}", file("${params.GTdir}/${it}.vcf")] }

        output:
        tuple id, file("${GT[0]}.eigenvec") into pca_ch

        script:
        (GTgen, vcf) = GT
        id = "$GTgen|n=$n"
        """
        # Compute PCs using all variants
        for i in {1..$n}; do echo -ne "S\$i\n" ; done > keep.txt  # subset n individuals
        plink2 --vcf $vcf --out ${GTgen}
        plink2 --pfile ${GTgen} --indep-pairwise 50 5 0.1 --keep keep.txt --out ${GTgen}
        plink2 --pfile ${GTgen} --extract ${GTgen}.prune.in --keep keep.txt --out ${GTgen}.pruned --make-pgen
        if [[ $n -ge 5000 ]]; then approx="approx"; else approx=""; fi
        plink2 --pfile ${GTgen}.pruned --pca \$approx --out ${GTgen} 
        """
    }

    pheno_ch.combine(pca_ch, by: 0).into{input_gemma_ch;input_mlm_ch;input_manova_ch}

} else {
   
    pheno_ch.spread(Channel.of("dummy")).into{input_gemma_ch;input_mlm_ch;input_manova_ch}

} 


/*
 * Run GEMMA, MLM and MANOVA
 */

process GEMMA {
  
    tag { par }

    input:
    tuple dummy, par, file(geno), file(kinship), file(pheno), file(ids), file(eigenval) from input_gemma_ch

    output:
    tuple par_ext, file(ids), file("gemma.assoc.txt") into gemma_ch

    script:
    q = par.split("\\|")[1].toInteger() // # 0-indexed
    pids = (1..q).join(' ')
    par_ext = par + "|GEMMA" 
    """
    gemma -lmm -g $geno -k $kinship -p $pheno -n $pids -outdir . -o gemma
    """
}

process MLM {

    input:
    tuple dummy, par, file(geno), file(kinship), file(pheno), file(ids), file(eigenval) from input_mlm_ch
    each t from grid.t

    output:
    tuple par_ext, file(ids), file("mlm.assoc.txt") into mlm_ch

    script:
    if (t == 'none'){
    par_ext = par + "|MLM"
    """
    mlm.R -p $pheno -g $geno -t $t -o mlm.assoc.txt
    """
    }else if (t == 'PCA'){
    par_ext = par + "|MLM_PCA"
    """
    mlm.R -p $pheno -g $geno -t $t -c $eigenval -o mlm.assoc.txt
    """
    }else if (t == 'GAMMA'){
    par_ext = par + "|MLM_GAMMA"
    """
    for i in {1..$q}; do
        gemma -vc 2 -p $pheno -k $kinship -n \$i -outdir . -o VC &> /dev/null
        grep -F "sigma2 estimates =" VC.log.txt | cut -d ' ' -f 7,9
    done > VC.txt
    mlm.R -p $pheno -g $geno -t $t -k $kinship -v VC.txt -o mlm.assoc.txt
    """
    }
}

process MANOVA { // identical to MLM with --manova

    input:
    tuple dummy, par, file(geno), file(kinship), file(pheno), file(ids), file(eigenval) from input_manova_ch
    each t from grid.t

    output:
    tuple par_ext, file(ids), file("manova.assoc.txt") into manova_ch

    script:
    if (t == 'none'){
    par_ext = par + "|MANOVA"
    """
    mlm.R -p $pheno -g $geno -t $t -o manova.assoc.txt --manova
    """
    } else if (t == 'PCA') {
    par_ext = par + "|MANOVA_PCA"
    """
    mlm.R -p $pheno -g $geno -t $t -o manova.assoc.txt -c $eigenval --manova
    """
    } else if (t == 'GAMMA') {
    par_ext = par + "|MANOVA_GAMMA"
    """
    for i in {1..$q}; do
        gemma -vc 2 -p $pheno -k $kinship -n \$i -outdir . -o VC &> /dev/null
        grep -F "sigma2 estimates =" VC.log.txt | cut -d ' ' -f 7,9
    done > VC.txt
    mlm.R -p $pheno -g $geno -t $t -k $kinship -v VC.txt -o manova.assoc.txt --manova
    """
    }
}

gemma_ch.concat(mlm_ch, manova_ch).set{tie_power_ch}

/*
 *  Compute TIE or power
 */

process tie_power {
    
    input:
    each m from grid.m
    tuple par, file(ids), file(assoc) from tie_power_ch
 
    output:
    file("res.txt") into out_ch  

    script:
    """
    tie.R -t $assoc -l ${params.level} -i $ids -m $m -o tieORpower.txt
    paste <(echo -e "$par|$m" | sed 's/|/\t/g') tieORpower.txt > res.txt
    """
}


out_ch.collectFile(name: "${params.out}", sort: { it.text }).set{pub_ch}

process end {

   publishDir "${params.dir}"

   input:
   file(sim) from pub_ch

   output:
   file(sim) into end_ch

   script:
   """
   sed -i "1 s/^/n\tq\tPTgen\tGTgen\ts\ths2\thg2\talphaG\tlambda\talphaH\tr\tmethod\tmtc\ttie\\n/" $sim
   """
}

