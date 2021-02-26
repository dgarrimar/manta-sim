/*
 * Simulation setting to benchmark multivariate methods using real genotype data (initial setup)
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.GTdir = "GT/2504_real"
params.dir = 'result'
params.out = 'simulation.tsv'
params.fx = "$baseDir/supp"
params.help = false

// Simulation params
params.n = 100 
params.q = 3
params.p = 100
params.PTgen = 'mvnorm'    
params.GTgen = 'real1000g' 
params.level = 0.05
params.m = 'none'
params.C = 'none'
params.t = 'none'
params.scale = false
params.k = 20
params.c = 10
params.delta = 0
params.hk = 1

// Generation: multivariate normal or copula
params.y_var = 'equal'
params.y_cor = 0

// Generation: copula
params.c_dist = 'unif-0-1'

// Generation: simplex (proportions)
params.p_loc = 1
params.p_sd = 0.1
params.p_dist = 'norm'

// Generation: multinomial
params.lambda = 1000


/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'Benchmark of GEMMA and MLM in a simulation using real GT data (initial approach)'
  log.info '======================================================================='
  log.info 'Simulations given a set of parameters'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run simulate.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info ' --n SAMPLE SIZE             total number of samples (default: 100)'
  log.info ' --q RESPONSES               number of response variables (default: 3)'
  log.info ' --p VARIANTS                number of variants to test (default: 100)'
  log.info ' --PTgen PHENO GENERATION    phenotype data generation: 'mvnorm', 'copula', 'simplex' or 'multinom' (default: mvnorm)'
  log.info ' --GTgen GENO GENERATION     genotype data generation: simPopStructure, simUnrelated, simRelated, simEmpirical (default: real1000g)'
  log.info ' --GTdir GT DIR              input genotype directory with files generated by simulateGT.nf (default: GT/2504_real)' 
  log.info ' --level SIG LEVEL           significance level (default: 0.05)'
  log.info ' --m MULTIPLE TESTING        perform multiple testing correction: none or any method in R:p.adjust (default: none)'
  log.info ' --scale SCALE               [WARN] scale the response variables for MLM (default: false)'
  log.info ' --C COVARIATES              include covariates: none, PCA (default: none)' 
  log.info ' --k NUMBER PC               number of genotype PCs (default: 20)'
  log.info ' --t TRANSFOMATION           response variable transformation: none, sqrt, log (default: none)'
  log.info ' --c NUMBER CHUNKS           Number of chunks (default: 10)'
  log.info ' --dir OUTPUT DIR            output directory (default: result)'
  log.info ' --out OUTPUT                output file (default: simulation.tsv)'
  log.info ' --delta DELTA               change in H1 (default: 0.1)'
  log.info ' --hk HETEROSKEDASTICITY     heteroskedasticity, 1 is homoskedastic (default: 1)'
  log.info ' --fx FUNCTIONS              path to helper functions and precomputed datasets (default: ./supp)'
  log.info ''
  log.info 'Additional parameters for PTgen = mvnorm:'
  log.info ' --y_var VARIANCE            variance of response variables: equal or unequal (default: equal)'
  log.info ' --y_cor CORRELATION         correlation of response variables (default: 0)'
  log.info ''
  log.info 'Additional parameters for PTgen = copula:'
  log.info ' --y_var VARIANCE            variance of response variables: equal or unequal (default: equal)'
  log.info ' --y_cor CORRELATION         correlation of response variables (default: 0)'
  log.info ' --c_dist DISTRIBUTION	 multivariate non-normal distribution definition (default: unif-0-1)'
  log.info ''
  log.info 'Additional parameters for PTgen = simplex:'
  log.info ' --p_loc LOCATION            location in the simplex of the generator model, 1 is centered (default: 1)'
  log.info ' --p_sd STDEV                standard deviation of the generator model (default: 0.1)'
  log.info ' --p_dist DISTRIBUTION	 distribution of the generator model (default: norm)'
  log.info ''
  log.info 'Additional parameters for PTgen = multinom:'
  log.info ' --lambda LAMBDA             lambda parameter (Poisson distribution) to generate multinomial distribution\'s size parameter (default: 1000)'
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
log.info "Phenotype data generation    : ${params.PTgen}"
log.info "Genotype data generation     : ${params.GTgen}"
log.info "Genotype input directory     : ${params.GTdir}"
log.info "Significance level           : ${params.level}"
log.info "Multiple testing correction  : ${params.m}"
log.info "Covariates                   : ${params.C}"
log.info "Transformation               : ${params.t}"
log.info "Scale responses [WARN]       : ${params.scale}"
log.info "Number of PCs                : ${params.k}"
log.info "Number of chunks             : ${params.c}"
log.info "Change in H1 (delta)         : ${params.delta}"
log.info "Heteroskedasticity           : ${params.hk}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info "Helper functions             : ${params.fx}"
log.info ''

if(params.PTgen == 'mvnorm'){
  log.info 'Additional parameters'
  log.info '---------------------'
  log.info "Variance of Y variables      : ${params.y_var}"
  log.info "Correlation of Y variables   : ${params.y_cor}"
} else if (params.PTgen == 'copula'){
  log.info 'Additional parameters'
  log.info '---------------------'
  log.info "Variance of Y variables      : ${params.y_var}"
  log.info "Correlation of Y variables   : ${params.y_cor}"
  log.info "Distribution definition      : ${params.c_dist}"
} else if (params.PTgen == 'simplex'){
  log.info 'Additional parameters'
  log.info '---------------------'
  log.info "Location in the simplex      : ${params.p_loc}"
  log.info "Generation stdev             : ${params.p_sd}"
  log.info "Generation distribution      : ${params.p_dist}"
} else if (params.PTgen == 'multinom'){
  log.info 'Additional parameters'
  log.info '---------------------'
  log.info "Lambda                       : ${params.lambda}"
}
log.info ''


/*
 *  Checks
 */ 

if (params.p%params.c != 0) {
    exit 1, "ERROR: p%c should be 0"
}

/*
 *  Expand parameters
 */

def grid = [:]
params.keySet().each{
  if(it in ['n','q','PTgen','GTgen','hs2','hg2','C','m','delta','hk','y_var','y_cor','c_dist','p_loc','p_sd','p_dist','lambda','t']){
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
        while(val.toFloat().round(4) <= end.toFloat().round(4)){
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
.map{

    id = [it[3] + "|n=" + it[0]]
    it = id + it

}.set{grid_ch}

/*
 *  Prepare genotype 
 */

process prepare {

    tag { id }

    input:
    each n from grid.n
    each GT from Channel.from(grid.GTgen).map { ["${it}", file("${params.GTdir}/${it}")] }

    output:
    tuple val(id), file("${GTgen}.pruned.bed"), file("${GTgen}.pruned.bim"), file("${GTgen}.pruned.fam") into geno_pruned_ch1, geno_pruned_ch2
    tuple val(id), file("${GTgen}.thin.bed"), file("${GTgen}.thin.bim"), file("${GTgen}.thin.fam") into geno_thin_ch

    script:
    (GTgen, prefix) = GT
    id = "$GTgen|n=$n"
    """
    # Subset n and prune
    if [[ $n -lt \$(wc -l ${prefix}.fam | cut -d ' ' -f 1) ]]; then
        head -n $n ${prefix}.fam > keep.txt
    else
        cut -f2 ${prefix}.fam > keep.txt
    fi

    plink2 --bfile $prefix --indep-pairwise 50 5 0.1 --keep keep.txt --out $GTgen --threads 1
    plink2 --bfile $prefix --extract ${GTgen}.prune.in --keep keep.txt --out ${GTgen}.pruned --make-bed --threads 1
    
    # Thin: subset n and p 
    plink2 --bfile $prefix --seed 1 --thin-count $params.p --keep keep.txt --out ${GTgen}.thin --make-bed --threads 1
    """
}

/*  
 *  Compute kinship
 */

process kinship {

    tag { id }

    input:
    tuple val(id), file(bed), file(bim), file(fam) from geno_pruned_ch1
 
    output:
    tuple val(id), file("${GTgen}.sXX.txt.gz") into kinship_ch

    script:
    GTgen = id.split("\\|")[0]
    """
    # Compute kinship
    sed -i 's/-9/1/' $fam    
    gemma -gk 2 -bfile \$(basename $bed | sed 's/.bed//') -outdir . -o $GTgen 
    gzip ${GTgen}.sXX.txt
    """
}


if("PCA" in grid.C) {

   /*
    *  Genotype PCA
    */

    process pca {

        tag { id }

        input:
        tuple val(id), file(bed), file(bim), file(fam) from geno_pruned_ch2

        output:
        tuple val(id), file("${GTgen}.eigenvec") into pca_ch

        script:
        GTgen = id.split("\\|")[0]
        """
        # Compute PCs using all variants
        if [[ \$(wc -l $fam | cut -d' ' -f2) -ge 5000 ]]; then approx="approx"; else approx=""; fi
        plink2 --bfile ${GTgen}.pruned --pca ${params.k} \$approx --out ${GTgen} --threads 1
        """
    }

    grid_ch.combine(geno_thin_ch, by:0).combine(kinship_ch, by: 0).combine(pca_ch, by: 0).set{gt2pt_ch}

} else {

    grid_ch.combine(geno_thin_ch, by:0).combine(kinship_ch, by: 0).spread(Channel.of("")).set{gt2pt_ch}

}

/*
 *  Simulate phenotype
 */ 

process simulate_test {

    tag { "$par|$C" }

    input:
    tuple id,n,q,PTgen,GTgen,file(bed),file(bim),file(fam),file(kinship),file(eigenval) from gt2pt_ch
    each c from Channel.fromList(1..params.c)
    each C from grid.C
    each d from grid.delta
    each hk from grid.hk
    each y_var from grid.y_var
    each y_cor from grid.y_cor
    each c_dist from grid.c_dist
    each p_loc from grid.p_loc
    each p_sd from grid.p_sd
    each p_dist from grid.p_dist      
    each lambda from grid.lambda
    each t from grid.t 
  
    output:
    tuple par_gemma, file('gemma.assoc.txt') optional true into gemma_v_ch
    tuple par_mlm, file('mlm.assoc.txt') into mlm_v_ch
    tuple par_manova, file('manova.assoc.txt') into manova_v_ch

    script:
    par = "$n|$q|$PTgen|$GTgen|$d|$hk|$y_var|$y_cor|$c_dist|$p_loc|$p_sd|$p_dist|$lambda|$t"
    par_gemma = "$par|GEMMA"
    par_mlm = "$par|MLM_$C"
    par_manova = "$par|MANOVA_$C" 
    pids = (1..q.toInteger()).join(' ')
    single = C
    if(grid.C instanceof List) {
       single = grid.C[0]
    }
    if(params.scale == true){scale = "--scale"} else {scale = ""}
    """ 
    # Manage chunks
    start=\$(( ($c-1)*(${params.p}/${params.c}) + 1 ))
    end=\$(( $c*(${params.p}/${params.c}) ))
 
    for (( v=\$start; v<=(\$end); v++ )); do
       # Extract single variant
       sed -n \${v}p $bim | awk '{print \$1"\t"\$4-1"\t"\$4}' > variant.bed
       plink2 --bfile \$(basename $bed | sed 's/.bed//') --extract bed0 variant.bed --out geno --make-bed --threads 1

       # Double check that is unique
       if [[ \$(cat geno.bim | wc -l) -gt 1 ]]; then continue; fi
 
       # Simulate phenotype
       simulatePT.R -s \$v -n $n -q $q -m $PTgen --geno geno -o pheno.txt -d $d -H $hk -v $y_var -c $y_cor -D $c_dist -p $p_loc -S $p_sd --p_dist $p_dist -l $lambda -t $t --fx ${params.fx}
    
       # Run MLM/MANOVA
       if [[ $C == 'PCA' ]]; then
          mlm.R -p pheno.txt -g geno -c $eigenval -k ${params.k} --mlm mlm_\$v.assoc.txt --manova manova_\$v.assoc.txt $scale 
       else 
          mlm.R -p pheno.txt -g geno --mlm mlm_\$v.assoc.txt --manova manova_\$v.assoc.txt $scale 
       fi

       # Run GEMMA once
       if [[ $single == $C ]]; then
          paste <(cut -f1-5 geno.fam) pheno.txt > tmpfile; mv tmpfile geno.fam
          (timeout 120 gemma -lmm -b geno -k $kinship -n $pids -outdir . -o gemma_\$v &> STATUS || exit 0)
          if [[ \$(grep ERROR STATUS) ]]; then
             touch gemma_\$v.assoc.txt
             continue
          else
             gemma -lmm -b geno -k $kinship -n $pids -outdir . -o gemma_\$v
             sed '1d' gemma_\$v.assoc.txt | awk '{print \$2"\t"\$NF}' > tmpfile; mv tmpfile gemma_\$v.assoc.txt
          fi
       fi

    done
    cat mlm_*.assoc.txt > mlm.assoc.txt
    cat manova_*.assoc.txt > manova.assoc.txt
    if [[ $single == $C ]]; then
       cat gemma_*.assoc.txt > gemma.assoc.txt
    fi
    """
}

gemma_v_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{gemma_ch}
mlm_v_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{mlm_ch}
manova_v_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{manova_ch}

gemma_ch.concat(mlm_ch, manova_ch).set{tie_power_ch}

 
/*
 *  Compute TIE/power
 */

process tie {

    tag{ "$par|$m" }   
 
    input:
    each m from grid.m
    tuple par, file(assoc) from tie_power_ch
 
    output:
    file("res.txt") into out_ch  

    script:
    """
    tie.R -t $assoc -l ${params.level} -m $m -o tieORpower.txt
    paste <(echo -e "$par|$m" | sed 's/|/\t/g') tieORpower.txt > res.txt
    """
}

out_ch.collectFile(name: "${params.out}", sort: { it.text }).set{pub_ch}

process end {

   publishDir "${params.dir}", mode: 'copy'

   input:
   file(sim) from pub_ch

   output:
   file(sim) into end_ch

   script:
   """
   sed -i "1 s/^/n\tq\tPTgen\tGTgen\tdelta\thk\ty_var\ty_cor\tc_dist\tp_loc\tp_sd\tp_dist\tlambda\ttransf\tmethod\tmtc\ttie\\n/" $sim
   """
}
