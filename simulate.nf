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
params.p = 100
params.PTgen = 'norm-0-1'    
params.GTgen = 'simPopStructure_chr22' 
params.level = 0.05
params.m = 'none'
params.hs2 = 0 
params.hg2 = 0
params.t = 'none'
params.k = 20
params.c = 10

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
  log.info ' --p VARIANTS                number of variants to test (default: 100)'
  log.info ' --PTgen PHENO GENERATION    phenotype data generation to build error copula  (default: norm-0-1)'
  log.info ' --GTgen GENO GENERATION     genotype data generation: simPopStructure, simUnrelated, simRelated, simEmpirical (default: simPopStructure)'
  log.info ' --GTdir GT DIR              input genotype directory with files generated by simulateGT.nf (default: data)' 
  log.info ' --level SIG LEVEL           significance level (default: 0.05)'
  log.info ' --m MULTIPLE TESTING        perform multiple testing correction: none or any method in R:p.adjust (default: none)'
  log.info ' --hs2 SNP HERITABILITY      average fraction of variance explained by causal variants across traits (default: 0.01)'
  log.info ' --hg2 REL HERITABILITY      average fraction of variance explained by relatedness across traits (default: 0.3)'
  log.info ' --t TRANSFORM               transformation of response variables: none, PCA, GAMMA (default: none)'
  log.info ' --k NUMBER PC               Number of PCs used when transformation = PCA (default: 20)'
  log.info ' --c NUMBER CHUNKS           Number of chunks (default: 10)'
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
log.info "Phenotype data generation    : ${params.PTgen}"
log.info "Genotype data generation     : ${params.GTgen}"
log.info "Genotype input directory     : ${params.GTdir}"
log.info "Significance level           : ${params.level}"
log.info "Multiple testing correction  : ${params.m}"
log.info "Causal variant heritability  : ${params.hs2}"
log.info "Relatedness heritability     : ${params.hg2}"
log.info "Transformation               : ${params.t}"
log.info "Number of PCs                : ${params.k}"
log.info "Number of chunks             : ${params.c}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
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
  if(it in ['n','q','PTgen','GTgen','hs2','hg2','t','m']){
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
.combine(grid2ch.hs2)
.combine(grid2ch.hg2)
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

    plink2 --bfile $prefix --indep-pairwise 50 5 0.1 --keep keep.txt --out $GTgen
    plink2 --bfile $prefix --extract ${GTgen}.prune.in --keep keep.txt --out ${GTgen}.pruned --make-bed
    
    # Thin: subset n and p 
    plink2 --bfile $prefix --seed 1 --thin-count $params.p --keep keep.txt --out ${GTgen}.thin --make-bed    
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
    tuple val(id), file("${GTgen}.sXX.txt") into kinship_ch

    script:
    GTgen = id.split("\\|")[0]
    """
    # Compute kinship
    sed -i 's/-9/1/' $fam    
    gemma -gk 2 -bfile \$(basename $bed | sed 's/.bed//') -outdir . -o $GTgen 
    """
}

if("PCA" in grid.t) {

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
        plink2 --bfile ${GTgen}.pruned --pca ${params.k} \$approx --out ${GTgen}
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

    tag { par }

    input:
    tuple id,n,q,PTgen,GTgen,hs2,hg2,file(bed),file(bim),file(fam),file(kinship),file(eigenval) from gt2pt_ch
    each c from Channel.fromList(1..params.c)
    each t from grid.t
   
    output:
    tuple par_gemma, file('gemma.assoc.txt') optional true into gemma_v_ch
    tuple par_mlm, file('mlm.assoc.txt') into mlm_v_ch
    tuple par_manova, file('manova.assoc.txt') into manova_v_ch

    script:
    par = "$n|$q|$PTgen|$GTgen|$hs2|$hg2"
    par_gemma = "$par|GEMMA"
    par_mlm = "$par|MLM_$t"
    par_manova = "$par|MANOVA_$t" 
    pids = (1..q.toInteger()).join(' ')
    single = t
    if(grid.t instanceof List) {
       single = grid.t[0]
    }
    """ 
    # Manage chunks
    start=\$(( ($c-1)*(${params.p}/${params.c}) + 1 ))
    end=\$(( $c*(${params.p}/${params.c}) ))
 
    for (( v=\$start; v<=(\$end); v++ )); do
       # Extract single variant
       sed -n \${v}p $bim | cut -f2 > variant.txt    
       plink2 --bfile \$(basename $bed | sed 's/.bed//') --extract variant.txt --out geno --make-bed   
   
       # Simulate phenotype
       simulatePT.R -s \$v -n $n -q $q --PTgen $PTgen --geno geno --kinship $kinship --hs2 $hs2 --hg2 $hg2 -o pheno.txt 
       paste <(cut -f1-5 geno.fam) pheno.txt > tmpfile; mv tmpfile geno.fam
    
       # Run GEMMA once
       if [[ $single == $t ]]; then
          gemma -lmm -b geno -k $kinship -n $pids -outdir . -o gemma_\$v
          sed '1d' gemma_\$v.assoc.txt | awk '{print \$2"\t"\$NF}' > tmpfile; mv tmpfile gemma_\$v.assoc.txt
       fi   
 
       # Run MLM/MANOVA with transformation
       if [[ $t == "GAMMA" ]]; then
          paste <(cut -f1-5 geno.fam) pheno.txt > tmpfile; mv tmpfile geno.fam
          for i in {1..$q}; do
             gemma -vc 2 -p pheno.txt -k $kinship -n \$i -outdir . -o VC &> /dev/null
             grep -F "sigma2 estimates =" VC.log.txt | cut -d ' ' -f 7,9
          done > VC.txt
       else
          touch VC.txt
       fi
       mlm.R -p pheno.txt -g geno -t $t -c $eigenval -n ${params.k} -k $kinship -v VC.txt --mlm mlm_\$v.assoc.txt --manova manova_\$v.assoc.txt --scale
    done
    for method in {gemma,mlm,manova}; do
       cat \${method}_*.assoc.txt > \${method}.assoc.txt
    done
    """
}

gemma_v_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.view().set{gemma_ch}
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
   sed -i "1 s/^/n\tq\tPTgen\tGTgen\ths2\thg2\tmethod\tmtc\ttie\\n/" $sim
   """
}


