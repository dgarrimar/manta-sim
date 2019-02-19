/*
 * Asymptotic distribution of Anderson's test statistic in complex models
 * Simulation setting
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.dir = 'result'
params.out = 'simulation.tsv'
params.help = false

// Simulation params
params.a = 2
params.b = 3
params.n = 100
params.u = 1
params.q = 3
params.t = 'none'
params.model = 'Y~A+B+AB.R'
params.gen = 'mvnorm'
params.sim = 10000
params.d = 0
params.hk = 1

// Generation: multivariate normal
params.y_var = 'equal'
params.y_cor = 0

// Generation: simplex (proportions)
params.p_loc = 1
params.p_sd = 0.1


/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'Asymptotic Distribution of Anderson\'s Test Statistic in Complex Models'
  log.info '======================================================================='
  log.info 'Simulations given a set of parameters'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run simulate.nf [options]'
  log.info ''
  log.info 'General parameters:'
  log.info ' --a LEVELS_A                levels of factor A (default: 2)'
  log.info ' --b LEVELS_B                levels of factor B (default: 3)'
  log.info ' --n SAMPLE_SIZE             total number of samples (default: 100)'
  log.info ' --u UNBALANCE               unbalance, 1 is balanced. (default: 1)'
  log.info ' --q RESPONSES               number of response variables (default: 3)'
  log.info ' --t TRANSFORM               transform response variables: none, sqrt or log (default: none)'
  log.info ' --model MODEL               R script with model definition (default: Y~A+B+AB.R)'
  log.info ' --gen GENERATION            data generation: mvnorm or simplex (default: mvnorm)'
  log.info ' --sim SIMULATIONS           number of simulations (default: 10000)'
  log.info ' --delta DELTA               change in H1 (default: 0.1)'
  log.info ' --hk HETEROSKEDASTICITY     heteroskedasticity, 1 is homoskedastic (default: 1)'
  log.info ' --dir DIRECTORY             output directory (default: result)'
  log.info ' --out OUTPUT                output file (default: simulation.tsv)'
  log.info ''
  log.info 'Additional parameters for generation = mvnorm:'
  log.info ' --y_var VARIANCE            variance of response variables: equal or unequal (default: equal)'
  log.info ' --y_cor CORRELATION         correlation of response variables (default: 0)'
  log.info ''
  log.info 'Additional parameters for generation = simplex:'
  log.info ' --p_loc LOCATION          location in the simplex of the generator model, 1 is centered (default: 1)'
  log.info ' --p_sd STDEV              standard deviation of the generator model (default: 0.1)'
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
log.info "Levels of factor A           : ${params.a}"
log.info "Levels of factor B           : ${params.b}"
log.info "Total sample size            : ${params.n}"
log.info "Unbalance degree             : ${params.u}"
log.info "Number of Y variables        : ${params.q}"
log.info "Transformation of Y          : ${params.t}"
log.info "Model definition             : ${params.model}"
log.info "Data generation              : ${params.gen}"
log.info "Number of simulations        : ${params.sim}"
log.info "Change in H1 (delta)         : ${params.d}"
log.info "Heteroskedasticity           : ${params.hk}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''
if(params.gen == 'mvnorm'){
  log.info 'Additional parameters'
  log.info '---------------------'
  log.info "Variance of Y variables      : ${params.y_var}"
  log.info "Correlation of Y variables   : ${params.y_cor}"
  log.info ''
} else if (params.gen == 'simplex'){
  log.info 'Additional parameters'
  log.info '---------------------'
  log.info "Location in the simplex      : ${params.p_loc}"
  log.info "Generation stdev             : ${params.p_sd}"
  log.info ''
}

/*
 *  Expand parameters
 */


def grid = [:]
params.keySet().each{
  if(it in ['a','b','n','u','q','d','hk','y_var','y_cor','p_loc','p_sd']){
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
 *  Simulate
 */

process simulation {
 
    input:
    each a from grid.a
    each b from grid.b
    each n from grid.n
    each u from grid.u
    each q from grid.q
    each d from grid.d
    each hk from grid.hk
    each y_var from grid.y_var
    each y_cor from grid.y_cor
    each p_loc from grid.p_loc
    each p_sd from grid.p_sd      

    output:
    
    file('sim.txt') into sim_ch

    script:
    """
    ${params.model} -a $a -b $b -n $n -u $u -q $q -d $d -H $hk -v $y_var -c $y_cor -l $p_loc -s $p_sd -S ${params.sim} -m ${params.gen} -t ${params.t} -o sim.txt
    """
}

sim_ch.collectFile(name: "${params.out}").set{pub_ch}

process end {

   publishDir "${params.dir}"     

   input:
   file(simulation) from pub_ch

   output:
   file(simulation) into end_ch

   script:

   if(params.gen == "mvnorm")
     """
     sed -i "1 s/^/a\tb\tn\tu\tq\tdelta\thk\ty_var\ty_cor\tt\tA\tB\tAB\tA_manova\tB_manova\tAB_manova\\n/" ${simulation}
     """
   else if (params.gen == "simplex")
     """
     sed -i "1 s/^/a\tb\tn\tu\tq\tdelta\thk\tp_loc\tp_sd\tt\tA\tB\tAB\tA_manova\tB_manova\tAB_manova\\n/" ${simulation}
     """
}

