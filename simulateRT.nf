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
params.out = 'runtime.tsv'
params.help = false

// Simulation params
params.r = 1
params.n = 100
params.q = 3
params.P = 1000

/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'Asymptotic Distribution of Anderson\'s Test Statistic in Complex Models'
  log.info '======================================================================='
  log.info 'Simulate and measure running time wrt the permutation'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run simulateRT.nf [options]'
  log.info ''
  log.info 'General parameters:'
  log.info ' --n SAMPLE_SIZE             number of samples (default: 100)'
  log.info ' --q RESPONSES               number of response variables (default: 3)'
  log.info ' --P PERMUTATIONS            number of permutations (default: 1000)'
  log.info ' --r REPLICATES              number of replicates (default: 1)'
  log.info ' --dir DIRECTORY             output directory (default: result)'
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
log.info "Number of Y variables        : ${params.q}"
log.info "Number of permutations       : ${params.P}"
log.info "Number of replicates         : ${params.r}"
log.info ''

/*
 *  Expand parameters
 */


def grid = [:]
params.keySet().each{
  if(it in ['n', 'q', 'P']){
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
    each n from grid.n
    each q from grid.q
    each r from params.r
    each P from grid.P

    output:
    
    file('rt.txt') into rt_ch

    script:
    """
    runtime.R -n $n -q $q -r $r -P $P -o rt.txt
    """
}

rt_ch.collectFile(name: "${params.out}").set{pub_ch}

process end {

   publishDir "${params.dir}"     

   input:
   file(rt) from pub_ch

   output:
   file(rt) into end_ch

   script:

   """
   sed -i '1 s/^/r\tn\tq\tP\tmlm\tadonis\\n/' ${rt}
   """
}

