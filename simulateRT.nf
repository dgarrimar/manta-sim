#!/bin/env nextflow

/*
 * Copyright (c) 2021, Diego Garrido-Martín
 *
 * Simulation setting to study the running time of the asymptotic
 * PERMANOVA test (as implemented in manta) with respect to 
 * permutations (as implemented in vegan::adonis)
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

// General params
params.dir = 'result'
params.out = 'runtime.tsv'
params.fx = "$baseDir/supp"
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
    log.info 'Asymptotic Distribution of PERMANOVA Statistic in Complex Models'
    log.info '======================================================================='
    log.info 'Simulate and measure the running time of asymptotic PERMANOVA'
    log.info 'with respect to permutations'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run simulateRT.nf [options]'
    log.info ''
    log.info 'General parameters:'
    log.info ' --n SAMPLE_SIZE             number of samples (default: 100)'
    log.info ' --q RESPONSES               number of response variables (default: 3)'
    log.info ' --P PERMUTATIONS            number of permutations (default: 1000)'
    log.info ' --r REPLICATES              number of replicates (default: 1)'
    log.info ' --fx FUNCTIONS              path to helper functions and precomputed datasets (default: ./supp)'
    log.info ' --dir DIRECTORY             output directory (default: result)'
    log.info ' --out OUTPUT                output file (default: runtime.tsv)'
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
log.info "Helper functions             : ${params.fx}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''


/*
 *  Expand parameters
 */

def grid = [:]
params.keySet().each {
    if (it in ['n', 'q', 'P']) {
        grid[it] = params[it]
    }
}

grid.keySet().each {
    if (grid[it] =~ /,/) {
        grid[it] = grid[it].tokenize(',')
    } else if (grid[it] =~ /:/) {
        def (start, end, step) = grid[it].tokenize(':')
        def seq = []
        def val = start.toFloat()
        while (val.round(4) <= end.toFloat().round(4)) {
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
    runtime.R -n $n -q $q -r $r -P $P --fx ${params.fx} -o rt.txt
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
    sed -i '1 s/^/r\tn\tq\tP\tmanta\tadonis\\n/' ${rt}
    """
}

