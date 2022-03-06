#!/bin/env nextflow

/*
 * Copyright (c) 2021, Diego Garrido-Martín
 * 
 * Simulation setting I to study the asymptotic distribution of the
 * PERMANOVA test statistic in complex models
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
params.out = 'simulation.tsv'
params.fx = "$baseDir/supp"
params.help = false

// Simulation params
params.a = 2
params.b = 3
params.n = 100
params.u = 1
params.q = 3
params.t = 'none'
params.which = 'B'
params.model = 'Y~A+B+AB.R'
params.gen = 'mvnorm'
params.sim = 10000
params.delta = 0
params.hk = 1
params.adonis = 0
params.cqf = 'farebrother-SS'
params.ir = 0

// Generation: multivariate normal
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

// Model: Y~C+A.R
params.r = 0
params.C_mean = 0
params.C_var = 1


/*
 *  Print usage and help
 */

if (params.help) {
    log.info ''
    log.info 'Asymptotic Distribution of PERMANOVA Test Statistic in Complex Models'
    log.info '======================================================================='
    log.info 'Simulations given a set of parameters (I)'
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
    log.info ' --which WHICH               which factor changes in H1: A, B or AB (default: B)'
    log.info ' --model MODEL               R script with model definition (default: Y~A+B+AB.R)'
    log.info ' --gen GENERATION            data generation: mvnorm, simplex, multinom or copula (default: mvnorm)'
    log.info ' --sim SIMULATIONS           number of simulations (default: 10000)'
    log.info ' --delta DELTA               change in H1 (default: 0)'
    log.info ' --hk HETEROSKEDASTICITY     heteroskedasticity, 1 is homoskedastic (default: 1)'
    log.info ' --adonis PERMUTATIONS       should permutation test be performed? Specify number of permutations (default: 0)'
    log.info ' --cqf CompQuadForm          method employed to compute the quadratic form (default: farebrother-SS)'
    log.info ' --ir IRREPRODUCIBLE         should random seed be set to Sys.time() (default: 0)'
    log.info ' --fx FUNCTIONS              path to helper functions and precomputed datasets (default: ./supp)'
    log.info ' --dir DIRECTORY             output directory (default: result)'
    log.info ' --out OUTPUT                output file (default: simulation.tsv)'
    log.info ''
    log.info 'Additional parameters for generation = mvnorm:'
    log.info ' --y_var VARIANCE            variance of response variables: equal or unequal (default: equal)'
    log.info ' --y_cor CORRELATION         correlation of response variables (default: 0)'
    log.info ''
    log.info 'Additional parameters for generation = copula:'
    log.info ' --y_var VARIANCE            variance of response variables: equal or unequal (default: equal)'
    log.info ' --y_cor CORRELATION         correlation of response variables (default: 0)'
    log.info ' --c_dist DISTRIBUTION       multivariate non-normal distribution definition (default: unif-0-1)'
    log.info ''
    log.info 'Additional parameters for generation = simplex:'
    log.info ' --p_loc LOCATION            location in the simplex of the generator model, 1 is centered (default: 1)'
    log.info ' --p_sd STDEV                standard deviation of the generator model (default: 0.1)'
    log.info ' --p_dist DISTRIBUTION       distribution of the generator model (default: norm)'
    log.info ''
    log.info 'Additional parameters for generation = multinom:'
    log.info ' --lambda LAMBDA             multinomial distribution\'s size parameter (default: 1000)'
    log.info ''
    log.info 'Additional parameters for model = Y~C+A.R'
    log.info ' --r CORRELATION             (if generation = simplex: amount of noise added to the) correlation between C and Y1 (default: 0)'
    log.info ' --C_mean MEAN               mean of covariate C (default: 0)'
    log.info ' --C_var VARIANCE            variance of covariate C (default: 1)'
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
log.info "Which factor changes in H1   : ${params.which}"
log.info "Model definition             : ${params.model}"
log.info "Data generation              : ${params.gen}"
log.info "Number of simulations        : ${params.sim}"
log.info "Change in H1 (delta)         : ${params.delta}"
log.info "Heteroskedasticity           : ${params.hk}"
log.info "Adonis permutation test      : ${params.adonis}"
log.info "CompQuadForm                 : ${params.cqf}"
log.info "Irreproducible               : ${params.ir}"
log.info "Helper functions             : ${params.fx}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''
if (params.gen == 'mvnorm') {
    log.info 'Additional parameters'
    log.info '---------------------'
    log.info "Variance of Y variables      : ${params.y_var}"
    log.info "Correlation of Y variables   : ${params.y_cor}"
} else if (params.gen == 'copula') {
    log.info 'Additional parameters'
    log.info '---------------------'
    log.info "Variance of Y variables      : ${params.y_var}"
    log.info "Correlation of Y variables   : ${params.y_cor}"
    log.info "Distribution definition      : ${params.c_dist}"
} else if (params.gen == 'simplex') {
    log.info 'Additional parameters'
    log.info '---------------------'
    log.info "Location in the simplex      : ${params.p_loc}"
    log.info "Generation stdev             : ${params.p_sd}"
    log.info "Generation distribution      : ${params.p_dist}"
} else if (params.gen == 'multinom') {
    log.info 'Additional parameters'
    log.info '---------------------'
    log.info "Lambda                       : ${params.lambda}"
}
if(params.model == 'Y~C+A.R') {
    if(params.gen == 'simplex') {
        log.info "(Noise added to) cor(C, Y1)  : ${params.r}"
    }
    log.info "cor(C, Y[,1])                : ${params.r}"
    log.info "Mean of covariate C          : ${params.C_mean}"
    log.info "Variance of covariate C      : ${params.C_var}"
}
log.info ''


/*
 *  Expand parameters
 */

def grid = [:]
params.keySet().each {
    if (it in ['a','b','n','u','q','delta','hk','y_var','y_cor','p_loc','p_sd','p_dist','c_dist','lambda','r','C_mean','C_var','t']) {
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
    each a from grid.a
    each b from grid.b
    each n from grid.n
    each u from grid.u
    each q from grid.q
    each d from grid.delta
    each hk from grid.hk
    each y_var from grid.y_var
    each y_cor from grid.y_cor
    each p_loc from grid.p_loc
    each p_sd from grid.p_sd
    each p_dist from grid.p_dist      
    each lambda from grid.lambda
    each c_dist from grid.c_dist
    each r from grid.r
    each C_mean from grid.C_mean
    each C_var from grid.C_var
    each t from grid.t

    output:
    
    file('sim.txt') into sim_ch

    script:
    """
    ${params.model} -a $a -b $b -n $n -u $u -q $q -d $d -H $hk -v $y_var -c $y_cor -p $p_loc -s $p_sd --p_dist $p_dist -l $lambda -D $c_dist -r $r --C_mean $C_mean --C_var $C_var --adonis ${params.adonis} -S ${params.sim} -i ${params.ir} -m ${params.gen} -t $t -w ${params.which} -Q ${params.cqf} -o sim.txt --fx ${params.fx}
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

    if(params.gen == "mvnorm" & params.model == "Y~A+B+AB.R")
        """
        sed -i "1 s/^/a\tb\tn\tu\tq\tdelta\thk\ty_var\ty_cor\tt\tA\tB\tAB\tA_manova\tB_manova\tAB_manova\\n/" ${simulation}
        """
    else if (params.gen == "copula" & params.model == "Y~A+B+AB.R")
        """
        sed -i "1 s/^/a\tb\tn\tu\tq\tdelta\thk\ty_var\ty_cor\tc_dist\tt\tA\tB\tAB\tA_manova\tB_manova\tAB_manova\\n/" ${simulation}
        """
    else if (params.gen == "simplex" & params.model == "Y~A+B+AB.R")
        """
        sed -i "1 s/^/a\tb\tn\tu\tq\tdelta\thk\tp_loc\tp_dist\tp_sd\tt\tA\tB\tAB\tA_manova\tB_manova\tAB_manova\\n/" ${simulation}
        """
    else if (params.gen == "multinom" & params.model == "Y~A+B+AB.R")
        """
        sed -i "1 s/^/a\tb\tn\tu\tq\tdelta\thk\tp_loc\tlambda\tt\tA\tB\tAB\tA_manova\tB_manova\tAB_manova\\n/" ${simulation}
        """
    else if(params.gen == "mvnorm" & params.model == "Y~C+A.R")
        """
        sed -i "1 s/^/a\tC_mean\tC_var\tn\tu\tq\tdelta\tr\thk\ty_var\ty_cor\tt\tC\tA\tC_manova\tA_manova\\n/" ${simulation}
        """
    else if (params.gen == "copula" & params.model == "Y~C+A.R")
        """
        sed -i "1 s/^/a\tC_mean\tC_var\tn\tu\tq\tdelta\tr\thk\ty_var\ty_cor\tc_dist\tt\tC\tA\tC_manova\tA_manova\\n/" ${simulation}
        """
    else if (params.gen == "simplex" & params.model == "Y~C+A.R")
        """
        sed -i "1 s/^/a\tC_mean\tC_var\tn\tu\tq\tdelta\tr\thk\tp_loc\tp_dist\tp_sd\tt\tC\tA\tC_manova\tA_manova\\n/" ${simulation}
        """
    else if (params.gen == "multinom" & params.model == "Y~C+A.R")
        """
        sed -i "1 s/^/a\tC_mean\tC_var\tn\tu\tq\tdelta\tr\thk\tp_loc\tlambda\tt\tC\tA\tC_manova\tA_manova\\n/" ${simulation}
        """
}
