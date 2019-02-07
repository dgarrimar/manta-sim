/*
 * Asymptotic distribution of Anderson's test statistic in complex models
 * Simulation setting for H1
 * v1.0
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

params.q = "3,5"
params.r = "10,20"
params.d = "0.01,0.1"
params.s = 100000
params.model = "Y~A+B+AB.H1.R"
params.transf = "None"
params.gen = "mvnorm"
params.pub = "result"
params.out = "simulation.tsv"

println """\
         =============================
         S I M U L A T E - N F   
         =============================
         q: ${params.q}
         r: ${params.r}
         d: ${params.d}
         s: ${params.s}
         model: ${params.model}
         transf: ${params.transf}
         gen: ${params.gen}
         pub: ${params.pub}
         out: ${params.out}
         """
         .stripIndent()

/*
 *  Expand parameters
 */

par_in = [params.q, params.r, params.d]
par_l = par_in.size()
par_out = [""]* par_l

i = 0
par_in.each {
    if (it =~ /,/){
        par_out[i] = it.tokenize(',')
    } else if (it =~ /:/) {
        def (start, end, step) = it.tokenize(':')
        if(step?.trim()){ // not empty
            par_out[i] = (start.toInteger()..end.toInteger()).step(step.toInteger()).toList()
        }else{
            par_out[i] = (start.toInteger()..end.toInteger()).toList()
        }
    } else {
	par_out[i] = it
    }
i = i + 1
}

/*
 *  Create channels given parameters
 */


q_ch = Channel.from(par_out[0])
r_ch = Channel.from(par_out[1])
d_ch = Channel.from(par_out[2])

q_ch.combine(r_ch).combine(d_ch).set{grid_ch}


/*
 *  Simulate
 */

process simulation {
 
    tag {"q=$q, r=$r, d=$d"}
  
    input:
    set val(q), val(r), val(d) from grid_ch
    
    output:    
    set val(q), val(r), val(d), file("${q}.${r}.${d}.out") into res_ch
    
    script:
    """
    ${params.model} -q $q -r $r -d $d -s ${params.s} -m ${params.gen} -t ${params.transf} -o "${q}.${r}.${d}.out"
    """
}

res_ch.collectFile(name: "${params.out}"){(it[0..it.size-2] + it[it.size-1].text).join("\t")}.set{pub_ch}

process end {

   tag { "" }
   publishDir "${params.pub}"     

   input:
   file(simulation) from pub_ch

   output:
   file(simulation) into end_ch

   script:
   """
   echo "DONE"
   """
}

