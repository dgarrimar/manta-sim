/*
 * Asymptotic distribution of Anderson's test statistic in complex models
 * Simulation setting for H0
 * v1.0
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

params.q = "3,5"
params.r = "10,20"
params.s = 100000
params.model = "Y~A+B+AB.H0.R"
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
         s: ${params.s}
         model: ${params.model}
         transf: ${params.transf}
         gen: ${params.gen}
         pub: ${params.pub}
         out: ${params.out}
         """
         .stripIndent()

/*
 *  Create channels given parameters
 */

q_ch = Channel.from(params.q.tokenize(","))
r_ch = Channel.from(params.r.tokenize(","))

q_ch.combine(r_ch).set{grid_ch}


/*
 *  Simulate
 */

process simulation {
 
    tag {"q=$q, r=$r"}
  
    input:
    set val(q), val(r) from grid_ch
    
    output:    
    set val(q), val(r), file("${q}.${r}.out") into res_ch
    
    script:
    """
    ${params.model} -q $q -r $r -s ${params.s} -m ${params.gen} -t ${params.transf} -o "${q}.${r}.out"
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

