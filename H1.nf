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
 *  Create channels given parameters
 */

q_ch = Channel.from(params.q.tokenize(","))
r_ch = Channel.from(params.r.tokenize(","))
d_ch = Channel.from(params.d.tokenize(","))

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

