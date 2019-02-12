/*
 * Asymptotic distribution of Anderson's test statistic in complex models
 * Benchmarking setting 
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

params.q = 3
params.r = 50
params.d = 3
params.s = 1000000
params.k = 100
params.model = "Y~A+B+AB.bm.R"
params.transf = "None"
params.gen = "mvnorm"
params.pub = "result"

println """\
         =============================
         B E N C H M A R K - N F   
         =============================
         q: ${params.q}
         r: ${params.r}
         d: ${params.d}
         s: ${params.s}
         k: ${params.k}
         model: ${params.model}
         transf: ${params.transf}
         gen: ${params.gen}
         pub: ${params.pub}
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
	par_out[i] = it.toInteger()
    }
i = i + 1
}

/*
 *  Create channels given parameters
 */


q_ch = Channel.from(par_out[0])
r_ch = Channel.from(par_out[1])
d_ch = Channel.from(par_out[2])

q_ch.combine(r_ch).combine(d_ch).into{grid_ch;grid_ch2}

/*
 *  Simulate
 */

process getEmpirical {
 
    tag {"q=$q, r=$r, d=$d, chunk ${chunk}"}
  
    input:
    set val(q), val(r), val(d) from grid_ch
    each chunk from Channel.from(1..params.k)

    output:    
    set val(q), val(r), val(d), file("${q}.${r}.${d}.emp_$chunk") into emp_ch
    
    script:
    """
    s=\$(echo ${params.s}/${params.k} | bc -l)
    ${params.model} -q $q -r $r -d $d -s \$s -m ${params.gen} -t ${params.transf} -k $chunk -o "${q}.${r}.${d}.emp_$chunk"
    """
}

process getProposed {

    clusterOptions="-pe smp 10"
 
    tag {"q=$q, r=$r, d=$d"}

    input:
    set val(q), val(r), val(d) from grid_ch2
 
    output:
    set val(q), val(r), val(d), file("${q}_${r}_${d}.pro.RData") into pro_ch

    script:
    """
    ${params.model} -q $q -r $r -d $d -s ${params.s} -m ${params.gen} -t ${params.transf} -k 0 -c 10 -o "${q}_${r}_${d}.pro"
    """
}

emp_ch
.map{ [it[0..it.size-2].join("_"), it[it.size-1] ] }
.collectFile(sort: {it.name})
.map{[it.name, it]}
.join(pro_ch.map{ [it[0..it.size-2].join("_"), it[it.size-1] ] })
.set{pub_ch}

process end {

   tag { "q=" + par.tokenize("_")[0] + "," +
         "r=" + par.tokenize("_")[1] + "," +
         "d=" + par.tokenize("_")[2]
       }

   publishDir "${params.pub}/${par}"     

   input:
   set val(par), file(sim), file(sim2) from pub_ch

   output:
   set file(sim), file(sim2) into end_ch
   
   script:
   """
   echo "DONE"
   """
}


