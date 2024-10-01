rule GFAtoVG:
    input:
        os.path.join(gfa),
    output:
        vg=os.path.join(graphpp_folderr, "index", "graph.vg")
    conda:
        "../envs/vg.yml"
    shell:
        """
        vg convert -g {input} -v > {output.vg}
        """
rule graphppIndex:
    input:
        os.path.join(graphpp_folderr, "index", "graph.vg")
    output:
        os.path.join(graphpp_folderr, "index", "graph.gcsa"),
    conda:
        "../envs/vg.yml"
    log:
        time=os.path.join(bench_folder, "graphpp", "index.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} vg index -g {output} -t {threads} {input}
        """

rule graphppQuery:
    input:
        g=os.path.join(graphpp_folderr, "index", "graph.gcsa"),
        exe_i=os.path.join(graphpp_folder, "graphpp"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
    output:
        os.path.join(graphpp_folderr, "query", "query_{q}.m"),
    conda:
        "../envs/vg.yml"
    log:
        time=os.path.join(bench_folder, "graphpp", "query_{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} {input.g} {input.f} > {output}
        """
