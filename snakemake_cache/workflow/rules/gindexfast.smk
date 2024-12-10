rule gindexFastIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_fast_folder, "gindex"),
    params:
        p=os.path.join(gindex_fast_folderr, "index", "graph"),
        t=threads
    output:
        b=os.path.join(gindex_fast_folderr, "index", "graph.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex_fast", "index", "index.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} -i {input.g} -o {params.p} -t {params.t}
        """

rule gindexFastQuery:
    input:
        exe_i=os.path.join(gindex_fast_folder, "gindexquery"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
        b=os.path.join(gindex_fast_folderr, "index", "graph.bwt"),
    params:
        p=os.path.join(gindex_fast_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_fast_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_fast", "query", "{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} -i {params.p} -q {input.f} > {output}
        """
