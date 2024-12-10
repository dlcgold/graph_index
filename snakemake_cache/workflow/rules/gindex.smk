rule gindexIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_folder, "gindex"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph"),
        t=threads
    output:
        b=os.path.join(gindex_folderr, "index", "graph.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex", "index", "index.time"),
    shell:
        """
          /usr/bin/time --verbose -o {log.time} {input.exe_i} index -i {input.g} -o {params.p} -t {params.t} -c {cache_l}
        """

rule gindexIndexNoCache:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_folder, "gindex"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph_nc"),
        t=threads
    output:
        b=os.path.join(gindex_folderr, "index", "graph_nc.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex_no_cache", "index", "index.time"),
    shell:
        """
          /usr/bin/time --verbose -o {log.time} {input.exe_i} index -i {input.g} -o {params.p} -t {params.t}
        """

rule gindexQueryCache:
    input:
        exe_i=os.path.join(gindex_folder, "gindex"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
        b=os.path.join(gindex_folderr, "index", "graph.bwt"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_cache_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_cache", "query", "{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} query -i {params.p} -q {input.f} > {output}
        """


rule gindexQueryFast:
    input:
        exe_i=os.path.join(gindex_folder, "gindex"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
        b=os.path.join(gindex_folderr, "index", "graph.bwt"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_fast_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_fast", "query", "{q}.time"),
    shell:
        """
          /usr/bin/time --verbose -o {log.time} {input.exe_i} query -i {params.p} -q {input.f} -n > {output}
        """

rule gindexQueryFull:
    input:
        exe_i=os.path.join(gindex_folder, "gindex"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
        b=os.path.join(gindex_folderr, "index", "graph_nc.bwt"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_full_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_full", "query", "{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} query -i {params.p} -q {input.f} -f > {output}
        """
