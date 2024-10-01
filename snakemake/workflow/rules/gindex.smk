rule gindexIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_folder, "gindex"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_folderr, "index", "graph.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex", "index.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} -i {input.g} -o {params.p} -t {threads}
        """

rule gindexFastIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_fast_folder, "gindex"),
    params:
        p=os.path.join(gindex_fast_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_fast_folderr, "index", "graph.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex_fast", "index.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} -i {input.g} -o {params.p} -t {threads}
        """

rule gindexCacheIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_cache_folder, "gindex"),
    params:
        p=os.path.join(gindex_cache_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_cache_folderr, "index", "graph.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex_cache", "index.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} index -i {input.g} -o {params.p} -t {threads} -c {cache_l}
        """


rule gindexQuery:
    input:
        exe_i=os.path.join(gindex_folder, "gindexquery"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
    params:
        p=os.path.join(gindex_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex", "query_{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} -i {params.p} -q {input.f} > {output}
        """

rule gindexFastQuery:
    input:
        exe_i=os.path.join(gindex_fast_folder, "gindexquery"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
    params:
        p=os.path.join(gindex_fast_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_fast_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_fast", "query_{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} -i {params.p} -q {input.f} > {output}
        """

rule gindexCacheQuery:
    input:
        exe_i=os.path.join(gindex_cache_folder, "gindex"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
    params:
        p=os.path.join(gindex_cache_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_cache_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_cache", "query_{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} query -i {params.p} -q {input.f} > {output}
        """
