rule gindexCacheIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_cache_folder, "gindex"),
    params:
        p=os.path.join(gindex_cache_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_cache_folderr, "index", "graph.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex_cache", "index", "index.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} index -i {input.g} -o {params.p} -t {threads} -c {cache_l}
        """

rule gindexCacheQuery:
    input:
        exe_i=os.path.join(gindex_cache_folder, "gindex"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
        b=os.path.join(gindex_cache_folderr, "index", "graph.bwt"),
    params:
        p=os.path.join(gindex_cache_folderr, "index", "graph"),
    output:
        b=os.path.join(gindex_cache_folderr, "query", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex_cache", "query", "{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} query -i {params.p} -q {input.f} > {output}
        """
