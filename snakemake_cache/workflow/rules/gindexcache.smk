rule gindexCacheIndex:
    input:
        g=os.path.join(gfa),
        exe_i=os.path.join(gindex_folder, "gindex"),
    params:
        p=os.path.join(gindex_folderr, "index", "{c}"),
        t=threads
    output:
        b=os.path.join(gindex_folderr, "index", "{c}.bwt"),
    log:
        time=os.path.join(bench_folder, "gindex", "index", "index_{c}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} index -i {input.g} -o {params.p} -t {params.t} -c {wildcards.c}
        """

rule gindexCacheQuery:
    input:
        exe_i=os.path.join(gindex_folder, "gindex"),
        f=os.path.join(input_folder, "q_{q}", "reads.fa"),
        b=os.path.join(gindex_folderr, "index", "{c}.bwt"),
    params:
        p=os.path.join(gindex_folderr, "index", "{c}"),
    output:
        b=os.path.join(gindex_folderr, "query_{c}", "query_{q}.m"),
    log:
        time=os.path.join(bench_folder, "gindex", "query_{c}", "{q}.time"),
    shell:
        """
         /usr/bin/time --verbose -o {log.time} {input.exe_i} query -i {params.p} -q {input.f} > {output}
        """
