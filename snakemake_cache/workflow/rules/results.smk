rule extractVerboseIndex:
    input:
        time = os.path.join(output_folder, "bench", "{tool}", "index", "index_{c}.time"),
    output:
        os.path.join(output_folder, "bench", "{tool}", "index", "index_{c}.csv")
    shell:
        """
        python workflow/scripts/time_verbose_extractor.py {wildcards.tool}_{wildcards.c}  index {n_queries} 0 < {input.time} > {output}
        """
rule extractVerboseQuery:
    input:
        time = os.path.join(output_folder, "bench", "{tool}", "query_{c}", "{q}.time"),
    output:
        os.path.join(output_folder, "bench", "{tool}", "query_{c}", "{q}.csv")
    shell:
        """
        python workflow/scripts/time_verbose_extractor.py {wildcards.tool}_{wildcards.c} query {n_queries} {wildcards.q} < {input.time} > {output}
        """

rule mergeIndexTime:
    input:
        expand(
            os.path.join(output_folder, "bench", "{tool}", "index", "index_{c}.csv"),
            tool = ["gindex"],
            c = CACHES
        ),
    output:
        os.path.join(results_folder, "index.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule mergeQueryTime:
    input:
        expand(
            os.path.join(output_folder, "bench", "{tool}", "query_{c}", "{q}.csv"),
            tool = ["gindex"],
            q=QUERIES,
            c=CACHES
        ),
    output:
        os.path.join(results_folder, "query.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """
rule plot:
    input:
        i = os.path.join(results_folder, "index.csv"),
        q = os.path.join(results_folder, "query.csv"),
    params:
        r = os.path.join(results_folder),
        i = os.path.join(gindex_folderr, "index"),
    output:
        os.path.join(results_folder, "query_mem.pdf"),
        os.path.join(results_folder, "query_time.pdf"),
        os.path.join(results_folder, "index.pdf"),
    conda: "../envs/results.yml"
    shell:
        """
        python workflow/scripts/plot.py {input.i} {input.q} {params.r} {n_queries} {params.i}
        """
