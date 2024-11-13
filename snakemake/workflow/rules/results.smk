rule extractVerboseIndex:
    input:
        time = os.path.join(output_folder, "bench", "{tool}", "index", "index.time"),
    output:
        os.path.join(output_folder, "bench", "{tool}", "index", "index.csv")
    shell:
        """
        python workflow/scripts/time_verbose_extractor.py {wildcards.tool} index {n_queries} 0 < {input.time} > {output}
        """
rule extractVerboseQuery:
    input:
        time = os.path.join(output_folder, "bench", "{tool}", "query", "{q}.time"),
    output:
        os.path.join(output_folder, "bench", "{tool}", "query", "{q}.csv")
    shell:
        """
        python workflow/scripts/time_verbose_extractor.py {wildcards.tool} query {n_queries} {wildcards.q} < {input.time} > {output}
        """

rule mergeIndexTime:
    input:
        expand(
            os.path.join(output_folder, "bench", "{tool}", "index", "index.csv"),
            tool = ["gindex", "gindex_fast", "gindex_cache", "gindex_merge", "graphpp"],
        ),
    output:
        os.path.join(results_folder, "results", "index.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule mergeQueryTime:
    input:
        expand(
            os.path.join(output_folder, "bench", "{tool}", "query", "{q}.csv"),
            tool = ["gindex", "gindex_fast", "gindex_cache", "gindex_merge", "graphpp"],
            q=QUERIES
        ),
    output:
        os.path.join(results_folder, "results", "query.csv"),
    conda: "../envs/csvkit.yml"
    shell:
        """
        csvstack {input} > {output}
        """

rule plot:
    input:
        i = os.path.join(results_folder, "results", "index.csv"),
        q = os.path.join(results_folder, "results", "query.csv"),
    params:
        os.path.join(results_folder, "results")
    output:
        os.path.join(results_folder, "results", "query_mem.pdf"),
        os.path.join(results_folder, "results", "query_time.pdf"),
        os.path.join(results_folder, "results", "index.pdf"),
    conda: "../envs/results.yml"
    shell:
        """
        python workflow/scripts/plot.py {input.i} {input.q} {params}
        """
