rule downloadGindex:
    output:
        d=directory(gindex_folder),
        exe_i=os.path.join(gindex_folder, "gindex"),
    conda: "../envs/compilation.yml"
    shell:
        """
        git clone --recursive https://github.com/dlcgold/graph_index.git {output.d}
        cd {output.d}
        git submodule update --init --recursive
        cd gfatools
        make
        cd ../mfmi
        make
        cd ..
        make
        """


rule generateReads:
    input:
        os.path.join(gfa),
    output:
        os.path.join(input_folder, "q_{q}", "reads.fa"),
    conda:
         "../envs/gen.yml"
    shell:
        """
        python workflow/scripts/read_gen.py {input} {n_queries} {wildcards.q} {wildcards.q} > {output}
        """
