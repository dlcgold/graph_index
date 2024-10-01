rule downloadGindex:
    output:
        d=directory(gindex_folder),
        exe_i=os.path.join(gindex_folder, "gindex"),
        exe_q=os.path.join(gindex_folder, "gindexquery")
    conda: "../envs/compilation.yml"
    shell:
        """
        git clone --recursive https://github.com/dlcgold/graph_index.git {output.d}
        cd {output.d}
        cd gfatools
        make
        cd ../mfmi
        make
        cd ..
        make
        """

rule downloadGindex_fast:
    output:
        d=directory(gindex_fast_folder),
        exe_i=os.path.join(gindex_fast_folder, "gindex"),
        exe_q=os.path.join(gindex_fast_folder, "gindexquery")
    conda: "../envs/compilation.yml"
    shell:
        """
        git clone --recursive https://github.com/dlcgold/graph_index.git {output.d}
        cd {output.d}
        git checkout fast
        cd gfatools
        make
        cd ../mfmi
        make
        cd ..
        make
        """

rule downloadGindex_cache:
    output:
        d=directory(gindex_cache_folder),
        exe_i=os.path.join(gindex_cache_folder, "gindex")
    conda: "../envs/compilation.yml"
    shell:
        """
        git clone --recursive https://github.com/dlcgold/graph_index.git {output.d}
        cd {output.d}
        git checkout cache
        git submodule update --init --recursive
        cd gfatools
        make
        cd ../mfmi
        make
        cd ..
        make
        """


rule downloadGraphpp:
    output:
        d=directory(graphpp_folder),
        exe_i=os.path.join(graphpp_folder, "graphpp"),
    conda: "../envs/compilation.yml"
    shell:
        """
        git clone https://github.com/ldenti/graphpp.git {output.d}
        cd {output.d}
        git apply ../../../patch/main.patch
        mkdir build ; cd build
        cmake ..
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
