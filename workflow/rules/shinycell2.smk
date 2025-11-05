from os.path import join, exists

seurat_object       = config['seurat_object']
run_dir             = config['run_dir']
project_title       = config['project_title']
marker_file         = config['marker_file'] if exists(config['marker_file']) else []
cluster_labels      = config['cluster_labels']
rmmeta              = config['rmmeta']
defaultreduction    = config['defaultreduction']
max_levels          = config['max_levels']
assaytouse          = config['assaytouse']
tmpdir              = config['tmpdir']


rule all:
    input:
        join(run_dir, "shinyapp.tar.gz")


rule shinycodes:
    input: 
        seurat_object       = seurat_object,
        marker_list         = marker_file,
        files_crumb         = join(run_dir, "wd", "sc1meta.rds"),
    output: join(run_dir, "wd", "server.R")
    container: "docker://rroutsong/shinycell2_builder:latest"
    params:
        rname = "shinycodes",
        project_title = project_title,
        cluster_labels = cluster_labels,
        rmmeta = rmmeta,
        defaultreduction = defaultreduction,
        max_levels = max_levels,
        assaytouse = assaytouse,
        outdir = join(run_dir, "wd"),
        marker_flag = lambda w, input: f"--markers {input.marker_list} " if exists(input.marker_list) else "",
    shell:
        """
        build_shinycell.R \\
            -j {input.seurat_object} {params.marker_flag}\\
            --proj {params.project_title} \\
            --cluster_labels {params.cluster_labels} \\
            --rmmeta {params.rmmeta} \\
            --defred {params.defaultreduction} \\
            -l {params.max_levels} \\
            -a {params.assaytouse} \\
            -o {params.outdir} \\
            --codesonly
        """

rule shinyfiles:
    input:
        seurat_object       = seurat_object,
        marker_list         = marker_file
    output:
        out_crumb           = join(run_dir, "wd", "sc1meta.rds"),
    container: "docker://rroutsong/shinycell2_builder:latest"
    params:
        rname = "shinyfiles",
        project_title = project_title,
        cluster_labels = cluster_labels,
        rmmeta = rmmeta,
        defaultreduction = defaultreduction,
        max_levels = max_levels,
        assaytouse = assaytouse,
        wd = join(run_dir, "wd"),
        marker_flag = lambda w, input: f"--markers {input.marker_list} " if exists(input.marker_list) else "",
        tmpdir = tmpdir
    shell:
        """
        build_shinycell.R \\
            -j {input.seurat_object} {params.marker_flag}\\
            --proj {params.project_title} \\
            --cluster_labels {params.cluster_labels} \\
            --rmmeta {params.rmmeta} \\
            --defred {params.defaultreduction} \\
            -l {params.max_levels} \\
            -a {assaytouse} \\
            -o {params.wd} \\
            --filesonly
        """


rule shinytar:
    input: 
        files_crumb         = join(run_dir, "wd", "sc1meta.rds"),
        codes_crumb         = join(run_dir, "wd", "server.R")
    output: join(run_dir, "shinyapp.tar.gz")
    params:
        rname = "shinytar",
        tmpdir = tmpdir,
        wd_dir = join(run_dir, "wd"),
    threads: 30
    shell:
        """
        cd {params.wd_dir}
        tar -cf - . | pigz -9 -p {threads} > {output}
        """