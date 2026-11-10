from os.path import join, exists
from textwrap import dedent

seurat_object       = config['seurat_object']
run_dir             = config['run_dir']
project_title       = config['project_title'].replace("\"", "")
marker_file         = config['marker_file'] if exists(config['marker_file']) else []
cluster_labels      = config['cluster_labels']
rmmeta              = config['rmmeta']
defaultreduction    = config['defaultreduction']
max_levels          = config['max_levels']
assaytouse          = config['assaytouse']
tmpdir              = config['tmpdir']


rule all:
    input:
        join(run_dir, "cellsightapp.tar.gz")


rule cellsightcodes:
    input:
        seurat_object       = seurat_object,
        marker_list         = marker_file,
        files_crumb         = join(run_dir, "wd", "sc1conf.qs2.rds"),
    output:
        server_out          = join(run_dir, "wd", "server.R"),
        ui_out              = join(run_dir, "wd", "ui.R"),
        obj_info            = join(run_dir, "wd", "seurat_inspector.info")
    container: config['images']['cellsight_builder']
    params:
        rname = "cellsightcodes",
        rmmeta = rmmeta,
        max_levels = max_levels,
        assaytouse = assaytouse,
        outdir = join(run_dir, "wd"),
        title = f'"{project_title}"',
        defred_flag = f"--defred {defaultreduction} " if defaultreduction else "",
        cluster_flag = f"--cluster_labels {cluster_labels} " if cluster_labels else "",
        marker_flag = lambda w, input: f"--markers {input.marker_list} " if input.marker_list and exists(input.marker_list) else ""
    shell:
        dedent("""
        seurat_inspector.R {input.seurat_object} > {output.obj_info}
        build_cellsight.R \\
            -j {input.seurat_object} {params.marker_flag} \\
            --proj {params.title} {params.cluster_flag} \\
            --rmmeta {params.rmmeta} {params.defred_flag} \\
            -l {params.max_levels} \\
            -a {params.assaytouse} \\
            -o {params.outdir} \\
            --codesonly
        """)

rule cellsightfiles:
    input:
        seurat_object       = seurat_object,
        marker_list         = marker_file
    output:
        conf_crumb          = join(run_dir, "wd", "sc1conf.qs2.rds"),
        meta_crumb          = join(run_dir, "wd", "sc1meta.qs2.rds"),
        gene_crumb          = join(run_dir, "wd", "sc1gene.qs2.rds"),
        dimr_crumb          = join(run_dir, "wd", "sc1dimr.qs2.rds"),
        def_crumb           = join(run_dir, "wd", "sc1def.qs2.rds"),
    container: config['images']['cellsight_builder']
    params:
        rname = "cellsightfiles",
        rmmeta = rmmeta,
        max_levels = max_levels,
        assaytouse = assaytouse,
        wd = join(run_dir, "wd"),
        title = f'"{project_title}"',
        defred_flag = f"--defred {defaultreduction} " if defaultreduction else "",
        cluster_flag = f"--cluster_labels {cluster_labels} " if cluster_labels else "",
        marker_flag = lambda w, input: f"--markers {input.marker_list} " if input.marker_list and exists(input.marker_list) else "",
    shell:
        dedent("""
        build_cellsight.R \\
            -j {input.seurat_object} {params.marker_flag} \\
            --proj {params.title} {params.cluster_flag} \\
            --rmmeta {params.rmmeta} {params.defred_flag} \\
            -l {params.max_levels} \\
            -a {assaytouse} \\
            -o {params.wd} \\
            --filesonly
        """)


rule cellsighttar:
    input:
        files_crumb         = join(run_dir, "wd", "sc1conf.qs2.rds"),
        codes_crumb         = join(run_dir, "wd", "server.R")
    output: join(run_dir, "cellsightapp.tar.gz")
    params:
        rname = "cellsighttar",
        tmpdir = tmpdir,
        wd_dir = join(run_dir, "wd"),
    threads: 30
    shell:
        dedent("""
        cd {params.wd_dir}
        tar -cf - . | pigz -9 -p {threads} > {output}
        """)
