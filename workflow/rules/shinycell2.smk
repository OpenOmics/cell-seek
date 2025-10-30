from os.path import join

seurat_object       = config['seurat_object']
run_dir             = config['run_dir']
project_title       = config['project_title']
marker_file         = config['marker_file']
cluster_labels      = config['cluster_labels']
rmmeta              = config['rmmeta']
defaultreduction    = config['defaultreduction']
max_levels          = config['max_levels']
assaytouse          = config['assaytouse']
tmpdir              = config['tmpdir']


rule shinycodes:
    input: seurat_object
    output: directory(temp(join(run_dir, "codes")))
    container: "docker://rroutsong/shinycell2_builder:latest"
    params:
        project_title = project_title,
        marker_file = marker_file,
        cluster_labels = cluster_labels,
        rmmeta = rmmeta,
        defaultreduction = defaultreduction,
        max_levels = max_levels,
        assaytouse = assaytouse
    shell:
        """
        build_shinycell.R \\
            -j {input} \\
            --proj {params.project_title} \\
            --markers {params.marker_file} \\
            --cluster_labels {params.cluster_labels} \\
            --rmmeta {params.rmmeta} \\
            --defred {params.defaultreduction} \\
            -l {params.max_levels} \\
            -a {assaytouse} \\
            --codesonly
        mv shinyApp/* {output}
        """

rule shinycodes:
    input: seurat_object
    output: temp(directory(join(run_dir, "files")))
    container: "docker://rroutsong/shinycell2_builder:latest"
    params:
        project_title = project_title,
        marker_file = marker_file,
        cluster_labels = cluster_labels,
        rmmeta = rmmeta,
        defaultreduction = defaultreduction,
        max_levels = max_levels,
        assaytouse = assaytouse,    
        tmpdir = tmpdir
    shell:
        """
        build_shinycell.R \\
            -j {input} \\
            --proj {params.project_title} \\
            --markers {params.marker_file} \\
            --cluster_labels {params.cluster_labels} \\
            --rmmeta {params.rmmeta} \\
            --defred {params.defaultreduction} \\
            -l {params.max_levels} \\
            -a {assaytouse} \\
            --filesonly
        mv shinyApp/* {output}
        """


rule shinytar:
    input: 
        files = join(run_dir, "files")
        codes = join(run_dir, "codes")
    output: join(run_dir, "shinyapp.tar.gz")
    params:
        tmpdir = tmpdir
    threads: 30
    shell:
        """
        cp -r {input.files}/* {params.tmpdir}
        cp -r {input.codes}/* {params.tmpdir}
        pigz -9 -p {threads} -c -f {params.tmpdir}/* > {output}
        """