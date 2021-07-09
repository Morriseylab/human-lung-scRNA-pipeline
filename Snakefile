
############################## Lung Mapping #######################################################
checkpoint preprocess:
    input:
        script = "scripts/preprocess/{dataset}.R",
        log = "logs/download_{dataset}.log"
    params:
        dset = "{dataset}"
    output:
        directory("seurat_objects/unmapped/{dataset}"),
    container:
        "docker://satijalab/seurat:4.0.1"
    shell:
        """
        mkdir -p seurat_objects/unmapped/{wildcards.dataset}
        Rscript {input.script} {params.dset} {output} > logs/preprocess_{wildcards.dataset}.Rout 2>&1
        """

rule map:
    input:
        script = "scripts/map.R",
        ref = "reference/lung/ref.Rds",
        idx = "reference/lung/idx.annoy",
        query = "seurat_objects/unmapped/{dataset}/{donor}.rds"
    params:
        ref_loc = "reference/lung",
        xfer = ["annotation.l1", "annotation.l2"]
    output:
        "seurat_objects/mapped/{dataset}/{donor}.rds"
    container:
        "docker://satijalab/azimuth:0.3.1"
    shell:
        """
        mkdir -p seurat_objects/mapped/
        Rscript {input.script} {input.query} {params.ref_loc} {params.xfer} {output} > logs/map_{wildcards.dataset}_{wildcards.donor}.Rout 2>&1
        """
