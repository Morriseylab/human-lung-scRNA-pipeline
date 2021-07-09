rule map:
    input:
        script = "scripts/PreProcessandMap.R",
        ref = "Periph_5sample_Intergrate/Seuratv4_Anno.RDS",
        query = "{donor}/{library}/"
    params:
        org = "human"
    output:
        "{donor}/{library}/seurat/seurat_refanno.rds"
    shell:
        """
        mkdir -p {donor}/{library}/seurat/
        Rscript {input.script} {donor} {params.org} {params.xfer} {output} > logs/map_{wildcards.dataset}_{wildcards.donor}.Rout 2>&1
        """
