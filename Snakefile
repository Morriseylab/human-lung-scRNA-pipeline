rule map:
    input:
        script = "scripts/PreProcessandMap.R",
        ref = "/home/bapoorva/Morrisey/Maria/lungMAP/Periph_5sample_Intergrate/Seuratv4_Anno.RDS"
    params:
        org = "human",
        input10x = "test/{donor}/{library}/filtered",
        outdir="test/{donor}/{library}/seurat/"
    output:
        "{donor}/{library}/seurat/seurat_refAnno.RDS"
    shell:
        """
        mkdir -p test/{wildcards.donor}/{wildcards.library}/seurat/
        Rscript {input.script} {wildcards.donor} {params.org} {params.input10x} {params.outdir} {input.ref}
        """
