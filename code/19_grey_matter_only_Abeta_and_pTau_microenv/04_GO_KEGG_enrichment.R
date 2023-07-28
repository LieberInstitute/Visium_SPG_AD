library("here")
library("sessioninfo")
library("spatialLIBD")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")

## Plot output directory
dir_plots <- here::here(
    "plots",
    "17_grey_matter_only_Abeta_microenv",
    "wholegenome"
)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

## Locate data directory
dir_rdata <- here::here(
    "code",
    "18_deploy_app_wholegenome_Abeta_microenv"
)

load(file.path(dir_rdata, "Visium_SPG_AD_modeling_results.Rdata"),
    verbose = TRUE
)
sce_pseudo <-
    readRDS(file.path(dir_rdata, "sce_pseudo_pathology_wholegenome.rds"))


## For sig_genes_extract_all() to work
sce_pseudo$spatialLIBD <- sce_pseudo$path_groups
sig_genes <- sig_genes_extract_all(
    n = nrow(sce_pseudo),
    modeling_results = modeling_results,
    sce_layer = sce_pseudo
)

subset_sig <- subset(sig_genes, model_type == "enrichment")
tapply(subset_sig$fdr, subset_sig$test, summary)
# $Ab_env
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000934 0.2086530 0.4394944 0.4690056 0.7245988 0.9999039
#
# $both
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.001212 0.923924 0.964169 0.912828 0.979182 0.999964
#
# $n_both
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.1718  0.9994  0.9994  0.9895  0.9994  0.9998
#
# $n_pTau
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       1       1       1       1       1       1
#
# $none
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02605 0.98223 0.98223 0.97743 0.98760 0.99999
#
# $pTau
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.9996  0.9996  0.9996  0.9996  0.9996  1.0000

dim(subset(subset_sig, fdr < 0.05))
# [1] 287  13

with(subset(subset_sig, fdr < 0.05), addmargins(table(sign(logFC), test)))
#    test
#     Ab_env both none Sum
# -1     265    9    2 276
# 1       10    0    1  11
# Sum    275    9    3 287

with(subset(subset_sig, fdr < 0.10), addmargins(table(sign(logFC), test)))
#    test
#     Ab_env both none Sum
# -1     901    9    2 912
# 1       76    1    1  78
# Sum    977   10    3 990

with(subset(subset_sig, fdr < 0.20), addmargins(table(sign(logFC), test)))
#    test
#     Ab_env both n_both none  Sum
# -1    1881   42      1    2 1926
# 1      246    5      0    2  253
# Sum   2127   47      1    4 2179

## Try with FDR < 0.10
fdr_genes <- subset(subset_sig, fdr < 0.05 & test == "Ab_env")
sigGene <- lapply(split(fdr_genes$ensembl, sign(fdr_genes$logFC)), function(ensembl) {
    bitr(ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
})
lengths(sigGene)
geneUniverse <- unique(bitr(unique(subset_sig$ensembl), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID)
length(geneUniverse)

## The following GO code was adapted from
## https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/05_GO_enrichment/01_GO_enrichment.R

## Run GO and KEGG enrichment analysis
go <- compareCluster(sigGene,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = "org.Hs.eg.db",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.3,
    qvalueCutoff = 0.3,
    readable = TRUE
)

kegg <- compareCluster(sigGene,
    fun = "enrichKEGG",
    universe = geneUniverse,
    organism = "hsa",
    pvalueCutoff = 0.3,
    qvalueCutoff = 0.3
)

## Save Rdata with gene ontology enrichment results
save(go, kegg, file = here::here(
    "processed-data",
    "17_grey_matter_only_Abeta_microenv",
    "wholegenome",
    "gene_ontology_enrichment_objects.Rdata"
))

## Plot BP, CC and MF
plot_go <- function(ont, title_p, path, filename) {
    dotplot_1 <- ggplot(filter(go, ONTOLOGY == ont), aes(Cluster, Description)) +
        theme_bw() +
        geom_point(aes(color = p.adjust, size = Count)) +
        scale_color_gradientn(
            colours = c("#f7ca64", "#46bac2", "#7e62a3"),
            trans = "log10",
            guide = guide_colorbar(reverse = TRUE, order = 1)
        ) +
        scale_size_continuous(range = c(2, 10)) +
        xlab("Cluster") +
        ylab("") +
        ggtitle(title_p)

    ggsave(filename = filename, path = path, dotplot_1, height = 6, width = 5)
}

plot_go(ont = "BP", title_p = "Biological Process", filename = "GOenrichment_BP.pdf", path = dir_plots)
# plot_go(ont = "CC", title_p = "Cellular Component", filename = "GOenrichment_CC.pdf", path = dir_plots) ## Has no results at FDR < 5%, only at < 20%
plot_go(ont = "MF", title_p = "Molecular Function", filename = "GOenrichment_MF.pdf", path = dir_plots)

## Plot KEGG
dotplot_1 <- ggplot(kegg, aes(Cluster, Description)) +
    theme_bw() +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
        colours = c("#f7ca64", "#46bac2", "#7e62a3"),
        trans = "log10",
        guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    xlab("Cluster") +
    ylab("") +
    ggtitle("KEGG")

ggsave(filename = "KEGGenrichment.pdf", path = dir_plots, dotplot_1, height = 6, width = 5)

## Print some results with FDR < 0.20
subset(as.data.frame(go)[, -10], p.adjust < 0.3)
#      Cluster ONTOLOGY         ID                                                 Description GeneRatio  BgRatio       pvalue   p.adjust     qvalue Count
# 4484      -1       CC GO:0030672                                   synaptic vesicle membrane   37/1750  95/8375 4.095356e-05 0.01255227 0.01255227    37
# 4485      -1       CC GO:0099501                                   exocytic vesicle membrane   37/1750  95/8375 4.095356e-05 0.01255227 0.01255227    37
# 4486      -1       CC GO:0008021                                            synaptic vesicle   54/1750 166/8375 2.744457e-04 0.05607841 0.05607841    54
# 4487      -1       CC GO:0070382                                            exocytic vesicle   55/1750 172/8375 3.929452e-04 0.06021886 0.06021886    55
# 4488      -1       CC GO:0030285             integral component of synaptic vesicle membrane   13/1750  25/8375 5.662937e-04 0.06942761 0.06942761    13
# 4489      -1       CC GO:0098793                                                  presynapse  106/1750 384/8375 8.083867e-04 0.08259017 0.08259017   106
# 4490      -1       CC GO:0030658                                  transport vesicle membrane   52/1750 166/8375 9.524999e-04 0.08341178 0.08341178    52
# 4491      -1       CC GO:0098563            intrinsic component of synaptic vesicle membrane   16/1750  37/8375 1.759721e-03 0.13483863 0.13483863    16
# 4492      -1       CC GO:0046540                                U4/U6 x U5 tri-snRNP complex   14/1750  33/8375 4.158867e-03 0.23775851 0.23775851    14
# 4493      -1       CC GO:0097526                              spliceosomal tri-snRNP complex   14/1750  33/8375 4.158867e-03 0.23775851 0.23775851    14
# 4494      -1       CC GO:0030136                                     clathrin-coated vesicle   40/1750 129/8375 4.266466e-03 0.23775851 0.23775851    40
# 4495      -1       CC GO:0005759                                        mitochondrial matrix   91/1750 342/8375 5.786139e-03 0.29557527 0.29557527    91
# 5886       1       BP GO:0050878                             regulation of body fluid levels    14/226 149/8068 6.964129e-05 0.16985950 0.16262819    14
# 5887       1       BP GO:0007599                                                  hemostasis    10/226  88/8068 1.624719e-04 0.16985950 0.16262819    10
# 5888       1       BP GO:0030198                           extracellular matrix organization    10/226  96/8068 3.340234e-04 0.16985950 0.16262819    10
# 5889       1       BP GO:0043062                        extracellular structure organization    10/226  96/8068 3.340234e-04 0.16985950 0.16262819    10
# 5890       1       BP GO:0045229               external encapsulating structure organization    10/226  96/8068 3.340234e-04 0.16985950 0.16262819    10
# 5891       1       BP GO:0002821             positive regulation of adaptive immune response     6/226  35/8068 3.713828e-04 0.16985950 0.16262819     6
# 5892       1       BP GO:0031349                     positive regulation of defense response    10/226 100/8068 4.649787e-04 0.16985950 0.16262819    10
# 5893       1       BP GO:0045945  positive regulation of transcription by RNA polymerase III     4/226  14/8068 4.810539e-04 0.16985950 0.16262819     4
# 5894       1       BP GO:0002819                      regulation of adaptive immune response     7/226  52/8068 5.638703e-04 0.16985950 0.16262819     7
# 5895       1       BP GO:0007596                                           blood coagulation     9/226  86/8068 6.389702e-04 0.16985950 0.16262819     9
# 5896       1       BP GO:0050778                      positive regulation of immune response    14/226 185/8068 6.701945e-04 0.16985950 0.16262819    14
# 5897       1       BP GO:0050817                                                 coagulation     9/226  87/8068 6.956703e-04 0.16985950 0.16262819     9
# 5898       1       BP GO:1903053             regulation of extracellular matrix organization     5/226  27/8068 8.045431e-04 0.18133163 0.17361192     5
# 5899       1       BP GO:0030168                                         platelet activation     7/226  57/8068 9.898449e-04 0.20396692 0.19528358     7
# 5900       1       BP GO:0045088                        regulation of innate immune response     9/226  92/8068 1.044199e-03 0.20396692 0.19528358     9
# 5901       1       BP GO:0018958                phenol-containing compound metabolic process     6/226  45/8068 1.475964e-03 0.24736439 0.23683352     6
# 5902       1       BP GO:0006536                                 glutamate metabolic process     4/226  19/8068 1.668367e-03 0.24736439 0.23683352     4
# 5903       1       BP GO:0009410                             response to xenobiotic stimulus    13/226 181/8068 1.672868e-03 0.24736439 0.23683352    13
# 5904       1       BP GO:0001775                                             cell activation    23/226 424/8068 1.729105e-03 0.24736439 0.23683352    23
# 5905       1       BP GO:0002526                                 acute inflammatory response     5/226  32/8068 1.790948e-03 0.24736439 0.23683352     5
# 5906       1       BP GO:0036473                  cell death in response to oxidative stress     7/226  63/8068 1.797484e-03 0.24736439 0.23683352     7
# 5907       1       BP GO:0050864                             regulation of B cell activation     6/226  47/8068 1.857344e-03 0.24736439 0.23683352     6
# 5908       1       BP GO:0001570                                              vasculogenesis     5/226  33/8068 2.063341e-03 0.26285168 0.25166148     5
# 5909       1       BP GO:0002684                positive regulation of immune system process    19/226 329/8068 2.160455e-03 0.26366999 0.25244495    19
# 5910       1       BP GO:0010574 regulation of vascular endothelial growth factor production     3/226  10/8068 2.249744e-03 0.26366999 0.25244495     3
# 5911       1       BP GO:0070527                                        platelet aggregation     5/226  34/8068 2.364779e-03 0.26649244 0.25514725     5
# 9213       1       MF GO:0005201                 extracellular matrix structural constituent     6/228  36/8223 4.126048e-04 0.19268645 0.19023254     6

## Note how only 2 have a p.adjust < 0.05: the synaptic / exocytic vesicle membrane CC

subset(as.data.frame(kegg)[, -9], p.adjust < 0.3)
#     Cluster       ID                                 Description GeneRatio  BgRatio       pvalue  p.adjust    qvalue Count
# 1        -1 hsa03018                             RNA degradation    24/793  61/3826 0.0006267409 0.1855153 0.1855153    24
# 297       1 hsa04510                              Focal adhesion    11/114 120/3826 0.0007619300 0.1306763 0.1285552    11
# 298       1 hsa04512                    ECM-receptor interaction     5/114  28/3826 0.0012212740 0.1306763 0.1285552     5
# 299       1 hsa04610         Complement and coagulation cascades     4/114  19/3826 0.0020498125 0.1462200 0.1438465     4
# 300       1 hsa00250 Alanine, aspartate and glutamate metabolism     4/114  24/3826 0.0050098690 0.2680280 0.2636773     4
# 301       1 hsa04810            Regulation of actin cytoskeleton    10/114 137/3826 0.0070715557 0.2760183 0.2715379    10
# 302       1 hsa04727                           GABAergic synapse     6/114  62/3826 0.0096913448 0.2760183 0.2715379     6
# 303       1 hsa00650                        Butanoate metabolism     3/114  16/3826 0.0108687384 0.2760183 0.2715379     3
# 304       1 hsa04940                    Type I diabetes mellitus     3/114  16/3826 0.0108687384 0.2760183 0.2715379     3
# 305       1 hsa04145                                   Phagosome     7/114  84/3826 0.0118344470 0.2760183 0.2715379     7
# 306       1 hsa05130       Pathogenic Escherichia coli infection     8/114 106/3826 0.0128980517 0.2760183 0.2715379     8

## Print some results with FDR < 0.05
options(max.print = 10000)
subset(as.data.frame(go), p.adjust < 0.2)
#     Cluster ONTOLOGY         ID                                                                                                  Description GeneRatio  BgRatio      pvalue  p.adjust    qvalue Count                  geneID
# 1         1       BP GO:0006690                                                                                  icosanoid metabolic process      2/10  41/8068 0.001105028 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 2         1       BP GO:0097193                                                                        intrinsic apoptotic signaling pathway      3/10 187/8068 0.001303985 0.1419006 0.1005438     3 ZNF385A/TNFRSF1A/CYP1B1
# 3         1       BP GO:0033559                                                                     unsaturated fatty acid metabolic process      2/10  47/8068 0.001450978 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 4         1       BP GO:0007259                                                                      receptor signaling pathway via JAK-STAT      2/10  55/8068 0.001982726 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 5         1       BP GO:0097696                                                                          receptor signaling pathway via STAT      2/10  57/8068 0.002128111 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 6         1       BP GO:0008630                                              intrinsic apoptotic signaling pathway in response to DNA damage      2/10  61/8068 0.002433676 0.1419006 0.1005438     2        ZNF385A/TNFRSF1A
# 7         1       BP GO:0071356                                                                   cellular response to tumor necrosis factor      2/10  94/8068 0.005687306 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 8         1       BP GO:0031348                                                                      negative regulation of defense response      2/10  95/8068 0.005805769 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 9         1       BP GO:0030198                                                                            extracellular matrix organization      2/10  96/8068 0.005925372 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 10        1       BP GO:0043062                                                                         extracellular structure organization      2/10  96/8068 0.005925372 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 11        1       BP GO:0045229                                                                external encapsulating structure organization      2/10  96/8068 0.005925372 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 12        1       BP GO:0031349                                                                      positive regulation of defense response      2/10 100/8068 0.006415138 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 13        1       BP GO:0034612                                                                            response to tumor necrosis factor      2/10 107/8068 0.007315565 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 14        1       BP GO:0097190                                                                                  apoptotic signaling pathway      3/10 347/8068 0.007552162 0.1419006 0.1005438     3 ZNF385A/TNFRSF1A/CYP1B1
# 15        1       BP GO:0009636                                                                                  response to toxic substance      2/10 114/8068 0.008270445 0.1419006 0.1005438     2             MT1F/CYP1B1
# 16        1       BP GO:0043122                                                            regulation of I-kappaB kinase/NF-kappaB signaling      2/10 123/8068 0.009576849 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 17        1       BP GO:0007249                                                                          I-kappaB kinase/NF-kappaB signaling      2/10 139/8068 0.012112969 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 18        1       BP GO:0008210                                                                                   estrogen metabolic process      1/10  10/8068 0.012332583 0.1419006 0.1005438     1                  CYP1B1
# 19        1       BP GO:0008298                                                                              intracellular mRNA localization      1/10  10/8068 0.012332583 0.1419006 0.1005438     1                 ZNF385A
# 20        1       BP GO:0010574                                                  regulation of vascular endothelial growth factor production      1/10  10/8068 0.012332583 0.1419006 0.1005438     1                  CYP1B1
# 21        1       BP GO:0034698                                                                                     response to gonadotropin      1/10  10/8068 0.012332583 0.1419006 0.1005438     1                  CYP1B1
# 22        1       BP GO:0042772                                          DNA damage response, signal transduction resulting in transcription      1/10  10/8068 0.012332583 0.1419006 0.1005438     1                 ZNF385A
# 23        1       BP GO:0043517                        positive regulation of DNA damage response, signal transduction by p53 class mediator      1/10  10/8068 0.012332583 0.1419006 0.1005438     1                 ZNF385A
# 24        1       BP GO:0001522                                                                                      pseudouridine synthesis      1/10  11/8068 0.013558281 0.1419006 0.1005438     1                   PUS7L
# 25        1       BP GO:0010573                                                                vascular endothelial growth factor production      1/10  11/8068 0.013558281 0.1419006 0.1005438     1                  CYP1B1
# 26        1       BP GO:0044849                                                                                                estrous cycle      1/10  11/8068 0.013558281 0.1419006 0.1005438     1                  CYP1B1
# 27        1       BP GO:0046427                                               positive regulation of receptor signaling pathway via JAK-STAT      1/10  11/8068 0.013558281 0.1419006 0.1005438     1                  CYP1B1
# 28        1       BP GO:0030220                                                                                           platelet formation      1/10  12/8068 0.014782610 0.1419006 0.1005438     1                 ZNF385A
# 29        1       BP GO:0035855                                                                                    megakaryocyte development      1/10  12/8068 0.014782610 0.1419006 0.1005438     1                 ZNF385A
# 30        1       BP GO:0061687                                                                         detoxification of inorganic compound      1/10  12/8068 0.014782610 0.1419006 0.1005438     1                    MT1F
# 31        1       BP GO:0097501                                                                                 stress response to metal ion      1/10  12/8068 0.014782610 0.1419006 0.1005438     1                    MT1F
# 32        1       BP GO:1902166 negative regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator      1/10  12/8068 0.014782610 0.1419006 0.1005438     1                 ZNF385A
# 33        1       BP GO:1904894                                                   positive regulation of receptor signaling pathway via STAT      1/10  12/8068 0.014782610 0.1419006 0.1005438     1                  CYP1B1
# 34        1       BP GO:0001819                                                                   positive regulation of cytokine production      2/10 160/8068 0.015842886 0.1419006 0.1005438     2            RIOK3/CYP1B1
# 35        1       BP GO:0036344                                                                                       platelet morphogenesis      1/10  13/8068 0.016005571 0.1419006 0.1005438     1                 ZNF385A
# 36        1       BP GO:0039531                       regulation of viral-induced cytoplasmic pattern recognition receptor signaling pathway      1/10  13/8068 0.016005571 0.1419006 0.1005438     1                   RIOK3
# 37        1       BP GO:0042572                                                                                    retinol metabolic process      1/10  13/8068 0.016005571 0.1419006 0.1005438     1                  CYP1B1
# 38        1       BP GO:0071294                                                                                cellular response to zinc ion      1/10  13/8068 0.016005571 0.1419006 0.1005438     1                    MT1F
# 39        1       BP GO:0071359                                                                                   cellular response to dsRNA      1/10  13/8068 0.016005571 0.1419006 0.1005438     1                   RIOK3
# 40        1       BP GO:0032103                                                         positive regulation of response to external stimulus      2/10 164/8068 0.016603444 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 41        1       BP GO:0032102                                                         negative regulation of response to external stimulus      2/10 167/8068 0.017184165 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 42        1       BP GO:0003176                                                                                     aortic valve development      1/10  14/8068 0.017227165 0.1419006 0.1005438     1                TNFRSF1A
# 43        1       BP GO:0042531                                              positive regulation of tyrosine phosphorylation of STAT protein      1/10  14/8068 0.017227165 0.1419006 0.1005438     1                TNFRSF1A
# 44        1       BP GO:1902165          regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator      1/10  14/8068 0.017227165 0.1419006 0.1005438     1                 ZNF385A
# 45        1       BP GO:1902742                                                                    apoptotic process involved in development      1/10  14/8068 0.017227165 0.1419006 0.1005438     1                TNFRSF1A
# 46        1       BP GO:0019369                                                                           arachidonic acid metabolic process      1/10  15/8068 0.018447395 0.1419006 0.1005438     1                  CYP1B1
# 47        1       BP GO:0032570                                                                                     response to progesterone      1/10  15/8068 0.018447395 0.1419006 0.1005438     1                  CYP1B1
# 48        1       BP GO:0061684                                                                                 chaperone-mediated autophagy      1/10  15/8068 0.018447395 0.1419006 0.1005438     1                    CTSA
# 49        1       BP GO:1905314                                                                                 semi-lunar valve development      1/10  15/8068 0.018447395 0.1419006 0.1005438     1                TNFRSF1A
# 50        1       BP GO:0014911                                                          positive regulation of smooth muscle cell migration      1/10  16/8068 0.019666261 0.1419006 0.1005438     1                  CYP1B1
# 51        1       BP GO:0030325                                                                                    adrenal gland development      1/10  16/8068 0.019666261 0.1419006 0.1005438     1                  CYP1B1
# 52        1       BP GO:0050687                                                             negative regulation of defense response to virus      1/10  16/8068 0.019666261 0.1419006 0.1005438     1                   RIOK3
# 53        1       BP GO:0006974                                                                     cellular response to DNA damage stimulus      3/10 498/8068 0.020241673 0.1419006 0.1005438     3  SUSD6/ZNF385A/TNFRSF1A
# 54        1       BP GO:0030199                                                                                 collagen fibril organization      1/10  17/8068 0.020883764 0.1419006 0.1005438     1                  CYP1B1
# 55        1       BP GO:0032728                                                            positive regulation of interferon-beta production      1/10  17/8068 0.020883764 0.1419006 0.1005438     1                   RIOK3
# 56        1       BP GO:0045601                                                               regulation of endothelial cell differentiation      1/10  17/8068 0.020883764 0.1419006 0.1005438     1                TNFRSF1A
# 57        1       BP GO:0046466                                                                             membrane lipid catabolic process      1/10  17/8068 0.020883764 0.1419006 0.1005438     1                  CYP1B1
# 58        1       BP GO:0006631                                                                                 fatty acid metabolic process      2/10 187/8068 0.021276768 0.1419006 0.1005438     2         TNFRSF1A/CYP1B1
# 59        1       BP GO:0010614                                                            negative regulation of cardiac muscle hypertrophy      1/10  18/8068 0.022099907 0.1419006 0.1005438     1                TNFRSF1A
# 60        1       BP GO:0033628                                                             regulation of cell adhesion mediated by integrin      1/10  18/8068 0.022099907 0.1419006 0.1005438     1                  CYP1B1
# 61        1       BP GO:0039528                              cytoplasmic pattern recognition receptor signaling pathway in response to virus      1/10  18/8068 0.022099907 0.1419006 0.1005438     1                   RIOK3
# 62        1       BP GO:0071280                                                                              cellular response to copper ion      1/10  18/8068 0.022099907 0.1419006 0.1005438     1                    MT1F
# 63        1       BP GO:1902254                           negative regulation of intrinsic apoptotic signaling pathway by p53 class mediator      1/10  18/8068 0.022099907 0.1419006 0.1005438     1                 ZNF385A
# 64        1       BP GO:0014741                                                                    negative regulation of muscle hypertrophy      1/10  19/8068 0.023314689 0.1419006 0.1005438     1                TNFRSF1A
# 65        1       BP GO:0016556                                                                                            mRNA modification      1/10  20/8068 0.024528114 0.1419006 0.1005438     1                   PUS7L
# 66        1       BP GO:0019748                                                                                  secondary metabolic process      1/10  20/8068 0.024528114 0.1419006 0.1005438     1                  CYP1B1
# 67        1       BP GO:0042509                                                       regulation of tyrosine phosphorylation of STAT protein      1/10  20/8068 0.024528114 0.1419006 0.1005438     1                TNFRSF1A
# 68        1       BP GO:1902230                       negative regulation of intrinsic apoptotic signaling pathway in response to DNA damage      1/10  20/8068 0.024528114 0.1419006 0.1005438     1                 ZNF385A
# 69        1       BP GO:0071276                                                                             cellular response to cadmium ion      1/10  21/8068 0.025740181 0.1419006 0.1005438     1                    MT1F
# 70        1       BP GO:0098751                                                                                        bone cell development      1/10  21/8068 0.025740181 0.1419006 0.1005438     1                 ZNF385A
# 71        1       BP GO:0007260                                                                     tyrosine phosphorylation of STAT protein      1/10  22/8068 0.026950893 0.1419006 0.1005438     1                TNFRSF1A
# 72        1       BP GO:0032608                                                                                   interferon-beta production      1/10  22/8068 0.026950893 0.1419006 0.1005438     1                   RIOK3
# 73        1       BP GO:0032648                                                                     regulation of interferon-beta production      1/10  22/8068 0.026950893 0.1419006 0.1005438     1                   RIOK3
# 74        1       BP GO:1901798                                             positive regulation of signal transduction by p53 class mediator      1/10  22/8068 0.026950893 0.1419006 0.1005438     1                 ZNF385A
# 75        1       BP GO:0006882                                                                                cellular zinc ion homeostasis      1/10  23/8068 0.028160251 0.1419006 0.1005438     1                    MT1F
# 76        1       BP GO:0043331                                                                                            response to dsRNA      1/10  23/8068 0.028160251 0.1419006 0.1005438     1                   RIOK3
# 77        1       BP GO:0071548                                                                                    response to dexamethasone      1/10  23/8068 0.028160251 0.1419006 0.1005438     1                  CYP1B1
# 78        1       BP GO:2000765                                                                        regulation of cytoplasmic translation      1/10  23/8068 0.028160251 0.1419006 0.1005438     1                 ZNF385A
# 79        1       BP GO:0006692                                                                                 prostanoid metabolic process      1/10  24/8068 0.029368255 0.1419006 0.1005438     1                TNFRSF1A
# 80        1       BP GO:0006693                                                                              prostaglandin metabolic process      1/10  24/8068 0.029368255 0.1419006 0.1005438     1                TNFRSF1A
# 81        1       BP GO:0043124                                                   negative regulation of I-kappaB kinase/NF-kappaB signaling      1/10  24/8068 0.029368255 0.1419006 0.1005438     1                   RIOK3
# 82        1       BP GO:0055069                                                                                         zinc ion homeostasis      1/10  24/8068 0.029368255 0.1419006 0.1005438     1                    MT1F
# 83        1       BP GO:1902253                                    regulation of intrinsic apoptotic signaling pathway by p53 class mediator      1/10  24/8068 0.029368255 0.1419006 0.1005438     1                 ZNF385A
# 84        1       BP GO:0001523                                                                                   retinoid metabolic process      1/10  25/8068 0.030574908 0.1419006 0.1005438     1                  CYP1B1
# 85        1       BP GO:0010043                                                                                         response to zinc ion      1/10  25/8068 0.030574908 0.1419006 0.1005438     1                    MT1F
# 86        1       BP GO:0046688                                                                                       response to copper ion      1/10  25/8068 0.030574908 0.1419006 0.1005438     1                    MT1F
# 87        1       BP GO:1901797                                             negative regulation of signal transduction by p53 class mediator      1/10  25/8068 0.030574908 0.1419006 0.1005438     1                 ZNF385A
# 88        1       BP GO:0010803                                               regulation of tumor necrosis factor-mediated signaling pathway      1/10  26/8068 0.031780211 0.1419006 0.1005438     1                TNFRSF1A
# 89        1       BP GO:0030219                                                                                megakaryocyte differentiation      1/10  26/8068 0.031780211 0.1419006 0.1005438     1                 ZNF385A
# 90        1       BP GO:0043516                                 regulation of DNA damage response, signal transduction by p53 class mediator      1/10  26/8068 0.031780211 0.1419006 0.1005438     1                 ZNF385A
# 91        1       BP GO:0046425                                                        regulation of receptor signaling pathway via JAK-STAT      1/10  26/8068 0.031780211 0.1419006 0.1005438     1                  CYP1B1
# 92        1       BP GO:0046685                                                                     response to arsenic-containing substance      1/10  26/8068 0.031780211 0.1419006 0.1005438     1                  CYP1B1
# 93        1       BP GO:0031347                                                                               regulation of defense response      2/10 233/8068 0.032074641 0.1419006 0.1005438     2          TNFRSF1A/RIOK3
# 94        1       BP GO:0016101                                                                                diterpenoid metabolic process      1/10  27/8068 0.032984165 0.1419006 0.1005438     1                  CYP1B1
# 95        1       BP GO:1902229                                regulation of intrinsic apoptotic signaling pathway in response to DNA damage      1/10  27/8068 0.032984165 0.1419006 0.1005438     1                 ZNF385A
# 96        1       BP GO:1903053                                                              regulation of extracellular matrix organization      1/10  27/8068 0.032984165 0.1419006 0.1005438     1                TNFRSF1A
# 97        1       BP GO:0061028                                                                         establishment of endothelial barrier      1/10  28/8068 0.034186772 0.1430534 0.1013606     1                TNFRSF1A
# 98        1       BP GO:1904892                                                            regulation of receptor signaling pathway via STAT      1/10  28/8068 0.034186772 0.1430534 0.1013606     1                  CYP1B1
# 99        1       BP GO:0003170                                                                                      heart valve development      1/10  29/8068 0.035388032 0.1430534 0.1013606     1                TNFRSF1A
# 100       1       BP GO:0032481                                                          positive regulation of type I interferon production      1/10  30/8068 0.036587947 0.1430534 0.1013606     1                   RIOK3
# 101       1       BP GO:0042771                        intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator      1/10  30/8068 0.036587947 0.1430534 0.1013606     1                 ZNF385A
# 102       1       BP GO:0071260                                                                     cellular response to mechanical stimulus      1/10  30/8068 0.036587947 0.1430534 0.1013606     1                TNFRSF1A
# 103       1       BP GO:0071320                                                                                    cellular response to cAMP      1/10  30/8068 0.036587947 0.1430534 0.1013606     1                  CYP1B1
# 104       1       BP GO:0002753                                                   cytoplasmic pattern recognition receptor signaling pathway      1/10  31/8068 0.037786519 0.1430534 0.1013606     1                   RIOK3
# 105       1       BP GO:0006081                                                                          cellular aldehyde metabolic process      1/10  31/8068 0.037786519 0.1430534 0.1013606     1                  CYP1B1
# 106       1       BP GO:0034308                                                                            primary alcohol metabolic process      1/10  31/8068 0.037786519 0.1430534 0.1013606     1                  CYP1B1
# 107       1       BP GO:0050688                                                                      regulation of defense response to virus      1/10  31/8068 0.037786519 0.1430534 0.1013606     1                   RIOK3
# 108       1       BP GO:0071385                                                                 cellular response to glucocorticoid stimulus      1/10  31/8068 0.037786519 0.1430534 0.1013606     1                  CYP1B1
# 109       1       BP GO:0043065                                                                     positive regulation of apoptotic process      2/10 258/8068 0.038697056 0.1430534 0.1013606     2         TNFRSF1A/CYP1B1
# 110       1       BP GO:0033627                                                                           cell adhesion mediated by integrin      1/10  32/8068 0.038983749 0.1430534 0.1013606     1                  CYP1B1
# 111       1       BP GO:0042698                                                                                              ovulation cycle      1/10  32/8068 0.038983749 0.1430534 0.1013606     1                  CYP1B1
# 112       1       BP GO:0043068                                                                 positive regulation of programmed cell death      2/10 263/8068 0.040081459 0.1430534 0.1013606     2         TNFRSF1A/CYP1B1
# 113       1       BP GO:0006721                                                                                  terpenoid metabolic process      1/10  33/8068 0.040179637 0.1430534 0.1013606     1                  CYP1B1
# 114       1       BP GO:0006805                                                                                 xenobiotic metabolic process      1/10  33/8068 0.040179637 0.1430534 0.1013606     1                  CYP1B1
# 115       1       BP GO:0046686                                                                                      response to cadmium ion      1/10  33/8068 0.040179637 0.1430534 0.1013606     1                    MT1F
# 116       1       BP GO:0071384                                                                 cellular response to corticosteroid stimulus      1/10  33/8068 0.040179637 0.1430534 0.1013606     1                  CYP1B1
# 117       1       BP GO:0006809                                                                            nitric oxide biosynthetic process      1/10  34/8068 0.041374187 0.1439091 0.1019669     1                  CYP1B1
# 118       1       BP GO:0001817                                                                            regulation of cytokine production      2/10 268/8068 0.041485269 0.1439091 0.1019669     2            RIOK3/CYP1B1
# 119       1       BP GO:0001816                                                                                          cytokine production      2/10 271/8068 0.042336770 0.1439091 0.1019669     2            RIOK3/CYP1B1
# 120       1       BP GO:0008631                                        intrinsic apoptotic signaling pathway in response to oxidative stress      1/10  35/8068 0.042567398 0.1439091 0.1019669     1                  CYP1B1
# 121       1       BP GO:0010611                                                                     regulation of cardiac muscle hypertrophy      1/10  36/8068 0.043759272 0.1439091 0.1019669     1                TNFRSF1A
# 122       1       BP GO:0014910                                                                   regulation of smooth muscle cell migration      1/10  36/8068 0.043759272 0.1439091 0.1019669     1                  CYP1B1
# 123       1       BP GO:2000379                                             positive regulation of reactive oxygen species metabolic process      1/10  36/8068 0.043759272 0.1439091 0.1019669     1                  CYP1B1
# 124       1       BP GO:0001885                                                                                 endothelial cell development      1/10  37/8068 0.044949811 0.1439091 0.1019669     1                TNFRSF1A
# 125       1       BP GO:0014743                                                                             regulation of muscle hypertrophy      1/10  37/8068 0.044949811 0.1439091 0.1019669     1                TNFRSF1A
# 126       1       BP GO:0030490                                                                                       maturation of SSU-rRNA      1/10  37/8068 0.044949811 0.1439091 0.1019669     1                   RIOK3
# 127       1       BP GO:0033363                                                                               secretory granule organization      1/10  37/8068 0.044949811 0.1439091 0.1019669     1                 ZNF385A
# 128       1       BP GO:0045600                                                              positive regulation of fat cell differentiation      1/10  37/8068 0.044949811 0.1439091 0.1019669     1                 ZNF385A
# 129       1       BP GO:0046209                                                                               nitric oxide metabolic process      1/10  37/8068 0.044949811 0.1439091 0.1019669     1                  CYP1B1
# 130       1       BP GO:2001057                                                                  reactive nitrogen species metabolic process      1/10  38/8068 0.046139015 0.1465801 0.1038595     1                  CYP1B1
# 131       1       BP GO:0034754                                                                           cellular hormone metabolic process      1/10  39/8068 0.047326887 0.1492061 0.1057202     1                  CYP1B1
# 132       1       BP GO:0006417                                                                                    regulation of translation      2/10 293/8068 0.048787099 0.1514968 0.1073432     2          ZNF385A/CYP1B1
# 133       1       BP GO:0032787                                                                        monocarboxylic acid metabolic process      2/10 293/8068 0.048787099 0.1514968 0.1073432     2         TNFRSF1A/CYP1B1
# 134       1       BP GO:0030330                                               DNA damage response, signal transduction by p53 class mediator      1/10  41/8068 0.049698637 0.1516601 0.1074590     1                 ZNF385A
# 135       1       BP GO:2000573                                                              positive regulation of DNA biosynthetic process      1/10  41/8068 0.049698637 0.1516601 0.1074590     1                  CYP1B1
# 136       1       BP GO:0010942                                                                            positive regulation of cell death      2/10 299/8068 0.050607577 0.1516601 0.1074590     2         TNFRSF1A/CYP1B1
# 137       1       BP GO:0061515                                                                                     myeloid cell development      1/10  42/8068 0.050882518 0.1516601 0.1074590     1                 ZNF385A
# 138       1       BP GO:0098586                                                                                   cellular response to virus      1/10  42/8068 0.050882518 0.1516601 0.1074590     1                   RIOK3
# 139       1       BP GO:0071407                                                                 cellular response to organic cyclic compound      2/10 303/8068 0.051835430 0.1516601 0.1074590     2            RIOK3/CYP1B1
# 140       1       BP GO:0014909                                                                                 smooth muscle cell migration      1/10  43/8068 0.052065072 0.1516601 0.1074590     1                  CYP1B1
# 141       1       BP GO:1902532                                                     negative regulation of intracellular signal transduction      2/10 306/8068 0.052763692 0.1516601 0.1074590     2           ZNF385A/RIOK3
# 142       1       BP GO:0010035                                                                              response to inorganic substance      2/10 307/8068 0.053074507 0.1516601 0.1074590     2             MT1F/CYP1B1
# 143       1       BP GO:0001676                                                                      long-chain fatty acid metabolic process      1/10  44/8068 0.053246299 0.1516601 0.1074590     1                  CYP1B1
# 144       1       BP GO:0032088                                               negative regulation of NF-kappaB transcription factor activity      1/10  44/8068 0.053246299 0.1516601 0.1074590     1                  CYP1B1
# 145       1       BP GO:0051591                                                                                             response to cAMP      1/10  44/8068 0.053246299 0.1516601 0.1074590     1                  CYP1B1
# 146       1       BP GO:0006720                                                                                 isoprenoid metabolic process      1/10  45/8068 0.054426202 0.1520846 0.1077597     1                  CYP1B1
# 147       1       BP GO:0033209                                                             tumor necrosis factor-mediated signaling pathway      1/10  45/8068 0.054426202 0.1520846 0.1077597     1                TNFRSF1A
# 148       1       BP GO:0032479                                                                   regulation of type I interferon production      1/10  46/8068 0.055604781 0.1520846 0.1077597     1                   RIOK3
# 149       1       BP GO:0032606                                                                                 type I interferon production      1/10  46/8068 0.055604781 0.1520846 0.1077597     1                   RIOK3
# 150       1       BP GO:0050729                                                                 positive regulation of inflammatory response      1/10  46/8068 0.055604781 0.1520846 0.1077597     1                TNFRSF1A
# 151       1       BP GO:0120254                                                                          olefinic compound metabolic process      1/10  46/8068 0.055604781 0.1520846 0.1077597     1                  CYP1B1
# 152       1       BP GO:0008625                                             extrinsic apoptotic signaling pathway via death domain receptors      1/10  47/8068 0.056782037 0.1522791 0.1078975     1                TNFRSF1A
# 153       1       BP GO:0062207                                                 regulation of pattern recognition receptor signaling pathway      1/10  47/8068 0.056782037 0.1522791 0.1078975     1                   RIOK3
# 154       1       BP GO:0097306                                                                                 cellular response to alcohol      1/10  47/8068 0.056782037 0.1522791 0.1078975     1                  CYP1B1
# 155       1       BP GO:2001021                                                       negative regulation of response to DNA damage stimulus      1/10  48/8068 0.057957973 0.1544300 0.1094215     1                 ZNF385A
# 156       1       BP GO:0002832                                                           negative regulation of response to biotic stimulus      1/10  49/8068 0.059132589 0.1555526 0.1102170     1                   RIOK3
# 157       1       BP GO:0003300                                                                                   cardiac muscle hypertrophy      1/10  49/8068 0.059132589 0.1555526 0.1102170     1                TNFRSF1A
# 158       1       BP GO:0014812                                                                                        muscle cell migration      1/10  50/8068 0.060305886 0.1555885 0.1102424     1                  CYP1B1
# 159       1       BP GO:0014897                                                                                  striated muscle hypertrophy      1/10  50/8068 0.060305886 0.1555885 0.1102424     1                TNFRSF1A
# 160       1       BP GO:0043502                                                                              regulation of muscle adaptation      1/10  50/8068 0.060305886 0.1555885 0.1102424     1                TNFRSF1A
# 161       1       BP GO:0014896                                                                                           muscle hypertrophy      1/10  51/8068 0.061477867 0.1555885 0.1102424     1                TNFRSF1A
# 162       1       BP GO:0045089                                                                positive regulation of innate immune response      1/10  51/8068 0.061477867 0.1555885 0.1102424     1                   RIOK3
# 163       1       BP GO:0072332                                                  intrinsic apoptotic signaling pathway by p53 class mediator      1/10  51/8068 0.061477867 0.1555885 0.1102424     1                 ZNF385A
# 164       1       BP GO:0034248                                                               regulation of cellular amide metabolic process      2/10 335/8068 0.062052696 0.1555885 0.1102424     2          ZNF385A/CYP1B1
# 165       1       BP GO:0050728                                                                 negative regulation of inflammatory response      1/10  52/8068 0.062648532 0.1555885 0.1102424     1                TNFRSF1A
# 166       1       BP GO:2000112                                                    regulation of cellular macromolecule biosynthetic process      2/10 340/8068 0.063710247 0.1555885 0.1102424     2          ZNF385A/CYP1B1
# 167       1       BP GO:0035270                                                                                 endocrine system development      1/10  53/8068 0.063817883 0.1555885 0.1102424     1                  CYP1B1
# 168       1       BP GO:0045446                                                                             endothelial cell differentiation      1/10  53/8068 0.063817883 0.1555885 0.1102424     1                TNFRSF1A
# 169       1       BP GO:0046683                                                                                 response to organophosphorus      1/10  53/8068 0.063817883 0.1555885 0.1102424     1                  CYP1B1
# 170       1       BP GO:0071345                                                                       cellular response to cytokine stimulus      2/10 341/8068 0.064043672 0.1555885 0.1102424     2         TNFRSF1A/CYP1B1
# 171       1       BP GO:0030856                                                                regulation of epithelial cell differentiation      1/10  55/8068 0.066152646 0.1597722 0.1132067     1                TNFRSF1A
# 172       1       BP GO:0042274                                                                           ribosomal small subunit biogenesis      1/10  57/8068 0.068482168 0.1632477 0.1156693     1                   RIOK3
# 173       1       BP GO:0043500                                                                                            muscle adaptation      1/10  57/8068 0.068482168 0.1632477 0.1156693     1                TNFRSF1A
# 174       1       BP GO:0010608                                                           post-transcriptional regulation of gene expression      2/10 355/8068 0.068777469 0.1632477 0.1156693     2          ZNF385A/CYP1B1
# 175       1       BP GO:0070301                                                                       cellular response to hydrogen peroxide      1/10  58/8068 0.069644966 0.1634282 0.1157973     1                  CYP1B1
# 176       1       BP GO:1901655                                                                                  cellular response to ketone      1/10  58/8068 0.069644966 0.1634282 0.1157973     1                  CYP1B1
# 177       1       BP GO:0006304                                                                                             DNA modification      1/10  59/8068 0.070806458 0.1642869 0.1164056     1                  CYP1B1
# 178       1       BP GO:0042742                                                                                defense response to bacterium      1/10  59/8068 0.070806458 0.1642869 0.1164056     1                TNFRSF1A
# 179       1       BP GO:0010507                                                                             negative regulation of autophagy      1/10  60/8068 0.071966645 0.1657330 0.1174303     1                    CTSA
# 180       1       BP GO:0098542                                                                           defense response to other organism      2/10 365/8068 0.072232303 0.1657330 0.1174303     2          TNFRSF1A/RIOK3
# 181       1       BP GO:0003158                                                                                      endothelium development      1/10  61/8068 0.073125528 0.1659387 0.1175760     1                TNFRSF1A
# 182       1       BP GO:0043534                                                                      blood vessel endothelial cell migration      1/10  61/8068 0.073125528 0.1659387 0.1175760     1                  CYP1B1
# 183       1       BP GO:0014074                                                                       response to purine-containing compound      1/10  62/8068 0.074283108 0.1667333 0.1181390     1                  CYP1B1
# 184       1       BP GO:0032355                                                                                        response to estradiol      1/10  62/8068 0.074283108 0.1667333 0.1181390     1                  CYP1B1
# 185       1       BP GO:0036473                                                                   cell death in response to oxidative stress      1/10  63/8068 0.075439387 0.1676807 0.1188103     1                  CYP1B1
# 186       1       BP GO:1901699                                                                       cellular response to nitrogen compound      2/10 377/8068 0.076456362 0.1676807 0.1188103     2            RIOK3/CYP1B1
# 187       1       BP GO:0071466                                                                     cellular response to xenobiotic stimulus      1/10  64/8068 0.076594366 0.1676807 0.1188103     1                  CYP1B1
# 188       1       BP GO:0002833                                                           positive regulation of response to biotic stimulus      1/10  65/8068 0.077748046 0.1676807 0.1188103     1                   RIOK3
# 189       1       BP GO:0046916                                                                    cellular transition metal ion homeostasis      1/10  65/8068 0.077748046 0.1676807 0.1188103     1                    MT1F
# 190       1       BP GO:0007584                                                                                         response to nutrient      1/10  66/8068 0.078900429 0.1676807 0.1188103     1                  CYP1B1
# 191       1       BP GO:0051384                                                                                   response to glucocorticoid      1/10  66/8068 0.078900429 0.1676807 0.1188103     1                  CYP1B1
# 192       1       BP GO:0098754                                                                                               detoxification      1/10  66/8068 0.078900429 0.1676807 0.1188103     1                    MT1F
# 193       1       BP GO:0034097                                                                                         response to cytokine      2/10 385/8068 0.079318515 0.1676807 0.1188103     2         TNFRSF1A/CYP1B1
# 194       1       BP GO:0001959                                                            regulation of cytokine-mediated signaling pathway      1/10  67/8068 0.080051515 0.1676807 0.1188103     1                TNFRSF1A
# 195       1       BP GO:0045766                                                                          positive regulation of angiogenesis      1/10  67/8068 0.080051515 0.1676807 0.1188103     1                  CYP1B1
# 196       1       BP GO:1904018                                                               positive regulation of vasculature development      1/10  67/8068 0.080051515 0.1676807 0.1188103     1                  CYP1B1
# 197       1       BP GO:0032101                                                                  regulation of response to external stimulus      2/10 390/8068 0.081125655 0.1676807 0.1188103     2          TNFRSF1A/RIOK3
# 198       1       BP GO:0008584                                                                                       male gonad development      1/10  68/8068 0.081201307 0.1676807 0.1188103     1                  CYP1B1
# 199       1       BP GO:0046546                                                           development of primary male sexual characteristics      1/10  68/8068 0.081201307 0.1676807 0.1188103     1                  CYP1B1
# 200       1       BP GO:2000278                                                                       regulation of DNA biosynthetic process      1/10  68/8068 0.081201307 0.1676807 0.1188103     1                  CYP1B1
# 201       1       BP GO:2001243                                                 negative regulation of intrinsic apoptotic signaling pathway      1/10  70/8068 0.083497012 0.1711607 0.1212761     1                 ZNF385A
# 202       1       BP GO:0045598                                                                       regulation of fat cell differentiation      1/10  71/8068 0.084642927 0.1711607 0.1212761     1                 ZNF385A
# 203       1       BP GO:0050731                                                     positive regulation of peptidyl-tyrosine phosphorylation      1/10  71/8068 0.084642927 0.1711607 0.1212761     1                TNFRSF1A
# 204       1       BP GO:0060041                                                                        retina development in camera-type eye      1/10  71/8068 0.084642927 0.1711607 0.1212761     1                  CYP1B1
# 205       1       BP GO:0060759                                                                  regulation of response to cytokine stimulus      1/10  72/8068 0.085787553 0.1711607 0.1212761     1                TNFRSF1A
# 206       1       BP GO:1901796                                                      regulation of signal transduction by p53 class mediator      1/10  72/8068 0.085787553 0.1711607 0.1212761     1                 ZNF385A
# 207       1       BP GO:1903531                                                                     negative regulation of secretion by cell      1/10  72/8068 0.085787553 0.1711607 0.1212761     1                TNFRSF1A
# 208       1       BP GO:0002221                                                               pattern recognition receptor signaling pathway      1/10  74/8068 0.088072940 0.1740389 0.1233154     1                   RIOK3
# 209       1       BP GO:0042445                                                                                    hormone metabolic process      1/10  74/8068 0.088072940 0.1740389 0.1233154     1                  CYP1B1
# 210       1       BP GO:0031960                                                                                   response to corticosteroid      1/10  75/8068 0.089213704 0.1746221 0.1237287     1                  CYP1B1
# 211       1       BP GO:0046661                                                                                     male sex differentiation      1/10  75/8068 0.089213704 0.1746221 0.1237287     1                  CYP1B1
# 212       1       BP GO:0016125                                                                                     sterol metabolic process      1/10  77/8068 0.091491381 0.1765698 0.1251087     1                  CYP1B1
# 213       1       BP GO:0051048                                                                             negative regulation of secretion      1/10  77/8068 0.091491381 0.1765698 0.1251087     1                TNFRSF1A
# 214       1       BP GO:2000377                                                      regulation of reactive oxygen species metabolic process      1/10  77/8068 0.091491381 0.1765698 0.1251087     1                  CYP1B1
# 215       1       BP GO:0042177                                                             negative regulation of protein catabolic process      1/10  82/8068 0.097163162 0.1866437 0.1322466     1                    CTSA
# 216       1       BP GO:0042542                                                                                response to hydrogen peroxide      1/10  84/8068 0.099422936 0.1892243 0.1340751     1                  CYP1B1
# 217       1       BP GO:0043433                                             negative regulation of DNA-binding transcription factor activity      1/10  84/8068 0.099422936 0.1892243 0.1340751     1                  CYP1B1
# 218       1       BP GO:0055076                                                                             transition metal ion homeostasis      1/10  85/8068 0.100550914 0.1904932 0.1349742     1                    MT1F
# 219       1       BP GO:0051129                                                       negative regulation of cellular component organization      2/10 444/8068 0.101483400 0.1913819 0.1356039     2          TNFRSF1A/RIOK3
# 220       1       BP GO:0014070                                                                          response to organic cyclic compound      2/10 450/8068 0.103834379 0.1936862 0.1372365     2            RIOK3/CYP1B1
# 221       1       BP GO:0007599                                                                                                   hemostasis      1/10  88/8068 0.103927220 0.1936862 0.1372365     1                 ZNF385A
# 222       1       BP GO:0009451                                                                                             RNA modification      1/10  89/8068 0.105050119 0.1936862 0.1372365     1                   PUS7L
# 223       1       BP GO:0034614                                                                 cellular response to reactive oxygen species      1/10  89/8068 0.105050119 0.1936862 0.1372365     1                  CYP1B1
# 224       1       BP GO:0043123                                                   positive regulation of I-kappaB kinase/NF-kappaB signaling      1/10  89/8068 0.105050119 0.1936862 0.1372365     1                TNFRSF1A
# 225       1       BP GO:0019752                                                                            carboxylic acid metabolic process      2/10 457/8068 0.106598170 0.1956669 0.1386400     2         TNFRSF1A/CYP1B1
# 226       1       BP GO:0045088                                                                         regulation of innate immune response      1/10  92/8068 0.108411219 0.1981143 0.1403741     1                   RIOK3
# 227       1       BP GO:0043436                                                                                    oxoacid metabolic process      2/10 465/8068 0.109783851 0.1995878 0.1414181     2         TNFRSF1A/CYP1B1
# 228       1       BP GO:0006082                                                                               organic acid metabolic process      2/10 466/8068 0.110184056 0.1995878 0.1414181     2         TNFRSF1A/CYP1B1
# 324       1       MF GO:0004180                                                                                    carboxypeptidase activity      1/10  15/8223 0.018102332 0.1625528 0.1326089     1                    CTSA
# 325       1       MF GO:0016866                                                                          intramolecular transferase activity      1/10  15/8223 0.018102332 0.1625528 0.1326089     1                   PUS7L
# 326       1       MF GO:0004497                                                                                       monooxygenase activity      1/10  22/8223 0.026448712 0.1625528 0.1326089     1                  CYP1B1
# 327       1       MF GO:0008138                                                       protein tyrosine/serine/threonine phosphatase activity      1/10  23/8223 0.027635825 0.1625528 0.1326089     1                  DUSP14
# 328       1       MF GO:0016836                                                                                         hydro-lyase activity      1/10  31/8223 0.037085929 0.1625528 0.1326089     1                  CYP1B1
# 329       1       MF GO:0019955                                                                                             cytokine binding      1/10  36/8223 0.042950187 0.1625528 0.1326089     1                TNFRSF1A
# 330       1       MF GO:0008238                                                                                        exopeptidase activity      1/10  41/8223 0.048782284 0.1625528 0.1326089     1                    CTSA
# 331       1       MF GO:0016835                                                                                 carbon-oxygen lyase activity      1/10  41/8223 0.048782284 0.1625528 0.1326089     1                  CYP1B1
# 332       1       MF GO:0020037                                                                                                 heme binding      1/10  42/8223 0.049944857 0.1625528 0.1326089     1                  CYP1B1
# 333       1       MF GO:0008236                                                                               serine-type peptidase activity      1/10  44/8223 0.052266169 0.1625528 0.1326089     1                    CTSA
# 334       1       MF GO:0017171                                                                                    serine hydrolase activity      1/10  47/8223 0.055738566 0.1625528 0.1326089     1                    CTSA
# 335       1       MF GO:0002039                                                                                                  p53 binding      1/10  48/8223 0.056893484 0.1625528 0.1326089     1                 ZNF385A
# 336       1       MF GO:0017018                                                                                  myosin phosphatase activity      1/10  48/8223 0.056893484 0.1625528 0.1326089     1                  DUSP14
# 337       1       MF GO:0046906                                                                                         tetrapyrrole binding      1/10  48/8223 0.056893484 0.1625528 0.1326089     1                  CYP1B1
# 338       1       MF GO:0005506                                                                                             iron ion binding      1/10  54/8223 0.063796339 0.1625837 0.1326341     1                  CYP1B1
# 339       1       MF GO:0004725                                                                        protein tyrosine phosphatase activity      1/10  56/8223 0.066087165 0.1625837 0.1326341     1                  DUSP14
# 340       1       MF GO:0016705        oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen      1/10  64/8223 0.075200105 0.1625837 0.1326341     1                  CYP1B1
# 341       1       MF GO:0004722                                                                protein serine/threonine phosphatase activity      1/10  65/8223 0.076333578 0.1625837 0.1326341     1                  DUSP14
# 342       1       MF GO:0008270                                                                                             zinc ion binding      2/10 402/8223 0.082737528 0.1625837 0.1326341     2            ZNF385A/MT1F
# 343       1       MF GO:0002020                                                                                             protease binding      1/10  71/8223 0.083108193 0.1625837 0.1326341     1                   RIOK3
# 344       1       MF GO:0003730                                                                                          mRNA 3'-UTR binding      1/10  73/8223 0.085356441 0.1625837 0.1326341     1                 ZNF385A
# 345       1       MF GO:0016853                                                                                           isomerase activity      1/10  94/8223 0.108665269 0.1975732 0.1611782     1                   PUS7L
# 346       1       MF GO:0016829                                                                                               lyase activity      1/10  99/8223 0.114135578 0.1984967 0.1619315     1                  CYP1B1

sort(table(subset(as.data.frame(go), p.adjust < 0.2)$geneID), decreasing = TRUE)
#          CYP1B1                TNFRSF1A                 ZNF385A                   RIOK3
#              81                      39                      31                      24
# TNFRSF1A/CYP1B1                    MT1F          TNFRSF1A/RIOK3                    CTSA
#              19                      13                      10                       7
#    RIOK3/CYP1B1                   PUS7L                  DUSP14          ZNF385A/CYP1B1
#               6                       5                       4                       4
#     MT1F/CYP1B1 ZNF385A/TNFRSF1A/CYP1B1  SUSD6/ZNF385A/TNFRSF1A            ZNF385A/MT1F
#               2                       2                       1                       1
#   ZNF385A/RIOK3        ZNF385A/TNFRSF1A
#               1                       1

## Looks like 2 genes dominate the results, so I'd rather we focus on those that GO results
## https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP1B1&keywords=CYP1B1
## https://pubmed.ncbi.nlm.nih.gov/?term=CYP1B1+alzheimer%27s+disease
##
## https://www.genecards.org/cgi-bin/carddisp.pl?gene=TNFRSF1A&keywords=TNFRSF1A
## https://pubmed.ncbi.nlm.nih.gov/?term=TNFRSF1A+alzheimer%27s+disease
##
## https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZNF385A&keywords=ZNF385A
## https://pubmed.ncbi.nlm.nih.gov/?term=ZNF385A+alzheimer%27s+disease ## leads to only 1 paper from 2017
##
## https://www.genecards.org/cgi-bin/carddisp.pl?gene=RIOK3&keywords=RIOK3
## https://pubmed.ncbi.nlm.nih.gov/?term=RIOK3+alzheimer%27s+disease ## leads to only 1 paper from 2020

subset(as.data.frame(kegg), p.adjust < 0.3)
#    Cluster       ID                                       Description GeneRatio  BgRatio       pvalue  p.adjust    qvalue                                   geneID Count
# 1       -1 hsa03018                                   RNA degradation      7/99  61/3826 0.0008717714 0.1935332 0.1935332 9652/23016/5214/11044/56915/246175/84186     7
# 2        1 hsa00140                      Steroid hormone biosynthesis       1/4  10/3826 0.0104179353 0.1521979 0.1237973                                     1545     1
# 3        1 hsa04913                           Ovarian steroidogenesis       1/4  13/3826 0.0135273819 0.1521979 0.1237973                                     1545     1
# 4        1 hsa04215                      Apoptosis - multiple species       1/4  15/3826 0.0155962714 0.1521979 0.1237973                                     7132     1
# 5        1 hsa05204             Chemical carcinogenesis - DNA adducts       1/4  17/3826 0.0176619049 0.1521979 0.1237973                                     1545     1
# 6        1 hsa00380                             Tryptophan metabolism       1/4  18/3826 0.0186935017 0.1521979 0.1237973                                     1545     1
# 7        1 hsa00980      Metabolism of xenobiotics by cytochrome P450       1/4  20/3826 0.0207542577 0.1521979 0.1237973                                     1545     1
# 8        1 hsa04060            Cytokine-cytokine receptor interaction       1/4  28/3826 0.0289648417 0.1558377 0.1267579                                     7132     1
# 9        1 hsa04064                      NF-kappa B signaling pathway       1/4  31/3826 0.0320304596 0.1558377 0.1267579                                     7132     1
# 10       1 hsa04978                                Mineral absorption       1/4  31/3826 0.0320304596 0.1558377 0.1267579                                     4494     1
# 11       1 hsa04920                   Adipocytokine signaling pathway       1/4  37/3826 0.0381399135 0.1558377 0.1267579                                     7132     1
# 12       1 hsa04668                             TNF signaling pathway       1/4  52/3826 0.0532870220 0.1558377 0.1267579                                     7132     1
# 13       1 hsa05142                                    Chagas disease       1/4  53/3826 0.0542904273 0.1558377 0.1267579                                     7132     1
# 14       1 hsa05145                                     Toxoplasmosis       1/4  53/3826 0.0542904273 0.1558377 0.1267579                                     7132     1
# 15       1 hsa04380                        Osteoclast differentiation       1/4  60/3826 0.0612919494 0.1558377 0.1267579                                     7132     1
# 16       1 hsa04931                                Insulin resistance       1/4  65/3826 0.0662691820 0.1558377 0.1267579                                     7132     1
# 17       1 hsa04936                           Alcoholic liver disease       1/4  67/3826 0.0682545228 0.1558377 0.1267579                                     7132     1
# 18       1 hsa04210                                         Apoptosis       1/4  71/3826 0.0722157046 0.1558377 0.1267579                                     7132     1
# 19       1 hsa05164                                       Influenza A       1/4  72/3826 0.0732040234 0.1558377 0.1267579                                     7132     1
# 20       1 hsa05160                                       Hepatitis C       1/4  75/3826 0.0761642427 0.1558377 0.1267579                                     7132     1
# 21       1 hsa05152                                      Tuberculosis       1/4  77/3826 0.0781337788 0.1558377 0.1267579                                     7132     1
# 22       1 hsa04217                                       Necroptosis       1/4  78/3826 0.0791173650 0.1558377 0.1267579                                     7132     1
# 23       1 hsa04071                    Sphingolipid signaling pathway       1/4  82/3826 0.0830438414 0.1558377 0.1267579                                     7132     1
# 24       1 hsa05418            Fluid shear stress and atherosclerosis       1/4  83/3826 0.0840234954 0.1558377 0.1267579                                     7132     1
# 25       1 hsa04142                                          Lysosome       1/4  84/3826 0.0850023642 0.1558377 0.1267579                                     5476     1
# 26       1 hsa05417                         Lipid and atherosclerosis       1/4  95/3826 0.0957182188 0.1559649 0.1268614                                     7132     1
# 27       1 hsa05207     Chemical carcinogenesis - receptor activation       1/4  96/3826 0.0966876981 0.1559649 0.1268614                                     1545     1
# 28       1 hsa05206                               MicroRNAs in cancer       1/4  98/3826 0.0986243178 0.1559649 0.1268614                                     1545     1
# 29       1 hsa04150                            mTOR signaling pathway       1/4 102/3826 0.1024882137 0.1559649 0.1268614                                     7132     1
# 30       1 hsa05167   Kaposi sarcoma-associated herpesvirus infection       1/4 103/3826 0.1034522435 0.1559649 0.1268614                                     7132     1
# 31       1 hsa05130             Pathogenic Escherichia coli infection       1/4 106/3826 0.1063396738 0.1559649 0.1268614                                     7132     1
# 32       1 hsa04932                 Non-alcoholic fatty liver disease       1/4 115/3826 0.1149601292 0.1570929 0.1277789                                     7132     1
# 33       1 hsa05163                   Human cytomegalovirus infection       1/4 118/3826 0.1178197067 0.1570929 0.1277789                                     7132     1
# 34       1 hsa05166           Human T-cell leukemia virus 1 infection       1/4 118/3826 0.1178197067 0.1570929 0.1277789                                     7132     1
# 35       1 hsa05170          Human immunodeficiency virus 1 infection       1/4 122/3826 0.1216216907 0.1573928 0.1280228                                     7132     1
# 36       1 hsa05168                  Herpes simplex virus 1 infection       1/4 130/3826 0.1291887777 0.1624087 0.1321028                                     7132     1
# 37       1 hsa05171                    Coronavirus disease - COVID-19       1/4 139/3826 0.1376431938 0.1682306 0.1368383                                     7132     1
# 38       1 hsa05131                                       Shigellosis       1/4 164/3826 0.1608048274 0.1784613 0.1451599                                     7132     1
# 39       1 hsa05165                    Human papillomavirus infection       1/4 164/3826 0.1608048274 0.1784613 0.1451599                                     7132     1
# 40       1 hsa05132                              Salmonella infection       1/4 168/3826 0.1644669314 0.1784613 0.1451599                                     7132     1
# 41       1 hsa04010                            MAPK signaling pathway       1/4 169/3826 0.1653805815 0.1784613 0.1451599                                     7132     1
# 42       1 hsa05208 Chemical carcinogenesis - reactive oxygen species       1/4 170/3826 0.1662934822 0.1784613 0.1451599                                     1545     1
# 43       1 hsa05010                                 Alzheimer disease       1/4 270/3826 0.2538684550 0.2597724 0.2112981                                     7132     1
# 44       1 hsa05014                     Amyotrophic lateral sclerosis       1/4 270/3826 0.2538684550 0.2597724 0.2112981                                     7132     1

## Note that gene with ENTREZ ID 7132 is this one https://www.ncbi.nlm.nih.gov/gene/7132
## "TNFRSF1A TNF receptor superfamily member 1A [ Homo sapiens (human) ]"

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.2.2 (2022-10-31)
#  os       macOS Ventura 13.0.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2023-04-05
#  rstudio  2023.03.0+386 Cherry Blossom (desktop)
#  pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  AnnotationDbi          * 1.60.2    2023-03-12 [1] Bioconductor
#  AnnotationHub            3.6.0     2022-11-01 [1] Bioconductor
#  ape                      5.7-1     2023-03-13 [1] CRAN (R 4.2.0)
#  aplot                    0.1.10    2023-03-08 [1] CRAN (R 4.2.0)
#  attempt                  0.3.1     2020-05-03 [1] CRAN (R 4.2.0)
#  beachmat                 2.14.0    2022-11-01 [1] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [1] CRAN (R 4.2.0)
#  benchmarkme              1.0.8     2022-06-12 [1] CRAN (R 4.2.0)
#  benchmarkmeData          1.0.4     2020-04-23 [1] CRAN (R 4.2.0)
#  Biobase                * 2.58.0    2022-11-01 [1] Bioconductor
#  BiocFileCache            2.6.1     2023-02-19 [1] Bioconductor
#  BiocGenerics           * 0.44.0    2022-11-01 [1] Bioconductor
#  BiocIO                   1.8.0     2022-11-01 [1] Bioconductor
#  BiocManager              1.30.20   2023-02-24 [1] CRAN (R 4.2.0)
#  BiocNeighbors            1.16.0    2022-11-01 [1] Bioconductor
#  BiocParallel             1.32.5    2022-12-25 [1] Bioconductor
#  BiocSingular             1.14.0    2022-11-01 [1] Bioconductor
#  biocthis                 1.9.3     2023-03-10 [1] Github (lcolladotor/biocthis@910b96b)
#  BiocVersion              3.16.0    2022-09-20 [1] Bioconductor
#  Biostrings               2.66.0    2022-11-01 [1] Bioconductor
#  bit                      4.0.5     2022-11-15 [1] CRAN (R 4.2.2)
#  bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.2.0)
#  bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.2.0)
#  blob                     1.2.4     2023-03-17 [1] CRAN (R 4.2.0)
#  brio                     1.1.3     2021-11-30 [1] CRAN (R 4.2.0)
#  bslib                    0.4.2     2022-12-16 [1] CRAN (R 4.2.2)
#  cachem                   1.0.7     2023-02-24 [1] CRAN (R 4.2.0)
#  callr                    3.7.3     2022-11-02 [1] CRAN (R 4.2.2)
#  cli                      3.6.0     2023-01-09 [1] CRAN (R 4.2.0)
#  clusterProfiler        * 4.6.2     2023-03-05 [1] Bioconductor
#  codetools                0.2-19    2023-02-01 [1] CRAN (R 4.2.0)
#  colorout                 1.2-2     2022-03-01 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.2.0)
#  config                   0.3.1     2020-12-17 [1] CRAN (R 4.2.0)
#  cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.2.0)
#  crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.2.0)
#  curl                     5.0.0     2023-01-12 [1] CRAN (R 4.2.0)
#  data.table               1.14.8    2023-02-17 [1] CRAN (R 4.2.2)
#  DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.2.0)
#  dbplyr                   2.3.1     2023-02-24 [1] CRAN (R 4.2.0)
#  DelayedArray             0.24.0    2022-11-01 [1] Bioconductor
#  DelayedMatrixStats       1.20.0    2022-11-01 [1] Bioconductor
#  devtools               * 2.4.5     2022-10-11 [1] CRAN (R 4.2.0)
#  digest                   0.6.31    2022-12-11 [1] CRAN (R 4.2.0)
#  doParallel               1.0.17    2022-02-07 [1] CRAN (R 4.2.0)
#  DOSE                     3.24.2    2022-11-23 [1] Bioconductor
#  dotCall64                1.0-2     2022-10-03 [1] CRAN (R 4.2.1)
#  downloader               0.4       2015-07-09 [1] CRAN (R 4.2.0)
#  dplyr                    1.1.0     2023-01-29 [1] CRAN (R 4.2.0)
#  dqrng                    0.3.0     2021-05-01 [1] CRAN (R 4.2.0)
#  DropletUtils             1.18.1    2022-11-23 [1] Bioconductor
#  DT                       0.27      2023-01-17 [1] CRAN (R 4.2.0)
#  edgeR                    3.40.2    2023-01-22 [1] Bioconductor
#  ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.2.0)
#  enrichplot               1.18.3    2022-12-07 [1] Bioconductor
#  evaluate                 0.20      2023-01-17 [1] CRAN (R 4.2.0)
#  ExperimentHub            2.6.0     2022-11-01 [1] Bioconductor
#  fansi                    1.0.4     2023-01-22 [1] CRAN (R 4.2.0)
#  farver                   2.1.1     2022-07-06 [1] CRAN (R 4.2.1)
#  fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.2.0)
#  fastmatch                1.1-3     2021-07-23 [1] CRAN (R 4.2.0)
#  fgsea                    1.24.0    2022-11-01 [1] Bioconductor
#  fields                   14.1      2022-08-12 [1] CRAN (R 4.2.0)
#  filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.2.0)
#  foreach                  1.5.2     2022-02-02 [1] CRAN (R 4.2.0)
#  fs                       1.6.1     2023-02-06 [1] CRAN (R 4.2.0)
#  generics                 0.1.3     2022-07-05 [1] CRAN (R 4.2.0)
#  GenomeInfoDb           * 1.34.9    2023-02-02 [1] Bioconductor
#  GenomeInfoDbData         1.2.9     2022-11-02 [1] Bioconductor
#  GenomicAlignments        1.34.1    2023-03-09 [1] Bioconductor
#  GenomicRanges          * 1.50.2    2022-12-18 [1] Bioconductor
#  ggbeeswarm               0.7.1     2022-12-16 [1] CRAN (R 4.2.2)
#  ggforce                  0.4.1     2022-10-04 [1] CRAN (R 4.2.1)
#  ggfun                    0.0.9     2022-11-21 [1] CRAN (R 4.2.0)
#  ggplot2                * 3.4.1     2023-02-10 [1] CRAN (R 4.2.0)
#  ggplotify                0.1.0     2021-09-02 [1] CRAN (R 4.2.0)
#  ggraph                   2.1.0     2022-10-09 [1] CRAN (R 4.2.0)
#  ggrepel                  0.9.3     2023-02-03 [1] CRAN (R 4.2.0)
#  ggtree                   3.6.2     2022-11-10 [1] Bioconductor
#  glue                     1.6.2     2022-02-24 [1] CRAN (R 4.2.0)
#  GO.db                    3.16.0    2022-11-02 [1] Bioconductor
#  golem                    0.4.0     2023-03-12 [1] CRAN (R 4.2.0)
#  GOSemSim                 2.24.0    2022-11-01 [1] Bioconductor
#  graphlayouts             0.8.4     2022-11-24 [1] CRAN (R 4.2.2)
#  gridExtra                2.3       2017-09-09 [1] CRAN (R 4.2.0)
#  gridGraphics             0.5-1     2020-12-13 [1] CRAN (R 4.2.0)
#  gson                     0.1.0     2023-03-07 [1] CRAN (R 4.2.0)
#  gtable                   0.3.2     2023-03-17 [1] CRAN (R 4.2.2)
#  HDF5Array                1.26.0    2022-11-01 [1] Bioconductor
#  HDO.db                   0.99.1    2022-11-02 [1] Bioconductor
#  here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.2.0)
#  hms                      1.1.2     2022-08-19 [1] CRAN (R 4.2.0)
#  htmltools                0.5.4     2022-12-07 [1] CRAN (R 4.2.0)
#  htmlwidgets              1.6.2     2023-03-17 [1] CRAN (R 4.2.2)
#  httpuv                   1.6.9     2023-02-14 [1] CRAN (R 4.2.0)
#  httr                     1.4.5     2023-02-24 [1] CRAN (R 4.2.0)
#  igraph                   1.4.1     2023-02-24 [1] CRAN (R 4.2.0)
#  interactiveDisplayBase   1.36.0    2022-11-01 [1] Bioconductor
#  IRanges                * 2.32.0    2022-11-01 [1] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [1] CRAN (R 4.2.1)
#  iterators                1.0.14    2022-02-05 [1] CRAN (R 4.2.0)
#  jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.2.0)
#  jsonlite                 1.8.4     2022-12-06 [1] CRAN (R 4.2.0)
#  KEGGREST                 1.38.0    2022-11-01 [1] Bioconductor
#  knitr                    1.42      2023-01-25 [1] CRAN (R 4.2.0)
#  labeling                 0.4.2     2020-10-20 [1] CRAN (R 4.2.0)
#  later                    1.3.0     2021-08-18 [1] CRAN (R 4.2.0)
#  lattice                  0.20-45   2021-09-22 [1] CRAN (R 4.2.2)
#  lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.2.0)
#  lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.2.1)
#  limma                    3.54.2    2023-02-28 [1] Bioconductor
#  locfit                   1.5-9.7   2023-01-02 [1] CRAN (R 4.2.0)
#  lubridate                1.9.2     2023-02-10 [1] CRAN (R 4.2.0)
#  magick                   2.7.4     2023-03-09 [1] CRAN (R 4.2.0)
#  magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.2.0)
#  maps                     3.4.1     2022-10-30 [1] CRAN (R 4.2.0)
#  MASS                     7.3-58.3  2023-03-07 [1] CRAN (R 4.2.0)
#  Matrix                   1.5-3     2022-11-11 [1] CRAN (R 4.2.0)
#  MatrixGenerics         * 1.10.0    2022-11-01 [1] Bioconductor
#  matrixStats            * 0.63.0    2022-11-18 [1] CRAN (R 4.2.0)
#  memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.2.0)
#  mime                     0.12      2021-09-28 [1] CRAN (R 4.2.0)
#  miniUI                   0.1.1.1   2018-05-18 [1] CRAN (R 4.2.0)
#  munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.2.0)
#  nlme                     3.1-162   2023-01-31 [1] CRAN (R 4.2.0)
#  org.Hs.eg.db           * 3.16.0    2022-11-02 [1] Bioconductor
#  paletteer                1.5.0     2022-10-19 [1] CRAN (R 4.2.0)
#  patchwork                1.1.2     2022-08-19 [1] CRAN (R 4.2.0)
#  pillar                   1.8.1     2022-08-19 [1] CRAN (R 4.2.0)
#  pkgbuild                 1.4.0     2022-11-27 [1] CRAN (R 4.2.2)
#  pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.2.0)
#  pkgload                  1.3.2     2022-11-16 [1] CRAN (R 4.2.2)
#  plotly                   4.10.1    2022-11-07 [1] CRAN (R 4.2.0)
#  plyr                     1.8.8     2022-11-11 [1] CRAN (R 4.2.0)
#  png                      0.1-8     2022-11-29 [1] CRAN (R 4.2.0)
#  polyclip                 1.10-4    2022-10-20 [1] CRAN (R 4.2.0)
#  prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.2.0)
#  processx                 3.8.0     2022-10-26 [1] CRAN (R 4.2.0)
#  profvis                  0.3.7     2020-11-02 [1] CRAN (R 4.2.0)
#  promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.2.0)
#  prompt                   1.0.1     2022-03-01 [1] Github (gaborcsardi/prompt@7ef0f2e)
#  ps                       1.7.2     2022-10-26 [1] CRAN (R 4.2.0)
#  purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.2.0)
#  qvalue                   2.30.0    2022-11-01 [1] Bioconductor
#  R.cache                  0.16.0    2022-07-21 [1] CRAN (R 4.2.0)
#  R.methodsS3              1.8.2     2022-06-13 [1] CRAN (R 4.2.0)
#  R.oo                     1.25.0    2022-06-12 [1] CRAN (R 4.2.0)
#  R.utils                  2.12.2    2022-11-11 [1] CRAN (R 4.2.0)
#  R6                       2.5.1     2021-08-19 [1] CRAN (R 4.2.0)
#  ragg                     1.2.5     2023-01-12 [1] CRAN (R 4.2.0)
#  rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.2.0)
#  RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.2.0)
#  Rcpp                     1.0.10    2023-01-22 [1] CRAN (R 4.2.0)
#  RCurl                    1.98-1.10 2023-01-27 [1] CRAN (R 4.2.0)
#  rematch2                 2.1.2     2020-05-01 [1] CRAN (R 4.2.0)
#  remotes                  2.4.2     2021-11-30 [1] CRAN (R 4.2.0)
#  reshape2                 1.4.4     2020-04-09 [1] CRAN (R 4.2.0)
#  restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.2.0)
#  rhdf5                    2.42.0    2022-11-01 [1] Bioconductor
#  rhdf5filters             1.10.0    2022-11-01 [1] Bioconductor
#  Rhdf5lib                 1.20.0    2022-11-01 [1] Bioconductor
#  rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.2.0)
#  rlang                    1.1.0     2023-03-14 [1] CRAN (R 4.2.0)
#  rmarkdown                2.20      2023-01-19 [1] CRAN (R 4.2.0)
#  rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.2.0)
#  Rsamtools                2.14.0    2022-11-01 [1] Bioconductor
#  RSQLite                  2.3.0     2023-02-17 [1] CRAN (R 4.2.2)
#  rsthemes                 0.3.1     2022-03-01 [1] Github (gadenbuie/rsthemes@bbe73ca)
#  rstudioapi               0.14      2022-08-22 [1] CRAN (R 4.2.0)
#  rsvd                     1.0.5     2021-04-16 [1] CRAN (R 4.2.0)
#  rtracklayer              1.58.0    2022-11-01 [1] Bioconductor
#  S4Vectors              * 0.36.2    2023-02-26 [1] Bioconductor
#  sass                     0.4.5     2023-01-24 [1] CRAN (R 4.2.0)
#  ScaledMatrix             1.6.0     2022-11-01 [1] Bioconductor
#  scales                   1.2.1     2022-08-20 [1] CRAN (R 4.2.0)
#  scater                   1.26.1    2022-11-13 [1] Bioconductor
#  scatterpie               0.1.8     2022-09-03 [1] CRAN (R 4.2.1)
#  scuttle                  1.8.4     2023-01-22 [1] Bioconductor
#  sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.2.0)
#  shadowtext               0.1.2     2022-04-22 [1] CRAN (R 4.2.0)
#  shiny                    1.7.4     2022-12-15 [1] CRAN (R 4.2.2)
#  shinyWidgets             0.7.6     2023-01-08 [1] CRAN (R 4.2.0)
#  SingleCellExperiment   * 1.20.0    2022-11-01 [1] Bioconductor
#  spam                     2.9-1     2022-08-07 [1] CRAN (R 4.2.0)
#  sparseMatrixStats        1.10.0    2022-11-01 [1] Bioconductor
#  SpatialExperiment      * 1.8.1     2023-03-05 [1] Bioconductor
#  spatialLIBD            * 1.11.12   2023-03-17 [1] Github (LieberInstitute/spatialLIBD@4b5d6e5)
#  statmod                  1.5.0     2023-01-06 [1] CRAN (R 4.2.0)
#  stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.2.0)
#  stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.2.0)
#  styler                   1.9.1     2023-03-04 [1] CRAN (R 4.2.0)
#  SummarizedExperiment   * 1.28.0    2022-11-01 [1] Bioconductor
#  suncalc                  0.5.1     2022-09-29 [1] CRAN (R 4.2.0)
#  systemfonts              1.0.4     2022-02-11 [1] CRAN (R 4.2.0)
#  testthat               * 3.1.7     2023-03-12 [1] CRAN (R 4.2.0)
#  textshaping              0.3.6     2021-10-13 [1] CRAN (R 4.2.0)
#  tibble                   3.2.0     2023-03-08 [1] CRAN (R 4.2.0)
#  tidygraph                1.2.3     2023-02-01 [1] CRAN (R 4.2.0)
#  tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.2.0)
#  tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.2.0)
#  tidytree                 0.4.2     2022-12-18 [1] CRAN (R 4.2.0)
#  timechange               0.2.0     2023-01-11 [1] CRAN (R 4.2.0)
#  treeio                   1.22.0    2022-11-01 [1] Bioconductor
#  tweenr                   2.0.2     2022-09-06 [1] CRAN (R 4.2.1)
#  urlchecker               1.0.1     2021-11-30 [1] CRAN (R 4.2.0)
#  usethis                * 2.1.6     2022-05-25 [1] CRAN (R 4.2.0)
#  utf8                     1.2.3     2023-01-31 [1] CRAN (R 4.2.0)
#  vctrs                    0.6.0     2023-03-16 [1] CRAN (R 4.2.2)
#  vipor                    0.4.5     2017-03-22 [1] CRAN (R 4.2.0)
#  viridis                  0.6.2     2021-10-13 [1] CRAN (R 4.2.0)
#  viridisLite              0.4.1     2022-08-22 [1] CRAN (R 4.2.0)
#  withr                    2.5.0     2022-03-03 [1] CRAN (R 4.2.0)
#  xfun                     0.37      2023-01-31 [1] CRAN (R 4.2.0)
#  XML                      3.99-0.13 2022-12-04 [1] CRAN (R 4.2.0)
#  xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.2.0)
#  XVector                  0.38.0    2022-11-01 [1] Bioconductor
#  yaml                     2.3.7     2023-01-23 [1] CRAN (R 4.2.0)
#  yulab.utils              0.0.6     2022-12-20 [1] CRAN (R 4.2.0)
#  zlibbioc                 1.44.0    2022-11-01 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
