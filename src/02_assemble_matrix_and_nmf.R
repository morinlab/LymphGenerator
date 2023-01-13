# Used to construct UMAP plots and phylogeny tree
library(lsa)
library(ggtree)
library(umap)
# Main packages
library(GAMBLR)
library(NMF)
library(ComplexHeatmap)
library(data.table)
library(tidyverse)

setwd("/projects/rmorin/projects/DLBCL_classification/")


###############################
#### Data definitions #########
###############################

metadata.genomes <- read_tsv(
    "data/metadata.tsv"
)

maf <- fread_maf(
    "data/maf.maf"
)

grch37_ashm_regions <- read_tsv(
      "data/somatic_hypermutation_locations_GRCh37.txt"
  ) %>%
  dplyr::mutate(
      name = paste(
          gene,
          region,
          sep = "-"
      )
  )

lymphoma_genes <- read_tsv(
    "data/lymphoma_genes.tsv"
)


#### Consistent color palette
clinical_colours <- get_gambl_colours("clinical")[c("POS", "NEG")]
COLOUR_2 <- list(
    pathology = get_gambl_colours("pathology")[
      unique(metadata.genomes$pathology)
    ],
    COO = get_gambl_colours("COO")[
      unique(metadata.genomes$COO)[
        !is.na(unique(metadata.genomes$COO))
      ]
    ],
    lymphgen =  get_gambl_colours()[
      unique(metadata.genomes$lymphgen)[
        !is.na(unique(metadata.genomes$lymphgen))
      ]
    ],
    time_point = c(
      "NA" = "#ACADAF",
      get_gambl_colours()[c("A", "B")]
    ),
    DHITSIG = c(
      "IND" = "#1E344B",
      "POS" = "#B82C2C",
      "NEG" = "#A17323",
      "NA" = "#ACADAF"
    ),
    genetic_subgroup = get_gambl_colours()[
      unique(metadata.genomes$genetic_subgroup)[
        !is.na(unique(metadata.genomes$genetic_subgroup))
        ]
      ],
    KMT2D = clinical_colours,
    BCL2_SV = clinical_colours,
    MYC_SV = clinical_colours,
    MYC = clinical_colours
)

# what are the represented pathologies?
ggplot(
  metadata.genomes %>% group_by(pathology) %>% mutate(count = n()),
  aes(x = forcats::fct_infreq(pathology),
      fill = pathology)
) +
  geom_bar() +
  geom_text(aes(y = count + 5, label = count), vjust = 0, size = 8) +
  scale_fill_manual(values = COLOUR_2$pathology) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 450),
    breaks = c(0, 100, 200, 300, 400)
  ) +
  xlab("") + ylab("Number of tumors") +
  theme_Morons(base_size = 18) +
  theme(legend.position = "none")

#########################################################################
### Construct aSHM matrix ###############################################
#### --------------------------------------------------------------------

ashm_counts <- get_ashm_count_matrix(
    regions_bed = grch37_ashm_regions,
    maf_data = maf,
    these_samples_metadata = metadata.genomes %>%
      mutate(sample_id = Tumor_Sample_Barcode)
)

# any samples without ashm? add them to matrix
ashm_counts <- complete_missing_from_matrix(
    ashm_counts,
    metadata.genomes$Tumor_Sample_Barcode
)

ha_bottom <- HeatmapAnnotation(
    df = metadata.genomes %>%
          select(
            Tumor_Sample_Barcode,
            DHITSIG,
            COO,
            genetic_subgroup,
            lymphgen,
            pathology
          ) %>%
          column_to_rownames("Tumor_Sample_Barcode"),
    col = COLOUR_2,
    simple_anno_size = unit(6, "mm"),
    gap = unit(0.25 * 6, "mm"),
    annotation_name_gp = gpar(fontsize = 12),
    annotation_legend_param = list(
          nrow = 1,
          direction = "horizontal"
          )
)

hmap_legend_param <- list(
    title = "Number of mutations",
    legend_direction = "horizontal",
    nrow = 2
)

hMap_ashm <- Heatmap(
    as.matrix(t(ashm_counts)),
    show_column_names = FALSE,
    col = structure(
        c(
            "white",
            (RColorBrewer::brewer.pal(9, "OrRd")[5:9])
        )
    ),
    bottom_annotation = ha_bottom,
    heatmap_legend_param = hmap_legend_param
    )

# draw plot
ComplexHeatmap::draw(
    hMap_ashm,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)

# Binarizing the matrix
ashm_counts <-
    ashm_counts %>%
    rownames_to_column("Tumor_Sample_Barcode") %>%
    base::merge(
        .,
        metadata.genomes %>% select(Tumor_Sample_Barcode, pathology)
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")

# calculate average N of mutations in each aSHM site per pathology
ashm_aggregated <-
    aggregate(
        . ~ pathology,
        data = ashm_counts,
        FUN = mean
    ) %>%
    relocate(
        pathology,
        .after = last_col()
    )

# convert ashm counts and averages to long format
ashm_counts_long <- ashm_counts %>%
    rownames_to_column("Tumor_Sample_Barcode") %>%
    melt(
      .,
      id.vars = c("Tumor_Sample_Barcode", "pathology"),
      variable.name = "Region",
      value.name = "N_mut"
    )

ashm_average_long <- ashm_aggregated %>%
    melt(
        .,
        id.vars = c("pathology"),
        variable.name = "Region",
        value.name = "Average_mut"
    )

# merge counts/averages together
ashm_binary_long <- base::merge(
    ashm_counts_long,
    ashm_average_long,
    by = c("pathology", "Region")
  )

ashm_binary <-
    ashm_binary_long %>%
    mutate(
        Average_mut = ifelse(
            pathology=="DLBCL",
            Average_mut + 3,
            Average_mut + 1)
    ) %>%
    mutate(
        Feature = ifelse(
            N_mut > Average_mut,
            1,
            0)
    ) %>%
    select(
        Tumor_Sample_Barcode,
        Region,
        Feature
    ) %>%
    spread(
        .,
        Region,
        Feature
    ) %>%
    arrange(
        match(
            Tumor_Sample_Barcode,
            metadata.genomes$Tumor_Sample_Barcode)
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")


hMap_binary <-
  Heatmap(
    as.matrix(t(ashm_binary)),
    show_column_names = FALSE,
    cluster_rows = FALSE,
    col = structure(c("white", (rev(
      heat.colors(1)
    )))),
    bottom_annotation = ha_bottom,
    heatmap_legend_param = hmap_legend_param
  )

# draw plot
ComplexHeatmap::draw(
    hMap_binary,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)

ASHM.MATRIX <- ashm_binary
ASHM.MATRIX <- ASHM.MATRIX %>%
    select(
        -c(
            `ST6GAL1-TSS-1`,
            `ST6GAL1-intronic-1`,
            `PAX5-distal-enhancer-2`,
            `ZCCHC7-intron-2`
        )
)

###############################
#### Collect the hotspot regions and construct hotspot matrix
#### ---------------------------------------------------------------------
# ensure genes driven by aSHM are not double-counted in the hotspot matrix
genes_for_hotspots <- sort(
  c(
    "FOXO1",
    "MEF2B",
    "MYD88",
    "STAT6",
    "CD79B",
    "EZH2",
    "CCND3",
    "TP53",
    "CREBBP",
    "NOTCH1",
    "NOTCH2",
    "TNFRSF14"
  )
)

hotspot_long <- annotate_hotspots(
    maf %>% filter(
        Hugo_Symbol %in% genes_for_hotspots),
    recurrence_min = 10
)

hotspot_long <- review_hotspots(hotspot_long) %>%
    filter(hot_spot==TRUE) %>%
    select(Tumor_Sample_Barcode, Hugo_Symbol, hot_spot)

hotspot_binary <- hotspot_long %>%
    mutate(hot_spot=1) %>%
    distinct(
        Tumor_Sample_Barcode,
        Hugo_Symbol,
        .keep_all = TRUE) %>%
    spread(
        .,
        Hugo_Symbol,
        hot_spot,
        fill = 0
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")

# are there any samples without hotspot mutations?
hotspot_binary <- complete_missing_from_matrix(
    hotspot_binary,
    list_of_samples = metadata.genomes$Tumor_Sample_Barcode
)

colnames(hotspot_binary) <- paste0(colnames(hotspot_binary), "HOTSPOT")

hMap_hotspot <-
    Heatmap(
      as.matrix(t(hotspot_binary)),
      show_column_names = FALSE,
      col = structure(c(
        "white",
        (rev(heat.colors(19)))
        )
      ),
      bottom_annotation = ha_bottom,
      heatmap_legend_param = hmap_legend_param
)

# draw plot
ComplexHeatmap::draw(
    hMap_hotspot,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)


HOTSPOT.MATRIX <- hotspot_binary


####################################################################
#### Collect the ssm mutations and construct ssm matrix ############
#### ---------------------------------------------------------------

# define my lymphoma genes

#GAMBLR lymphoma genes, unless not throughout FALSE, unless not known
gamblr_genes <- lymphoma_genes %>%
    select(Gene) %>%
    mutate(lymphoma_genes=1)

chapuy_genes <- read_tsv("data/chapuy_features.tsv") %>%
    dplyr::filter(Description == "Mutation") %>%
    select(Gene=Name) %>%
    mutate(chapuy=1)

lymphgen_genes <-
    annotables::grch37 %>%
    right_join(
      .,
      read_tsv(
        "data/lymphgen_genes.txt",
        col_names = FALSE,
        col_types = "i"
      ),
      by = c("entrez" = "X1")
    ) %>%
    select(Gene=symbol) %>%
    distinct %>%
    mutate(lymphgen=1)

hmrn_features <- read_tsv("data/hmrn_features.tsv") %>%
    separate(Gene, into = "Gene") %>%
    mutate(hmrn=1)

mutsig_genes <-
    read_tsv(
      "data/clustering_genomes_all/clustering_genomes_all.sig_genes.txt"
    ) %>%
    dplyr::filter(q<=0.1) %>%
    select(Gene=gene) %>%
    mutate(mutsig=1)

THESE_GENES <-
    list(gamblr_genes,
        chapuy_genes,
        lymphgen_genes,
        hmrn_features,
        mutsig_genes
    ) %>%
    plyr::join_all(
        .,
        by = "Gene",
        type = "full"
    ) %>%
    arrange(Gene) %>%
    replace(is.na(.), 0) %>%
    mutate(sums = rowSums(select(., -"Gene"))) %>%
    dplyr::filter(sums>1) %>%
    dplyr::filter(
        !grepl("HIST", Gene)|
        (grepl("HIST", Gene) & sums==5)
    ) %>%
    pull(Gene)

# coding mutations for genomes
ssm_binary <-
    get_coding_ssm_status(
        gene_symbols = THESE_GENES,
        these_samples_metadata = metadata.genomes %>%
            mutate(sample_id=Tumor_Sample_Barcode),
        maf_data = maf,
        include_hotspots = FALSE,
        include_silent = FALSE
    ) %>%
    column_to_rownames("sample_id") %>%
    select(
        sort(
          tidyselect::peek_vars()
        )
    )

nfkbiz_mat <- maf %>%
    dplyr::filter(
        Hugo_Symbol=="NFKBIZ" & Variant_Classification=="3'UTR"
    ) %>%
    select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    distinct() %>%
    mutate(
        Hugo_Symbol=paste0(Hugo_Symbol, "_3UTR"),
        mutated=1
    ) %>%
    spread(
        .,
        Hugo_Symbol,
        mutated,
        fill = 0
    ) %>%
    mutate(these_names = Tumor_Sample_Barcode) %>%
    column_to_rownames("these_names")

nfkbiz_mat <- complete_missing_from_matrix(
        nfkbiz_mat,
        list_of_samples = metadata.genomes$Tumor_Sample_Barcode
    ) %>%
    select(-Tumor_Sample_Barcode)

# add NFKBIZ 3'UTR status
ssm_binary <- cbind(
        ssm_binary,
        nfkbiz_mat
    ) %>%
    select(
      sort(
        tidyselect::peek_vars()
      )
    )

ssm_binary <- complete_missing_from_matrix(
    ssm_binary,
    list_of_samples = metadata.genomes$Tumor_Sample_Barcode
)

SSM.MATRIX <- ssm_binary

###############################
#### Collect SV data and construct SV matrix
#### ----------------------------------------------------

annotated_sv <- read_tsv(
    "data/annotated_svs.tsv"
)

annotated_sv <- annotated_sv %>%
    dplyr::filter(
      gene %in% THESE_GENES |
      partner %in% THESE_GENES
    )

this_gene_expression <- get_gene_expression(
      metadata = metadata.genomes %>% mutate(sample_id = Tumor_Sample_Barcode),
      hugo_symbols = c("FOXP1", "MYC", "CIITA", "BCL6"),
      join_with = "genome"
    )

ciita_samples <- annotated_sv %>%
    dplyr::distinct(
        tumour_sample_id,
        fusion,
        .keep_all=TRUE
    ) %>%
    dplyr::filter(fusion=="BCL6-CIITA") %>%
    pull(tumour_sample_id)

this_gene_expression %>%
    dplyr::select(sample_id, pathology, BCL6, CIITA, MYC) %>%
    dplyr::filter(pathology=="DLBCL") %>%
    melt(., id.vars = c("sample_id", "pathology")) %>%
    ggplot(., aes(x =variable,
                  y = value,
                  colour = sample_id %in% ciita_samples,
                  size = sample_id %in% ciita_samples)) +
    ggbeeswarm::geom_quasirandom() +
    facet_wrap( ~ variable, scales = "free") +
    theme_Morons(base_size = 18) +
    scale_colour_discrete(guide = "none") +
    theme(legend.position = "none",
          axis.title.x = element_blank())


# this sample was found to have a complex MYC-BCL6-CIITA event and
# all 3 genes should be annotated
extra_sv <- data.frame(
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  "14-41461T", "BCL6", NA, "CIITA", NA
)
names(extra_sv) <- names(annotated_sv)
annotated_sv <- rbind(annotated_sv, extra_sv)

sv_long <- annotated_sv %>%
    dplyr::filter(
        !fusion %in% c(
            "CIITA-BCL6",
            "FOXP1-BCL6"
        )
    ) %>%
    dplyr::filter(
        !partner=="NA"
    ) %>%
    select(tumour_sample_id, gene) %>%
    distinct(tumour_sample_id, gene) %>%
    group_by(gene) %>%
    mutate(count = n()) %>%
    dplyr::filter(
        count > nrow(metadata.genomes) * 0.002
    )

# how many SV's are there now for each gene?
sv_long %>%
    ggplot(.,
          aes(x = forcats::fct_infreq(gene),
              fill = gene)) +
    geom_bar() +
    geom_text(aes(y = count + 2, label = count), vjust = 0, size = 8) +
    xlab("") + ylab("Number of tumors") +
    theme_Morons(base_size = 18) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# make a binary matrix
sv_binary <- sv_long %>%
  select(-count) %>%
  mutate(mutated = 1) %>%
  spread(
      ., key = gene,
      value = mutated,
      fill = 0
  ) %>%
  rename("Tumor_Sample_Barcode" = "tumour_sample_id")

colnames(sv_binary)[2:ncol(sv_binary)] <-
    paste0(colnames(sv_binary)[2:ncol(sv_binary)], "_SV")

# add manually picked SV and FISH results
collated_data <- read_tsv(
    "data/collated_data.tsv"
)

collated_data <- collated_data %>%
    select(
        Tumor_Sample_Barcode,
        manta_BCL6_sv,
        MYC_SV_any,
        BCL2_SV_any,
        manta_IRF4_sv
    ) %>%
    mutate(
        MYC_SV = ifelse(MYC_SV_any == "POS", 1, 0), # binarize features
        BCL6_SV = ifelse(manta_BCL6_sv == "POS", 1, 0),
        BCL2_SV = ifelse(BCL2_SV_any == "POS", 1, 0),
        IRF4_SV = ifelse(manta_IRF4_sv == "POS", 1, 0)
    ) %>%
    select(-c(
        manta_BCL6_sv,
        MYC_SV_any,
        BCL2_SV_any,
        manta_IRF4_sv
        )
    ) %>% # harmonize names
    unique() %>%
    as.data.frame(.)

# Split MYC SVs by partner
# which samples are FISH-positive, but do not have partner information?
collated_data_test <- collated_data

colnames(collated_data_test)[2:ncol(collated_data_test)] <-
    paste0(
      colnames(collated_data_test)[2:ncol(collated_data_test)],
      "_collated"
    )

myc_others <- full_join(
        collated_data_test,
        sv_binary
    ) %>%
    select(
        Tumor_Sample_Barcode,
        starts_with("MYC")
    ) %>%
    dplyr::filter(MYC_SV_collated > MYC_SV) %>%
    pull(Tumor_Sample_Barcode)


# separate WGS-identified MYC Tx by partner
myc_sv <-
    annotated_sv %>%
    dplyr::filter(
        !fusion %in%  c("CIITA-BCL6", "FOXP1-BCL6")
    ) %>%
    dplyr::filter(!partner == "NA") %>%
    select(tumour_sample_id, gene, partner) %>%
    distinct() %>%
    dplyr::filter(gene=="MYC") %>%
    mutate(
        partner = ifelse(
            ! partner %in% c("IGH", "IGK", "IGL", "BCL6"),
            "OTHER",
            partner
        ),
        partner = paste0(gene, "_", partner)
    ) %>%
    select(-gene) %>%
    mutate(mutated = 1) %>%
    distinct() %>% # sample 10-15025T has MYC Tx to RFTN1 and CCNL1
    spread(., partner, mutated) %>%
    replace(is.na(.), 0)

sv_binary <- full_join(
        sv_binary %>% select(-MYC_SV),
        myc_sv,
        by = c("Tumor_Sample_Barcode" = "tumour_sample_id")
    ) %>%
    replace(is.na(.), 0)

sv_binary <- sv_binary %>%
  left_join(collated_data %>% select(-MYC_SV),
            by = 'Tumor_Sample_Barcode',
            suffix = c(".X", ".Y")) %>%
  na_if(., 0) %>%
  split.default(gsub('.[XY]', '', names(.))) %>%
  map_dfc( ~ if (ncol(.x) == 1)
    .x
    else
      mutate(.x,!!sym(gsub('.X', '', names(
        .x
      )[1])) := coalesce(!!!syms(names(
        .x
      ))))) %>%
  select(!contains(".")) %>%
  column_to_rownames("Tumor_Sample_Barcode") %>%
  replace(is.na(.), 0)


# are there any samples without SVs?
sv_binary <- complete_missing_from_matrix(
    sv_binary,
    list_of_samples = metadata.genomes$Tumor_Sample_Barcode
    )

SV.MATRIX <- sv_binary

hMap_sv <-
  Heatmap(
    as.matrix(t(sv_binary)),
    show_column_names = FALSE,
    col = structure(c("white", (rev(
      heat.colors(1)
    )))),
    bottom_annotation = ha_bottom,
    heatmap_legend_param = hmap_legend_param
  )

# draw plot
ComplexHeatmap::draw(
    hMap_sv,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)


###############################
#### Collect the CNVs and construct cnv matrix
#### --------------------------------------------------------


gistic_peaks <- read_tsv(
    "data/all_lesions.conf_90.txt"
    )

gistic_peaks <- gistic_peaks %>%
    separate(
        `Wide Peak Limits`,
        into=c(
          "chromosome",
          "start",
          "end",
          NA
        )
    ) %>%
    dplyr::select(
        `Unique Name`,
        Descriptor,
        chromosome,
        start,
        end
    ) %>%
    dplyr::filter(!grepl("CN", `Unique Name`)) %>%
    dplyr::mutate(
        end = as.numeric(end),
        start = as.numeric(start),
        length = end-start
    ) %>%
    dplyr::filter(length <= 10000000) %>%
    as.data.table()

setkey(gistic_peaks, chromosome, start, end)

oncogenes_bed <- annotables::grch38 %>%
    dplyr::filter(
        symbol %in% c(THESE_GENES, "NSD2", "MIR17HG"),
        ! grepl("CHR", chr)
    ) %>%
    dplyr::mutate(
        symbol = ifelse(symbol=="NSD2", "WHSC1", symbol),
        chr = paste0("chr", chr)
    ) %>%
    dplyr::select(chr, start, end, symbol) %>%
    dplyr::rename("chromosome" = "chr") %>%
    as.data.table()

setkey(oncogenes_bed, chromosome, start, end)

recurrent_cnv <- foverlaps(
        gistic_peaks,
        oncogenes_bed,
        nomatch = NULL
    ) %>%
    dplyr::arrange(symbol) %>%
    group_by(symbol) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::filter(!n == 2) %>%
    dplyr::rename("Name" = "Unique Name") %>%
    dplyr::mutate(
        feature = ifelse(
            grepl("Amplification", Name),
            paste0(symbol, "_AMP"),
            paste0(symbol, "_LOSS")
        )
    )


# get segments for all samples
seg_data <- read_tsv(
    "data/cnvs.seg"
)

colnames(seg_data)[1] = "sample"

# adjust ploidy
seg_data <- adjust_ploidy(
    this_seg = seg_data,
    projection = "hg38"
)

# drop low-level events
seg_data <- seg_data %>%
    dplyr::mutate(CN = (2 * 2 ^ log.ratio)) %>%
    dplyr::filter(
        abs(log.ratio) > 0.56 &
        !CN %in% c(1, 2, 3)
    ) %>%
    as.data.table()

setkey(seg_data, chrom, start, end)

cnv_coordinates <- oncogenes_bed %>%
    dplyr::filter(
        symbol %in% c(recurrent_cnv$symbol)
    )

setkey(cnv_coordinates, chromosome, start, end)

cnv_binary <- foverlaps(
        cnv_coordinates,
        seg_data
    ) %>%
    dplyr::rename("Tumor_Sample_Barcode" = "sample") %>%
    dplyr::select(Tumor_Sample_Barcode, symbol, CN) %>%
    base::merge(
        .,
        as.data.frame(c(recurrent_cnv$feature)) %>%
            `names<-`("CNV") %>%
            separate(
                CNV,
                c("symbol", "direction")
            ) %>%
            dplyr::mutate(
                direction = ifelse(
                    symbol == "ETS1",
                    "AMP",
                    direction
                  )
            ),
        by = "symbol"
    ) %>%
    dplyr::filter(
        (CN < 2 & direction == "LOSS") |
        (CN > 2 & direction %in% c("AMP", "GAIN"))
    ) %>%
    dplyr::mutate(
        mutated = 1,
        symbol = paste0(symbol, "_", direction)
    ) %>%
    dplyr::select(-c(CN, direction)) %>%
    dplyr::distinct(
        symbol,
        Tumor_Sample_Barcode,
        mutated
    ) %>%
    spread(
        .,
        symbol,
        mutated,
        fill = 0
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode") %>%
    `names<-`(gsub("GAIN", "AMP", colnames(.))) %>%
    as.data.frame(.)


# are there any samples without SVs?
cnv_binary <- complete_missing_from_matrix(
    cnv_binary,
    list_of_samples = metadata.genomes$Tumor_Sample_Barcode
)

CNV.MATRIX <- cnv_binary


hMap_cnv <-
  Heatmap(
    as.matrix(t(cnv_binary)),
    show_column_names = FALSE,
    col = structure(c("white", (rev(
      heat.colors(1)
    )))),
    bottom_annotation = ha_bottom,
    heatmap_legend_param = hmap_legend_param
  )

# draw plot
ComplexHeatmap::draw(
    hMap_cnv,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)

CNV.MATRIX <- CNV.MATRIX %>%
    dplyr::select(
        -c(
          SPIB_AMP,
          NOL9_AMP
        )
    )


###############################
### bind all matrices together and. prepare for NMF input
### double check equal samples in each matrix
dim(ASHM.MATRIX)
dim(SV.MATRIX)
dim(SSM.MATRIX)
dim(HOTSPOT.MATRIX)
dim(CNV.MATRIX)

# DLBCL-only with no FATs and only some CNVs
data <- cbind(
    ASHM.MATRIX,
    SV.MATRIX,
    SSM.MATRIX,
    HOTSPOT.MATRIX,
    CNV.MATRIX[, c(
        "TP53_LOSS",
        "CDKN2A_LOSS",
        "CD58_LOSS",
        "MIR17HG_AMP")]
    ) %>%
    t

# keep all but tese features
consider_ashm_only <- c(
    "MYC_AMP", # consider only TSS aSHM
    "BCL2_AMP", "BCL2_SV", # consider only TSS aSHM
    "BCL6-MUTorAMP", # consider only aSHM/SV
    "EBF1", "SGK1",  # consider only TSS aSHM
    "KMT2D", "KMT2D-MUTorLOSS", # because it is overrepresented
    "BCL11A_AMP",  # redundant
    "ZCCHC7-intronic",
    "LPP-intronic",
    "PAX5-TSS-1",
    "PAX5-intron-1",
    "PAX5-distal-enhancer-3",
    "ZCCHC7-intron-3",
    "DTX1",
    "MKI67",  # physiologic or non-informative for clustering
    "IGLL5",
    "IGLL5-TSS",
    "IGLL5-MUTorLOSS",
    "IL16",
    "HIST1H1B", # draw their own cluster without biological informance
    "HIST1H1E",
    "HIST1H2BK",
    "HIST1H3B",
    "HIST1H1C",
    "HIST1H1D",
    "FAT4", # late replicating genes
    "FAT1"
)

data <- massage_matrix_for_clustering(
        data,
        blacklisted_cnv_regex = "3UTR|SV|HOTSPOT|MYC|BCL2|TP53BP1|intronic",
        drop_these_features = consider_ashm_only,
        min_feature_percent = 0.002
    )

data_dlbcl <- data[,
                  intersect(
                    pull(metadata.genomes %>%
                            dplyr::filter(
                                pathology=="DLBCL"),
                        Tumor_Sample_Barcode),
                    colnames(data)
                  )]

# squishing together some features
bcl2_feat <- c("BCL2-intronic", "BCL2-TSS", "BCL2")
myc_feat <- c("MYC-TSS", "MYC")
bcl6_feat <- c("BCL6_SV", "BCL6")
pim1_feat <- c("PIM1", "PIM1-TSS")
btg2_feat <- c("BTG2", "BTG2-intronic")

# BCL2
data_dlbcl <- data_dlbcl %>%
    as.data.frame %>%
    t() %>%
    as.data.frame %>%
    dplyr::mutate(across({{ bcl2_feat }}, ~replace(., . == 0, NA))) %>%
    dplyr::mutate(BCL2_any = coalesce(!!! syms (bcl2_feat))) %>%
    dplyr::select(-{{bcl2_feat}}) %>%
    dplyr::select(sort(tidyselect::peek_vars())) %>%
    t

# MYC
data_dlbcl <- data_dlbcl %>%
    as.data.frame %>%
    t() %>%
    as.data.frame %>%
    dplyr::mutate(across({{ myc_feat }}, ~replace(., . == 0, NA))) %>%
    dplyr::mutate(MYC_aSHM = coalesce(!!! syms (myc_feat))) %>%
    dplyr::select(-{{myc_feat}}) %>%
    dplyr::select(sort(tidyselect::peek_vars())) %>%
    t

# BCL6
data_dlbcl <- data_dlbcl %>%
    as.data.frame %>%
    t() %>%
    as.data.frame %>%
    dplyr::mutate(across({{ bcl6_feat }}, ~replace(., . == 0, NA))) %>%
    dplyr::mutate(BCL6_any = coalesce(!!! syms (bcl6_feat))) %>%
    dplyr::select(-{{bcl6_feat}}) %>%
    dplyr::select(sort(tidyselect::peek_vars())) %>%
    t


# PIM1
data_dlbcl <- data_dlbcl %>%
    as.data.frame %>%
    t() %>%
    as.data.frame %>%
    dplyr::mutate(across({{ pim1_feat }}, ~replace(., . == 0, NA))) %>%
    dplyr::mutate(PIM1_any = coalesce(!!! syms (pim1_feat))) %>%
    dplyr::select(-{{pim1_feat}}) %>%
    dplyr::select(sort(tidyselect::peek_vars())) %>%
    t


# BTG2
data_dlbcl <- data_dlbcl %>%
    as.data.frame %>%
    t() %>%
    as.data.frame %>%
    dplyr::mutate(across({{ btg2_feat }}, ~replace(., . == 0, NA))) %>%
    dplyr::mutate(BTG2_any = coalesce(!!! syms (btg2_feat))) %>%
    dplyr::select(-{{btg2_feat}}) %>%
    dplyr::select(sort(tidyselect::peek_vars())) %>%
    t

data_dlbcl <- data_dlbcl %>%
  replace(is.na(.), 0)

# Ensure after feature squishing there are none with less than 2% frequency
data_dlbcl <- massage_matrix_for_clustering(
    data_dlbcl,
    blacklisted_cnv_regex = "3UTR|SV|HOTSPOT|MYC|BCL2|TP53BP1|intronic",
    drop_these_features = consider_ashm_only,
    min_feature_percent = 0.002
    )


write_tsv(
    data_dlbcl %>%
    as.data.frame %>%
    rownames_to_column("Feature"),
    "results/nmf_matrix.tsv"
)


# Estimating factorization rank using different techniques
nrun <- 20
seed <- 123456
run_nmf <- function(
        method
    ){
      return(
        nmf(
          data_dlbcl,
          3:12,
          nrun = nrun,
          method = method,
          seed = seed,
          .opt='vp72'
        )
      )
}
estim_lee <- run_nmf(method="lee")
estim_brunet <- run_nmf(method="brunet")
estim_snmf_r <- run_nmf(method="snmf/r")
estim_snmf_l <- run_nmf(method="snmf/l")
estim_offset <- run_nmf(method="offset")

plot(estim_snmf_l)

# After analysis of different methods,
# the Brunet 7-cluster solution is the most stable
num_clusters <- 7

res.multirun <- nmf(
        data_dlbcl,
        num_clusters,
        nrun = nrun * 10,
        method = 'brunet',
        seed = seed,
        .opt = 'vp72'
    )

res.multirun
fit(res.multirun)

initial_order <- c(
    "MP3", "BNZ", "EGB",
    "aSEL", "aSCI", "ETB",
    "MCaP"
)

COLOUR <- c(
    list(
      cluster = get_gambl_colours("lymphgenerator")[initial_order] %>%
          `names<-`(paste0("cluster_", 1:num_clusters)),
      cluster_name = get_gambl_colours("lymphgenerator")[initial_order],
      "MYC partner" = c(
        BCL6 = "#264653",
        IGH = "#E76F51",
        IGK = "#F4A261",
        IGL = "#2A9D8F",
        "0" = "#ACADAF",
        OTHER = "#E9C46A"
      )
    ),
    COLOUR_2)

# get all features
w <- basis(res.multirun)
nmfw <- as.data.frame(w) %>%
    mutate_if(is.character, as.numeric)

# get all sample clusters
h <- NMF::coef(res.multirun) %>% as.data.frame()
h$cluster <- paste0("cluster_", 1:num_clusters)
h <- h %>% dplyr::select(cluster, everything())
h <- melt(h, id.vars=c("cluster")) %>%
    `names<-`(c("cluster", "pid", "value"))


# get which cluster samples belong to
NMF.data <- h %>%
    group_by(pid) %>%
    dplyr::filter(value == max(value)) %>%
    dplyr::select(pid, cluster)
NMF.data <- NMF.data %>%
    dplyr::distinct(pid, .keep_all = TRUE)

# make annotations to display on feature heatmap
anno_col_2 <- base::merge(
        metadata.genomes,
        NMF.data,
        by.x="Tumor_Sample_Barcode",
        by.y="pid"
    ) %>%
    dplyr::mutate(
        DHITSIG=case_when(
          DHITSIG == "DHITsigNeg" ~ "NEG",
          DHITSIG == "DHITsigPos" ~ "POS",
          DHITSIG == "DHITsig-IND" ~ "IND",
          TRUE ~ DHITSIG
        ),
        COO=toupper(COO)
    ) %>%
    base::merge(
        .,
        SSM.MATRIX %>%
          rownames_to_column("Tumor_Sample_Barcode") %>%
          dplyr::select(Tumor_Sample_Barcode, KMT2D, MYC) %>%
          mutate(
            KMT2D=ifelse(KMT2D==0, "NEG", "POS"),
            MYC=ifelse(MYC==0, "NEG", "POS")
          )
    ) %>%
    base::merge(
        .,
        SV.MATRIX %>%
          rownames_to_column("Tumor_Sample_Barcode") %>%
          dplyr::select(Tumor_Sample_Barcode, BCL2_SV) %>%
          mutate(
            BCL2_SV = ifelse(
              BCL2_SV==0,
              "NEG",
              "POS"
            )
          )
    ) %>%
    base::merge(
        .,
        SV.MATRIX %>%
          as.data.frame %>%
          transmute_all(funs(ifelse(. == 1,
                                    deparse(substitute(.)),
                                    NA)
                            )
                        ) %>%
          dplyr::mutate(across(everything(),~ gsub("MYC_", "", .))) %>%
          dplyr::select(starts_with("MYC")) %>%
          dplyr::mutate("MYC partner" = coalesce(
                                            MYC_BCL6,
                                            MYC_IGH,
                                            MYC_IGK,
                                            MYC_IGL,
                                            MYC_OTHER)
                                        ) %>%
          dplyr::select(`MYC partner`) %>%
          rownames_to_column("Tumor_Sample_Barcode")
    ) %>%
    relocate(seq_type, .after = last_col())


#########################################################
## Weird stuff to allow for zooming in of the heatmap
## Why does Complexheatmap not have it as a native option?
##########################################################

# extract most important features, while taking
# the feature with highest weight
# for a particular cluster if it was seen before
# for other cluster with lower weight

n_feat <- 15
FEATURES <- nmfw[,1] %>%
    as.data.frame() %>%
    `rownames<-`(rownames(nmfw)) %>%
    arrange(desc(.)) %>%
    head(., n_feat) %>%
    rownames_to_column(., var="Feature") %>%
    mutate(cluster=1)

for (i in 2:num_clusters){
    FEATURES <- rbind(
      as.data.frame(FEATURES),
      nmfw[,i] %>%
          as.data.frame() %>%
          `rownames<-`(rownames(nmfw)) %>%
          arrange(desc(.)) %>%
          head(., n_feat) %>%
          rownames_to_column(., var="Feature") %>%
          mutate(cluster=i)
      ) %>%
      group_by(Feature) %>%
      dplyr::filter(. == max(.)) %>%
      arrange(cluster)
}
FEATURES <- as.data.frame(FEATURES)


# data mungling to get proper cohorts (EBV+/- BL)
ALL.DATA <- data_dlbcl %>%
    as.data.frame(.) %>%
    t(.)
ALL.DATA <- ALL.DATA %>%
    as.data.frame(.) %>%
    mutate_if(
        is.character,
        as.numeric
    )
data <- ALL.DATA %>%
    rownames_to_column(., var="pid") %>%
    as.data.frame(.)


#MASTER.METADATA <- metadata.genomes


# annotations to be associated with feature heatmap
#ANNO <- metadata.genomes
ANNO <- base::merge(
      metadata.genomes %>%
        dplyr::rename("pid"="Tumor_Sample_Barcode"),
      NMF.data,
      by="pid"
    )

# breaks used to display clusters with different colors on heatmap
bk <- c(0, rep(0:num_clusters+0.5))

# colors to display on heatmap, corresponding to each break
my_palette <- c("white", unname(COLOUR$cluster))
# get each cluster and label the events for each feature with cluster number
data_2 <- data
# subset samples of each cluster
MY.LIST <- list()
for (i in 1:num_clusters){
    name <- paste("cluster", i, sep = "")
    MY.LIST[[i]] <-
        assign(
          name,
          data_2 %>%
            as.data.frame(.) %>%
            column_to_rownames(., var="pid") %>%
            t(.) %>%
            as.data.frame(.) %>%
            dplyr::select(
              ANNO %>%
              dplyr::filter(cluster==paste0("cluster_", i)) %>%
              dplyr::select(pid) %>%
              unlist(list(.)) %>%
              unname(.)
            )
        )
}

# assign cluster number
for(i in 1:num_clusters){
  MY.LIST[[i]][MY.LIST[[i]]>0] <- i
}

# bind them all together for plotting
data_2 <- do.call(
        cbind,
        MY.LIST
    ) %>%
    as.data.frame(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(., var="pid")


# specify where row breaks should be on heatmap
breaks <- 0
for (i in 1:max(FEATURES$cluster)){
    N <- (nrow(FEATURES %>% dplyr::filter(cluster==i)))
    breaks <- c(breaks, N+max(breaks))
}

# second, make a vector that will be supplied to ComplexHeatmap
my_vector <- NULL
for (i in 1:(length(breaks)-1)){
    my_vector <- c(
            my_vector,
            rep(i, breaks[i+1]-breaks[i])
        )
}

# prepare matrix for stacked barplots on the left of the heatmap
STACKED <- data.frame(matrix(NA, ncol=1, nrow=nrow(FEATURES)))[-1]

rownames(STACKED) <- FEATURES$Feature

for (i in 1:num_clusters) {
    STACKED <- cbind(
        STACKED,
        data_2[,c("pid", FEATURES$Feature)] %>%
            base::merge(
              .,
              ANNO %>%
                dplyr::select(pid, cluster),
              by="pid") %>%
            dplyr::arrange(cluster) %>%
            dplyr::filter(
              cluster==paste0("cluster_", i)
            ) %>%
            dplyr::select(-pid, -cluster) %>%
            summarise_all(funs(sum)) %>%
            t(.) %>%
            `colnames<-`(paste0("cluster_",i)) %>%
            as.data.frame(.) %>%
            mutate_all(~(./i)/nrow(data))
    )
}
# normalize it so the values are scaled and sum up to 1
m <- t(apply(STACKED, 1, function(x) x/sum(x)))

# Ensure for consistent ordering of samples on the heatmap
this_is_ordered_df <-
    (anno_col_2 %>%
      column_to_rownames("Tumor_Sample_Barcode")
    )[ order(
        match(
          rownames(
            anno_col_2 %>%
              column_to_rownames("Tumor_Sample_Barcode")
          ),
          colnames(
            t(
              base::merge(
                data_2,
                ANNO,
                by="pid"
              ) %>%
              column_to_rownames("pid") %>%
              dplyr::arrange(cluster) %>%
              dplyr::select(FEATURES$Feature)
            )
          )
        )
      ), ] %>%
    dplyr::select(
        COO, DHITSIG, lymphgen,
        pathology, time_point, BCL2_SV,
        `MYC partner`, KMT2D, cluster
    ) %>%
    dplyr::arrange(
        cluster, COO, DHITSIG,
        lymphgen, pathology, time_point,
        BCL2_SV, `MYC partner`, KMT2D
    )


# left annotation: stacked feature weights
ha <- rowAnnotation(
    `feature abundance` = anno_barplot(
              m,
              gp = gpar(
                fill = COLOUR[["cluster"]]
              ),
              bar_width = 1,
              width = unit(4, "cm"),
              axis_param = list(
                side = "bottom",
                labels_rot = 0
              )
            )
      )

# bottom annotation: tracks indicating metadata
ha_bottom <- HeatmapAnnotation(
    df = this_is_ordered_df %>%
        dplyr::select(-cluster),
    col = COLOUR,
    simple_anno_size = unit(3, "mm"),
    gap = unit(0.25*3, "mm"),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(
        nrow = 1,
        direction = "horizontal"
    )
  )

# top annotation: identified clusters
ha_top <- HeatmapAnnotation(
    df = this_is_ordered_df %>%
        dplyr::select(cluster),
    col = COLOUR["cluster"],
    simple_anno_size = unit(8, "mm"),
    gap = unit(0.25*8, "mm"),
    annotation_name_gp = gpar(fontsize = 15),
    annotation_legend_param = list(
        nrow = 1,
        direction = "horizontal"
    )
  )

ch6 <- ComplexHeatmap::Heatmap(
    t(
      base::merge(
        data_2,
        ANNO,
        by = "pid") %>%
      column_to_rownames("pid") %>%
      dplyr::arrange(cluster) %>%
      dplyr::select(FEATURES$Feature)
    )  %>%
    `rownames<-`(
      gsub('CNV','',rownames(.))
    ) %>%
    as.data.frame() %>%
    dplyr::select(
      rownames(this_is_ordered_df)
    ) %>%
    as.matrix(),
    col = my_palette,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_names_gp = gpar(fontsize = 10),
    show_heatmap_legend = FALSE,
    row_split = my_vector,
    row_title = NULL,
    left_annotation = ha,
    bottom_annotation = ha_bottom,
    top_annotation = ha_top,
    column_split = pull(this_is_ordered_df, cluster),
    column_title = initial_order
  )


cairo_pdf(
    "results/Dconstruct.pdf",
    width = 17,
    height = 22
)

draw(
    ch6,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)

dev.off()
##################################################################
#### End of the NMF results



##################################################################
### Now test differential mutation frequency for each cluster ####
##################################################################
translate_cluster <- function(data, input_col, output_col) {

    mutate(data,
        !!output_col := case_when(
        !!sym(input_col) == "cluster_1" ~ "MP3",
        !!sym(input_col) == "cluster_2" ~ "BNZ",
        !!sym(input_col) == "cluster_3" ~ "EGB",
        !!sym(input_col) == "cluster_4" ~ "aSEL",
        !!sym(input_col) == "cluster_5" ~ "aSCI",
        !!sym(input_col) == "cluster_6" ~ "ETB",
        !!sym(input_col) == "cluster_7" ~ "MCaP"
        )
  )
}

res_meta <- anno_col_2 %>%
    translate_cluster(
      "cluster",
      "cluster_name"
    ) %>%
    mutate(sample_id = Tumor_Sample_Barcode)

res_mutmat <- data_dlbcl %>%
  as.data.frame %>%
  t %>%
  as.data.frame %>%
  rownames_to_column("sample_id")

my_ftest_function_all_possible <- function(
        this_meta,
        this_mutmat,
        this_comparison_column,
        this_max_q
    ){

    levels <- this_meta %>%
        pull({{ this_comparison_column }}) %>%
        unique()

    possible_combinations <- apply(
                                    combn(levels, 2),
                                    2,
                                    paste,
                                    collapse = '--'
                                  )

    message(possible_combinations)

    results <- list()

    for (this_combination in possible_combinations){
        message(
          paste(
            "Processsing the comparison",
            this_combination,
            "..."
          )
        )

        comparison_meta <- this_meta

        these_groups <- unlist(
          strsplit(
            this_combination,
            split = "--"
          )
        )

        print(these_groups)

        results[[this_combination]] <- prettyForestPlot(
          mutmat = this_mutmat,
          metadata = comparison_meta,
          comparison_column = this_comparison_column,
          comparison_values = these_groups,
          comparison_name = this_combination,
          genes = colnames(
            this_mutmat %>%
              dplyr::select(
                -{{ this_comparison_column }},
                -sample_id
              )
          ),
          max_q = this_max_q
        )
    }
    return(results)
}


my_ftest_function_groups <- function(
        this_meta,
        this_mutmat,
        this_comparison_column,
        this_max_q
    ){

    levels <- this_meta %>%
        pull({{ this_comparison_column }}) %>%
        unique()

    message(levels)

    results = list()

    for (this_level in levels){
        message(
          paste(
            "Processsing the comparison",
            this_level,
            "..."
          )
        )

        comparison_meta <- this_meta %>%
          dplyr::mutate(
              {{ this_comparison_column }} := ifelse(
                base::get(this_comparison_column) == this_level,
                this_level,
                paste0("non-", this_level)
              )
            )

        results[[this_level]] <- prettyForestPlot(
          mutmat = this_mutmat,
          metadata = comparison_meta,
          comparison_column = this_comparison_column,
          comparison_values = c(
            this_level,
            paste0("non-",this_level)
          ),
          comparison_name = this_level,
          genes = colnames(
            this_mutmat %>%
            dplyr::select(
              -{{this_comparison_column}},
              -sample_id
            )
          ),
          max_q = this_max_q
        )
    }

    return(results)
}

pairwise_comparisons_all <-  my_ftest_function_groups(
        res_meta,
        res_mutmat,
        "cluster_name",
        this_max_q = 0.1
    )

pairwise_comparisons_pairwise <-  my_ftest_function_all_possible(
        res_meta,
        res_mutmat,
        "cluster_name",
        this_max_q = 0.1
    )

arranged_plots <- sapply(
        pairwise_comparisons_all,
        `[`,
        5
    )

arranged_plots_pairwise <- sapply(
        pairwise_comparisons_pairwise,
        `[`,
        5
    )

cairo_pdf(
    file = "results/Dcounstruct_fisher.pdf",   # The directory you want to save the file in
    width = 10,
    height = 7.5,
    onefile = TRUE
)

for (i in 1:length(arranged_plots)){
  print(arranged_plots[i])
}

for (i in 1:length(arranged_plots_pairwise)){
  print(arranged_plots_pairwise[i])
}

dev.off()


#############################################################
# Reorder levels to place some LymphGen groups side-by-side #
#############################################################

reorder_clusters <- c("2", "3", "6", "4", "5", "1", "7")
reorder_levels <- c("BNZ", "EGB", "ETB", "aSEL", "aSCI", "MP3", "MCaP")

FEATURES_relevel <- FEATURES %>%
    dplyr::mutate(
        cluster = factor(
          cluster,
          levels = reorder_clusters
        )
    ) %>%
    dplyr::mutate(
        cluster = as.numeric(cluster)
    ) %>%
    dplyr::arrange(cluster)

# specify where row breaks should be on heatmap
breaks <- 0
for (i in 1:max(FEATURES_relevel$cluster)){
    N <- (nrow(FEATURES_relevel %>%
                dplyr::filter(cluster==i)
              )
          )
    breaks <- c(breaks, N+max(breaks))
}

# second, make a vector that will be supplied to ComplexHeatmap
my_vector <- NULL
for (i in 1:(length(breaks)-1)){
    my_vector <- c(
      my_vector,
      rep(i, breaks[i+1]-breaks[i])
    )
}

# prepare matrix for stacked barplots on the left of the heatmap
STACKED <- data.frame(matrix(NA, ncol = 1, nrow = nrow(FEATURES_relevel)))[-1]

rownames(STACKED) <- FEATURES_relevel$Feature

for (i in as.numeric(reorder_clusters)) {
    STACKED <- cbind(
        STACKED,
        data_2[,c("pid", FEATURES_relevel$Feature)] %>%
            base::merge(
              .,
              ANNO %>%
                dplyr::select(pid, cluster),
              by="pid"
            ) %>%
            dplyr::arrange(cluster) %>%
            dplyr::filter(
              cluster==paste0("cluster_", i)
            ) %>%
            dplyr::select(-pid, -cluster) %>%
            summarise_all(funs(sum)) %>%
            t(.) %>%
            `colnames<-`(paste0("cluster_", i)) %>%
            as.data.frame(.) %>%
            mutate_all(~(./i)/nrow(data))
        )
}

# normalize it so the values are scaled and sum up to 1
m <- t(apply(STACKED, 1, function(x) x/sum(x)))

# left annotation: stacked feature weights
ha_relevel <- rowAnnotation(
    `feature abundance` = anno_barplot(
      m,
      gp = gpar(
        fill = COLOUR[["cluster"]][paste0(
          "cluster_",
          reorder_clusters
        )]
      ),
      bar_width = 1,
      width = unit(4, "cm"),
      axis_param = list(
        side = "bottom",
        labels_rot = 0
      )
    )
  )

this_is_ordered_df_relevel <-
    (anno_col_2 %>% column_to_rownames("Tumor_Sample_Barcode"))[
      order(
        match(
          rownames(
            anno_col_2 %>%
            column_to_rownames("Tumor_Sample_Barcode")
          ),
          colnames(
            t(
              base::merge(
                data_2,
                ANNO,
                by = "pid"
              ) %>%
              column_to_rownames("pid") %>%
              dplyr::arrange(cluster) %>%
              dplyr::select(FEATURES_relevel$Feature)
            )
          )
        )
      ),] %>%
  dplyr::select(
    COO,
    DHITSIG,
    lymphgen,
    pathology,
    time_point,
    BCL2_SV,
    `MYC partner`,
    KMT2D,
    cluster
  ) %>%
  dplyr::mutate(
    cluster = factor(
      cluster,
      levels=c(
        paste0(
          "cluster_",
          reorder_clusters
        )
      )
    )
  ) %>%
  dplyr::arrange(
    cluster,
    COO,
    DHITSIG,
    lymphgen,
    pathology,
    time_point,
    BCL2_SV,
    `MYC partner`,
    KMT2D
  )

# bottom annotation: tracks indicating metadata
ha_bottom_relevel <- HeatmapAnnotation(
    df = this_is_ordered_df_relevel %>%
      dplyr::select(-cluster),
    col = COLOUR,
    simple_anno_size = unit(3, "mm"),
    gap = unit(0.25*3, "mm"),
    annotation_name_gp = gpar(fontsize = 12),
    annotation_legend_param = list(
      nrow = 2,
      direction = "horizontal")
    )

# top annotation: identified clusters
ha_top_relevel <- HeatmapAnnotation(
    df = this_is_ordered_df_relevel %>%
      dplyr::select(cluster),
    col = COLOUR["cluster"],
    simple_anno_size = unit(8, "mm"),
    gap = unit(0.25*8, "mm"),
    annotation_name_gp = gpar(fontsize=15),
    annotation_legend_param = list(
      nrow = 2,
      direction = "horizontal")
    )

ch6_relevel <- ComplexHeatmap::Heatmap(
    t(
      base::merge(
        data_2,
        ANNO,
        by = "pid") %>%
      column_to_rownames("pid") %>%
      dplyr::mutate(
        cluster = factor(
          cluster,
          levels = paste0(
            "cluster_",
            reorder_clusters
          )
        )
      ) %>%
      dplyr::arrange(cluster) %>%
      dplyr::select(FEATURES_relevel$Feature)
    ) %>%
    `rownames<-`(gsub('CNV','',rownames(.))) %>%
    as.data.frame() %>%
    dplyr::select(
      rownames(this_is_ordered_df_relevel)
    ) %>%
    as.matrix(),
    col = my_palette,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_names_gp = gpar(fontsize = 15),
    show_heatmap_legend = FALSE,
    row_split = my_vector,
    row_title = NULL,
    left_annotation = ha_relevel,
    bottom_annotation = ha_bottom_relevel,
    top_annotation = ha_top_relevel,
    column_split = pull(this_is_ordered_df_relevel, cluster),
    column_title = reorder_levels
  )

cairo_pdf(
    "results/Dconstruct_relevel.pdf",
    width = 17,
    height = 22
)
draw(
    ch6_relevel,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
)
dev.off()


#######################################################
# Save resulting labels ###############################
#######################################################

anno_col_2 %>%
    translate_cluster(
        "cluster",
        "Dconstruct_label"
    ) %>%
    as.data.frame %>%
    write_tsv("results/Dconstruct_labels.tsv")

########################################################
######### Exploratory plots ############################
########################################################

# umap plot
normalize_coefficients <- NMF::coef(res.multirun) %>%
    as.data.frame() %>%
    t %>%
    as.data.frame %>%
    `colnames<-`(initial_order)

for_umap_space <-
    (normalize_coefficients/rowSums(normalize_coefficients)) %>%
    as.data.frame() %>%
    rownames_to_column("Tumor_Sample_Barcode") %>%
    left_join(
        .,
        res_meta %>%
            select(Tumor_Sample_Barcode, cluster_name)
    ) %>%
    column_to_rownames("Tumor_Sample_Barcode")

set.seed(142)
umap_fit <-
    for_umap_space %>%
    select(where(is.numeric)) %>%
    scale() %>%
    umap::umap()

umap_df <- umap_fit$layout %>%
    as.data.frame()%>%
    rename(
        UMAP1 = "V1",
        UMAP2 = "V2"
    ) %>%
    rownames_to_column("Tumor_Sample_Barcode") %>%
    inner_join(
        for_umap_space %>%
            rownames_to_column("Tumor_Sample_Barcode"),
        by = "Tumor_Sample_Barcode"
    )

umap_plot <- umap_df %>%
    ggplot(
        aes(
            x = UMAP1,
            y = UMAP2,
            color = cluster_name
        )
    )+
    geom_point(
        size = 2,
        alpha = 0.75
    ) +
    labs(
        x = "UMAP1",
        y = "UMAP2"
    ) +
    scale_color_manual(
        values = get_gambl_colours("lymphgenerator"),
        name = "Genetic subgroup") +
    theme_Morons()

umap_plot


# Phylogeny tree
for_similarity <-
    for_umap_space %>%
    select(where(is.numeric)) %>%
    as.matrix %>%
    t

similarity_matrix <- cosine(for_similarity)

# dissimilarity matrix
my_dist_mat <- 1-similarity_matrix

# convert to distances
my_dist_mat2 <- as.dist(my_dist_mat)

# building tree
my_nj <- phangorn::upgma(my_dist_mat2)

d <- my_nj$tip.label %>%
    as.data.frame %>%
    `colnames<-`("Tumor_Sample_Barcode") %>%
    left_join(
        .,
        umap_df %>%
            select(
                Tumor_Sample_Barcode,
                cluster_name
            )
        )

tree <- left_join(
        my_nj,
        d,
        by=c('label' = "Tumor_Sample_Barcode")
    )

phylo_plot <- ggtree(
        tree,
        branch.length = 'none',
        aes(color=cluster_name)
    ) +
    scale_color_manual(
        values=get_gambl_colours("lymphgenerator")
    ) +
    labs(color = "Genetic subgroup") +
    theme_tree2()

phylo_plot
