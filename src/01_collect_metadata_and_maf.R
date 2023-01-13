library(GAMBLR)
library(data.table)
library(tidyverse)

setwd("/projects/rmorin/projects/DLBCL_classification/")


###############################
#### Handle metadata
only_b_dlbcls <- c(
    "16-37587_tumorB",
    "16-12281_tumorB"
)
b_better_dlbcls <- c(
    "02-24492_tumorB",
    "06-11677_tumorB",
    "06-22314_tumorB",
    "13-38657_tumorB",
    "14-20552_tumorB",
    "95-32141_tumorB"
)
b_dlbcl_patients <- c(
    "06-11677",
    "06-22314",
    "02-24492",
    "16-37587",
    "13-38657",
    "16-12281",
    "95-32141",
    "14-20552"
)

metadata.genomes <-
    get_gambl_metadata(seq_type_filter = "genome") %>%
    dplyr::filter(
        !pathology %in% c("CLL", "MCL", "MM", "UNSPECIFIED", "SCBC") &
        !cohort %in% c("DLBCL_ctDNA")
    ) %>%
    dplyr::filter(
        !cohort=="DLBCL_LSARP_Trios" |
        (cohort=="DLBCL_LSARP_Trios" & ! patient_id %in% b_dlbcl_patients & time_point == "A") |
        patient_id %in% b_dlbcl_patients & time_point == "B") %>%
    dplyr::select(
        Tumor_Sample_Barcode,
        COO_consensus,
        DHITsig_consensus,
        pathology,
        seq_type,
        time_point,
        cohort
    ) %>%
    base::merge(
        .,
        read_tsv(
            "data/GAMBL_all_the_things.lymphgen_calls.with_cnvs.with_sv.no_A53.tsv"
            ) %>%
            dplyr::rename(
                "Tumor_Sample_Barcode" = "Sample.Name",
                "lymphgen" = "Subtype.Prediction") %>%
            dplyr::select(Tumor_Sample_Barcode, lymphgen),
        all.x=TRUE) %>%
    dplyr::rename(
        "COO" = "COO_consensus",
        "DHITSIG" = "DHITsig_consensus"
    ) %>%
    tidy_lymphgen(
        .,
        lymphgen_column_in = "lymphgen",
        lymphgen_column_out = "lymphgen"
    ) %>%
    dplyr::mutate(
        pathology = ifelse(pathology == "BLL", "B-ALL", pathology),
        DHITSIG = case_when(
            DHITSIG == "DHITsigNeg" ~ "NEG",
            DHITSIG == "DHITsigPos" ~ "POS",
            DHITSIG == "DHITsig-IND" ~ "IND"
        )
    ) %>%
    dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)

# add inferred metadata
metadata.genomes <-
    GAMBLR:::collate_curated_sv_results(
        get_gambl_metadata()
    ) %>%
    dplyr::filter(
        Tumor_Sample_Barcode %in% metadata.genomes$Tumor_Sample_Barcode
    ) %>%
    dplyr::distinct(
        sample_id,
        patient_id,
        .keep_all = TRUE
    ) %>%
    dplyr::mutate(
        genetic_subgroup = ifelse(
            genetic_subgroup %in% c(
                "DGG-BL", "DLBCL-A", "DLBCL-B",
                "DLBCL-C", "IC-BL", "M53-BL"),
            BL_subgroup_HTMCP,
            genetic_subgroup
        )
    ) %>%
    dplyr::select(
        Tumor_Sample_Barcode,
        genetic_subgroup,
        the_11q_event,
        is_adult,
        EBV_status_inf
    ) %>%
    base::merge(
        metadata.genomes,
        .,
        by = "Tumor_Sample_Barcode"
    )

# ensure genetic_subgroup does not contain lymphgen
metadata.genomes <-
  metadata.genomes %>%
  dplyr::mutate(
    genetic_subgroup = ifelse(
        ! genetic_subgroup %in% c(
            "DGG-BL",
            "IC-BL",
            "M53-BL",
            "Q53-BL",
            "DLBCL-A",
            "DLBCL-B",
            "DLBCL-C",
            "cFL",
            "dFL"
            ),
    NA,
    genetic_subgroup
    )
)

maf_path <- "/projects/rmorin/projects/gambl-repos/gambl-kdreval/results/all_the_things/slms_3-1.0_vcf2maf-1.3/genome--projection/deblacklisted/maf/all_slms-3--grch37.maf"

# get the latest maf
maf <- fread_maf(
    maf_path
)

maf <- maf %>%
  dplyr::filter(
    Tumor_Sample_Barcode %in% metadata.genomes$Tumor_Sample_Barcode
)

maf <- maf[,1:45]

# ensure the samples that don't complete SLMS-3 are not here
metadata.genomes <- metadata.genomes %>%
  dplyr::filter(
    Tumor_Sample_Barcode %in% maf$Tumor_Sample_Barcode
)

# Collect the SV data
annotated_sv <- get_combined_sv() %>%
    dplyr::filter(VAF_tumour >= 0.1) %>%
    rename("SOMATIC_SCORE" = "SCORE") %>%
    annotate_sv(.)

annotated_sv <- annotated_sv %>%
    dplyr::filter(
        tumour_sample_id %in% metadata.genomes$Tumor_Sample_Barcode
    )

# To add manually picked SV and FISH results
collated_data <-
    collate_results(
        join_with_full_metadata = TRUE,
        from_cache = FALSE
    ) %>%
    dplyr::filter(
        Tumor_Sample_Barcode %in% metadata.genomes$Tumor_Sample_Barcode
    )


# get CNV segments for all samples
seg_data <- get_sample_cn_segments(
        multiple_samples = TRUE,
        sample_list = metadata.genomes$Tumor_Sample_Barcode,
        projection = "hg38"
    )

# Save flatfiles for future use
write_tsv(
    metadata.genomes,
    "data/metadata.tsv"
)

write_tsv(
    maf,
    "data/maf.maf"
)

write_tsv(
    collated_data,
    "data/collated_data.tsv"
)

write_tsv(
    annotated_sv,
    "data/annotated_svs.tsv"
)

write_tsv(
    seg_data,
    "data/cnvs.seg"
)
