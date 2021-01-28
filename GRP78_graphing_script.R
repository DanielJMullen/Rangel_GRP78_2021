## Load the BiocManager to install packages from Bioconductor:
if (!require("BiocManager", quietly = TRUE)){

  install.packages("BiocManager")

}

## Load devtools
if (!require("devtools", quietly = TRUE)){

  install.packages("devtools")

}

## Load rtracklayer:
if (!require('rtracklayer')){

  BiocManager::install("rtracklayer")
}

## Load ggplot2
if (!require('ggplot2')){

  BiocManager::install("ggplot2")
}

## Load TCGAbiolinks
if (!require('TCGAbiolinks')){

  BiocManager::install("TCGAbiolinks")
}

## Load EnvStats
if (!require('EnvStats')){

  BiocManager::install("EnvStats")
}

## Load ggpubr
if (!require('ggpubr')){

  BiocManager::install("ggpubr")
}

## Load rstatix
if (!require('rstatix')){

  BiocManager::install("rstatix")
}

## Set the stem for loading files:
stem='D:/'

# ## Read in gene info:
# gtf <- import(
#   paste(
#     stem,
#     'Reference_files/Gene_ensembl_data/gencode.v22.annotation.gtf',
#     sep=''
#   )
# )
# gencode_v22_genes <- as.data.frame(gtf)
# 
# ## Subset gene info to only the whole genes themselves:
# gencode_v22_genes_df <- gencode_v22_genes[gencode_v22_genes$type=='gene',]
# 
# ## Subset the gene IDs to remove the stuff after the periods (and the periods)
# gencode_v22_genes_df$gene_id <- sub(
#   '\\..*',
#   '',
#   gencode_v22_genes_df$gene_id
# )
# 
# ## Set the rownames equal to the genes ENSG:
# rownames(gencode_v22_genes_df) <- gencode_v22_genes_df$gene_id
# 
# ## Read in CpG info:
# hg38_cpg_info <- read.delim(
#   file= paste(
#     stem,
#     'Reference_files/Methylation_probe_annotation/hm450.hg38.manifest.tsv',
#     sep=''
#   ),
#   stringsAsFactors = FALSE
# )
# rownames(hg38_cpg_info) <- hg38_cpg_info$probeID

## Load the .tsv of cbioportal data:
## Get from: https://www.cbioportal.org/results/download?genetic_profile_ids_PROFILE_MUTATION_EXTENDED=luad_tcga_pan_can_atlas_2018_mutations&genetic_profile_ids_PROFILE_COPY_NUMBER_ALTERATION=luad_tcga_pan_can_atlas_2018_gistic&cancer_study_list=luad_tcga_pan_can_atlas_2018&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&data_priority=0&profileFilter=0&case_set_id=luad_tcga_pan_can_atlas_2018_cnaseq&gene_list=KRAS%250AEGFR&geneset_list=%20&tab_index=tab_visualize&Action=Submit
cBioportal_complete <- read.delim(
  "C:/Users/Danie/Downloads/LUAD_TCGA_PanCancerAtlas_cBioPortal_alterations_across_samples_1_27_21.tsv",
  sep='\t',
  header=TRUE,
  stringsAsFactors = FALSE
)

## Set working directory:
setwd('C:/Users/Danie/TCGA_biolinks_data')

## Get LUAD expression data

## Set up the Biolinks query:
LUAD_expression_query <- TCGAbiolinks::GDCquery(
  project = "TCGA-LUAD" ,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy= 'RNA-Seq',
  workflow.type="HTSeq - FPKM-UQ",
  legacy = FALSE
)

## Download the LUAD data
TCGAbiolinks::GDCdownload(LUAD_expression_query)

## Put the TCGA data in a dataframe:
LUAD_expression_data <- TCGAbiolinks::GDCprepare(
  query=LUAD_expression_query,
  summarizedExperiment= FALSE
)
LUAD_expression_data <- as.data.frame(LUAD_expression_data)

## Remove the stuff after the periods from the gene annotation:
LUAD_expression_data$X1 <- sub(
  '\\..*',
  '',
  LUAD_expression_data$X1
)

## Set the rownames to be the ensembl gene annotations, then remove that column:
rownames(LUAD_expression_data) <- LUAD_expression_data[,1]
LUAD_expression_data <- LUAD_expression_data[,-1]

## Log2 transform:
LUAD_expression_data <- log2(LUAD_expression_data+1)

## Select tumor and normal samples only
LUAD_tumor_expression <- LUAD_expression_data[
  as.numeric(
    substring(
      colnames(LUAD_expression_data),
      14,
      15
    )
  )<10
]

LUAD_normal_expression <- LUAD_expression_data[
  as.numeric(
    substring(
      colnames(LUAD_expression_data),
      14,
      15
    )
  )>=10
]

## Remove duplicated tumor samples:

## First get the tumor names and sort them
## alphanumerically and number them:
tumor_names <- colnames(LUAD_tumor_expression)
tumor_names <- sort(tumor_names)
names(tumor_names) <- c(
  1:length(tumor_names)
)

# Truncate the tumor_names to just the participant IDs:
tumor_names_trunc <- substring(
  tumor_names,
  1,
  12
)

## Get the deduplicated patient sample names:
tumor_names_trunc_dedup <- tumor_names_trunc[
  !duplicated(tumor_names_trunc)
]

tumor_names_dedup <- unname(
  tumor_names[
    names(tumor_names_trunc_dedup)
  ]
)

## Now get a dataset of tumor expression from the dedup samples:
LUAD_tumor_expression_dedup <- LUAD_tumor_expression[
  tumor_names_dedup
]

## Sort the normal samples for the heck of it:
LUAD_normal_expression <- LUAD_normal_expression[
  sort(
    colnames(
      LUAD_normal_expression
    )
  )
]

## Get LUAD clinical data
LUAD_clinical_query <- TCGAbiolinks::GDCquery(
  project = "TCGA-LUAD" ,
  data.category = "Clinical",
  file.type = "xml")

TCGAbiolinks::GDCdownload(
  LUAD_clinical_query
)

LUAD_clinical_data <- TCGAbiolinks::GDCprepare_clinic(
  LUAD_clinical_query,
  clinical.info = "patient"
)

LUAD_clinical_data <- unique(LUAD_clinical_data)

## Let's create a dataset of files to create box plots:
boxplot_df <- data.frame(
  'sample_names'= c(
    colnames(LUAD_normal_expression),
    colnames(LUAD_tumor_expression_dedup)
  ),
  'patient_names'= c(
    substring(
      colnames(LUAD_normal_expression),
      1,
      12
    ),
    substring(
      colnames(LUAD_tumor_expression_dedup),
      1,
      12
    )
  ),
  'patient_names_extended'= c(
    substring(
      colnames(LUAD_normal_expression),
      1,
      15
    ),
    substring(
      colnames(LUAD_tumor_expression_dedup),
      1,
      15
    )
  ),
  'GRP78_expression'= c(
    unlist(
      LUAD_normal_expression[
        'ENSG00000044574',
      ]
    ),
    unlist(
      LUAD_tumor_expression_dedup[
        'ENSG00000044574',
      ]
    )
  ),
  stringsAsFactors = FALSE
)

boxplot_df$basic_grouping <- c(
  rep(
    "LUAD_NTL",
    ncol(LUAD_normal_expression)
  ),
  rep(
    "LUAD_Tumor",
    ncol(LUAD_tumor_expression_dedup)
  )
)

## Add just the extended patient names as rownames to the dataset:
rownames(boxplot_df) <- boxplot_df$patient_names_extended

## Add info to the cbioportal mutational dataset on if a KRAS/EGFR alteration
## is present or not:
cBioportal_complete$KRAS_alteration_present <- ifelse(
  cBioportal_complete$KRAS!='no alteration',
  TRUE,
  FALSE
)

cBioportal_complete$EGFR_alteration_present <- ifelse(
  cBioportal_complete$EGFR!='no alteration',
  TRUE,
  FALSE
)

## Now identify samples in which a KRAS/EGFR mutation is present:
cBioportal_complete$KRAS_mutation_only_present <- ifelse(
  cBioportal_complete$KRAS..MUT!='no alteration',
  TRUE,
  FALSE
)

cBioportal_complete$EGFR_mutation_only_present <- ifelse(
  cBioportal_complete$EGFR..MUT!='no alteration',
  TRUE,
  FALSE
)

## Now identify samples in which a G12 KRAS mutation is present:
cBioportal_complete$KRAS_G12_mutation_only_present <- grepl(
  'G12',
  cBioportal_complete$KRAS..MUT
)

## Identify samples with dual KRAS/EGFR alterations and mutations:
cBioportal_complete$KRAS_EGFR_dual_alteration_present <- ifelse(
  cBioportal_complete$KRAS_alteration_present== TRUE,
  ifelse(
    cBioportal_complete$EGFR_alteration_present== TRUE,
    TRUE,
    FALSE
  ),
  FALSE
)

cBioportal_complete$KRAS_EGFR_dual_mutation_present <- ifelse(
  cBioportal_complete$KRAS_mutation_only_present== TRUE,
  ifelse(
    cBioportal_complete$EGFR_mutation_only_present== TRUE,
    TRUE,
    FALSE
  ),
  FALSE
)

## Find samples that have just a KRAS mutation only and no KRAS alteration:
cBioportal_complete$KRAS_mutation_no_EGFR_alterations <- ifelse(
  cBioportal_complete$KRAS_mutation_only_present==TRUE,
  ifelse(
    cBioportal_complete$EGFR_alteration_present==TRUE,
    FALSE,
    TRUE
  ),
  FALSE
)

cBioportal_complete$KRAS_G12_mutation_no_EGFR_alterations <- ifelse(
  cBioportal_complete$KRAS_G12_mutation_only_present==TRUE,
  ifelse(
    cBioportal_complete$EGFR_alteration_present==TRUE,
    FALSE,
    TRUE
  ),
  FALSE
)

cBioportal_complete$EGFR_mutation_no_KRAS_alterations <- ifelse(
  cBioportal_complete$EGFR_mutation_only_present==TRUE,
  ifelse(
    cBioportal_complete$KRAS_alteration_present==TRUE,
    FALSE,
    TRUE
  ),
  FALSE
)

## Let's get the information into the boxplot df and identify samples without
## the mutational data:
rownames(cBioportal_complete) <- cBioportal_complete$Sample.ID

## Let's add mutational info back into the dataframe:
## Add the mutational data into the boxplot_df:
boxplot_df$KRAS_alteration_present <- cBioportal_complete[
  rownames(boxplot_df),
  'KRAS_alteration_present'
]

boxplot_df$EGFR_alteration_present <- cBioportal_complete[
  rownames(boxplot_df),
  'EGFR_alteration_present'
]

boxplot_df$KRAS_mutation_present <- cBioportal_complete[
  rownames(boxplot_df),
  'KRAS_mutation_only_present'
]

boxplot_df$KRAS_mutation_no_EGFR_alterations <- cBioportal_complete[
  rownames(boxplot_df),
  'KRAS_mutation_no_EGFR_alterations'
]

boxplot_df$EGFR_mutation_no_KRAS_alterations <- cBioportal_complete[
  rownames(boxplot_df),
  'EGFR_mutation_no_KRAS_alterations'
]

boxplot_df$KRAS_G12_mutation_no_EGFR_alterations <- cBioportal_complete[
  rownames(boxplot_df),
  'KRAS_G12_mutation_no_EGFR_alterations'
]

## Let's assemble a new dataset with multiple elements for entries in multiple groups

## First group, the Normal samples:
new_boxplot_df <- boxplot_df[
  boxplot_df$basic_grouping=='LUAD_NTL',
  c(1:4)
]

## Then remove the remaining NA samples from the boxplot_df:
boxplot_df_no_NA <- boxplot_df[
  !is.na(boxplot_df$KRAS_alteration_present),
]

## Next let's add the 149 samples with a KRAS mutation and no EGFR mutation:
new_boxplot_df <- rbind(
  new_boxplot_df,
  boxplot_df_no_NA[
    boxplot_df_no_NA$KRAS_mutation_no_EGFR_alterations==TRUE,
    c(1:4)
  ]
)

## Add the 133 samples that are KRAS G12 mutants:
new_boxplot_df <- rbind(
  new_boxplot_df,
  boxplot_df_no_NA[
    boxplot_df_no_NA$KRAS_G12_mutation_no_EGFR_alterations==TRUE,
    c(1:4)
  ]
)

## Add the 62 samples with just an EGFR mutation, no other alteration:
new_boxplot_df <- rbind(
  new_boxplot_df,
  boxplot_df_no_NA[
    boxplot_df_no_NA$EGFR_mutation_no_KRAS_alterations==TRUE,
    c(1:4)
  ]
)

## Now add the 263 samples without a KRAS or EGFR alteration:
new_boxplot_df <- rbind(
  new_boxplot_df,
  boxplot_df_no_NA[
    (boxplot_df_no_NA$KRAS_alteration_present==FALSE & boxplot_df_no_NA$EGFR_alteration_present==FALSE),
    c(1:4)
  ]
)

## Now create a new grouping for this dataset:
new_boxplot_df$grouping <- c(
  rep(
    'Normal',
    59
  ),
  rep(
    'Tumor\nKRAS Mut\nEGFR WT',
    149
  ),
  rep(
    'Tumor\nKRAS G12 Mut\nEGFR WT',
    133
  ),
  rep(
    'Tumor\nKRAS WT\nEGFR Mut',
    62
  ),
  rep(
    'Tumor\nKRAS WT\nEGFR WT',
    263
  )
)

## From the old database, create a tumor vs normal grouping:
boxplot_df$old_grouping <- ifelse(
  substring(
    rownames(boxplot_df),
    14,
    15
  ) > 10,
  'Normal\n\n',
  'Tumor\n\n'
)

## In the new dataset, transform the advanced grouping to a factor:
new_boxplot_df$grouping <- factor(
  new_boxplot_df$grouping,
  levels= c(
    'Normal',
    'Tumor\nKRAS Mut\nEGFR WT',
    'Tumor\nKRAS G12 Mut\nEGFR WT',
    'Tumor\nKRAS WT\nEGFR Mut',
    'Tumor\nKRAS WT\nEGFR WT'
  )
)

boxplot_df$old_grouping <- factor(
  boxplot_df$old_grouping,
  levels= c(
    'Normal\n\n',
    'Tumor\n\n'
  )
)

## Add colors to plot:
grey_scale <- gray.colors(
  n= 5, 
  start = 0, 
  end = 1,
  rev = TRUE
)

advanced_color_fill <- c(
  'Normal'= 'white',
  'Tumor\nKRAS Mut\nEGFR WT'= 'orange',
  'Tumor\nKRAS G12 Mut\nEGFR WT'= 'red',
  'Tumor\nKRAS WT\nEGFR Mut'= 'lightblue',
  'Tumor\nKRAS WT\nEGFR WT'= 'green'
)

basic_color_fill <- c(
  'Normal\n\n'= grey_scale[1],
  'Tumor\n\n'= grey_scale[4]
)

## Add group comparisons:
my_comparisons <- list(
  c("Normal", 'Tumor\nKRAS Mut\nEGFR WT'),
  c("Normal", 'Tumor\nKRAS G12 Mut\nEGFR WT'),
  c("Normal", 'Tumor\nKRAS WT\nEGFR Mut'),
  c("Normal", 'Tumor\nKRAS WT\nEGFR WT')
)

my_comparisons_basic <- list(
  c("Normal\n\n", 'Tumor\n\n')
)

## Perform anova tests comparing the 4 different groups:
mutational_anova <- aov(
  GRP78_expression~grouping,
  data = new_boxplot_df
)
summary(mutational_anova)[[1]][["Pr(>F)"]][[1]]
TukeyHSD(mutational_anova)
mutational_anova_tukey <- TukeyHSD(mutational_anova)

mutational_anova_tukey_rstatix <- tukey_hsd(
  aov(
    GRP78_expression~grouping,
    data = new_boxplot_df
  )
)

overall_mutational_p <- formatC(
  summary(mutational_anova)[[1]][["Pr(>F)"]][[1]],
  format = "e",
  digits = 3
)

## Create the basic plot for comparing tumor vs normal:
test_plot_basic <- ggplot(
  boxplot_df,
  aes(x = old_grouping, y = GRP78_expression)
) +
  geom_boxplot(
    aes(fill = old_grouping)
  )

test_plot_basic +
  ggtitle(
    expression(
      paste(
        "T-test p=", 
        5.648 %*% 10^-49
      ), 
      sep=''
    )
  ) +
  ylab(
    paste(
      "GRP78 expression [log2(FPKM-UQ)]\n"
    )
  ) +
  xlab("Sample Grouping") +
  scale_y_continuous(
    lim=c(20.2,25.5),
    breaks=c(21:24)
  ) +
  guides(fill=FALSE) +
  stat_compare_means(
    comparisons = my_comparisons_basic,
    label= "p.signif",
    size= 8
  ) +
  EnvStats::stat_n_text(
    geom='label',
    size=6
  ) +
  theme_bw() +
  scale_fill_manual(values=basic_color_fill)  +
  theme(
    plot.title = element_text(hjust=0.5, size=20),
    panel.border = element_rect(colour = 'black', fill=NA, size=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18, colour = 'black'),
    axis.text.y = element_text(size=16, colour = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## Create an advanced plot for comparing different mutational groups:
test_plot <- ggplot(
  new_boxplot_df,
  aes(x = grouping, y = GRP78_expression)
) +
geom_boxplot(
  aes(fill = grouping)
)

test_plot +
  ggtitle(
    expression(
      paste(
        "Overall Anova p=", 
        1.139 %*% 10^-32
      ), 
      sep=''
    )
  ) +
  ylab(
    paste(
      "GRP78 expression [log2(FPKM-UQ)]"
    )
  ) +
  xlab("Sample Grouping") +
  scale_y_continuous(
    lim=c(20.2,25.5),
    breaks=c(21:24)
  ) +
  guides(fill=FALSE) +
  stat_compare_means(
    comparisons = my_comparisons,
    label= "p.signif",
    size= 8
  ) +
  EnvStats::stat_n_text(
    geom='label',
    size=6
    ) +
  theme_bw() +
  scale_fill_manual(values=advanced_color_fill)  +
  theme(
    plot.title = element_text(hjust=0.5, size=20),
    panel.border = element_rect(colour = 'black', fill=NA, size=1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=18, colour = 'black'),
    axis.text.y = element_text(size=16, colour = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
