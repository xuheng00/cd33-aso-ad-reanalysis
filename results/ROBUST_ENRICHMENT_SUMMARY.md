# Functional characterization of robust gene sets

Gene sets: robust response = meta-FDR < 0.05, concordant direction, I2 < 50%.
Over-representation: hypergeometric test against the 12,272-gene measured space; BH FDR.

## CD33-specific vs control-ASO overlap

         tissue cd33_n control_n overlap overlap_frac_cd33 overlap_frac_control
 Frontal_cortex    245       722      55             0.224                0.076
    Hippocampus    108       915      16             0.148                0.017
 jaccard
   0.060
   0.016

## GO-BP enrichment (CD33-specific, FDR < 0.05)

      tissue                        name set_size query_overlap expected
 Hippocampus     protein phosphorylation      403            10     3.55
 Hippocampus       inflammatory response      205             6     1.80
 Hippocampus               cell adhesion      337             8     2.97
 Hippocampus             phosphorylation      456             9     4.01
 Hippocampus oxidation-reduction process      458             9     4.03
 fold_enrichment           P        FDR
            2.82 0.002880686 0.02880686
            3.33 0.009440699 0.03238402
            2.70 0.009715206 0.03238402
            2.24 0.019034695 0.03904640
            2.23 0.019523201 0.03904640

## GO-BP enrichment (control-ASO background, FDR < 0.05, top 20)

         tissue                                      name set_size
 Frontal_cortex                 cellular response to cAMP       21
 Frontal_cortex               response to toxic substance       40
    Hippocampus                 regulation of cell growth       41
    Hippocampus                 defense response to virus       93
    Hippocampus   immunoglobulin mediated immune response       10
    Hippocampus negative regulation of peptidase activity       56
 query_overlap fold_enrichment         FDR
             7            5.67 0.023700146
             9            3.82 0.037559805
            13            4.25 0.001255942
            19            2.74 0.005751552
             5            6.71 0.035137242
            12            2.87 0.045124463

## KEGG enrichment (CD33-specific, FDR < 0.05)

         tissue       name set_size query_overlap fold_enrichment        FDR
 Frontal_cortex Alcoholism      134             9            3.36 0.04366456
