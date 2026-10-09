# Human AD external validation

## Included sample counts

    series         region  status Freq
1 GSE48350 Frontal_cortex      AD   21
2  GSE5281 Frontal_cortex      AD   23
3 GSE48350    Hippocampus      AD   19
4  GSE5281    Hippocampus      AD   10
5 GSE48350 Frontal_cortex Control   26
6  GSE5281 Frontal_cortex Control   11
7 GSE48350    Hippocampus Control   25
8  GSE5281    Hippocampus Control   13

## Cohort-level validation

    series         tissue n_mapped spearman_rho      cor_P down_N
1 GSE48350    Hippocampus       94   0.22186179 0.03162572     37
2 GSE48350 Frontal_cortex      215   0.08026516 0.24121745     66
3  GSE5281    Hippocampus       94  -0.01396525 0.89372442     37
4  GSE5281 Frontal_cortex      215   0.04940110 0.47116817     66
  down_human_direction     down_P up_N up_human_direction       up_P
1                 Down 0.10049295   57               Down 0.14368148
2                 Down 0.01098908  149               Down 0.05226770
3                   Up 0.95776315   57               Down 0.01703668
4                   Up 0.38753920  149                 Up 0.01597788
  direction_support down_expected up_expected   down_FDR     up_FDR
1             FALSE         FALSE        TRUE 0.16078873 0.19157530
2             FALSE         FALSE        TRUE 0.04543114 0.10453540
3              TRUE          TRUE        TRUE 0.95776315 0.04543114
4             FALSE          TRUE       FALSE 0.44290194 0.04543114

## Cross-cohort human meta-validation

          tissue n_mapped spearman_rho      cor_P down_N down_human_direction
1    Hippocampus       94    0.1134342 0.27633530     37                 Down
2 Frontal_cortex      215    0.1238125 0.07001071     66                 Down
      down_P up_N up_human_direction       up_P all_expected   down_FDR
1 0.26146562   57               Down 0.04590437        FALSE 0.26146562
2 0.01800406  149               Down 0.22469616        FALSE 0.07201622
      up_FDR
1 0.09180873
2 0.26146562
