# Gene-level cross-cohort meta-analysis

Fixed-effect inverse-variance meta-analysis was performed separately for hippocampus and frontal cortex. Surgery date was included in each cohort-specific model.

## Signal counts

          tissue           contrast meta_FDR05 concordant_FDR05
1 Frontal_cortex      CD33_specific        302              302
2    Hippocampus      CD33_specific        242              233
3 Frontal_cortex Control_ASO_effect        875              859
4    Hippocampus Control_ASO_effect       1451             1413
5 Frontal_cortex            Disease       2714             2693
6    Hippocampus            Disease        459              459
  low_heterogeneity_FDR05
1                     245
2                     108
3                     722
4                     915
5                    2333
6                     346

## CD33 meta-analysis

                                          tissue      contrast logFC_75ug
Hippocampus Disease.1747             Hippocampus       Disease  0.3424101
Hippocampus CD33_specific.1747       Hippocampus CD33_specific -1.5529770
Frontal_cortex Disease.1747       Frontal_cortex       Disease  0.7044291
Frontal_cortex CD33_specific.1747 Frontal_cortex CD33_specific -1.0888253
                                  logFC_50ug meta_logFC       meta_P
Hippocampus Disease.1747           0.4150842  0.3828293 1.246325e-07
Hippocampus CD33_specific.1747    -0.7431527 -1.0806524 5.822432e-56
Frontal_cortex Disease.1747        0.7780759  0.7451154 8.998765e-17
Frontal_cortex CD33_specific.1747 -0.5633474 -0.7962347 3.346662e-21
                                      meta_FDR       I2
Hippocampus Disease.1747          1.084745e-05  0.00000
Hippocampus CD33_specific.1747    7.145288e-52 97.05068
Frontal_cortex Disease.1747       6.001785e-15  0.00000
Frontal_cortex CD33_specific.1747 2.053512e-17 89.58418
