

##### OVERLAP BETWEEN SC AND BULK ATAC PEAK

module load bedtools/2.31.0

bedtools intersect -a data/compare_sc_bulk_atac/atacpeaks_DE_hg38.bed -b data/compare_sc_bulk_atac/sc_telohaec_tnf_atac_peaks.bed -wo > data/compare_sc_bulk_atac/atac_peak_overlap_amount.bed