#!/bin/bash

dos2unix *.txt
dos2unix *.tsv

#scRNA-seq output file filtered to have only lines for Hindgut (1) cluster to get sp48_hindgut_genes.tsv, from which only gene WHL IDs were kept.
tail -n +2 sp48_hindgut_genes.tsv | cut -f 4 | sort | uniq > sp48_hindgut_cells_whls.txt

#Using the footprinting output and filtering it to have only human TF name and putative CRM target.
cut -f 4,10 Sp48hpf_gut_kept_TFs.txt > Sp_48hpf_gut_hTF_pCRM.txt

#Filtering the DE expression results to only have WHL IDs and fold change, necessary for sign assignment script.
cut -f 1,3 hindgut_Sp48_lox_wt_gc_seq_bias_oneLox.tsv | tail -n +2 | sort | sed '1s/^/target_gene	fold_change\n/' > sorted_SpPdx1_affected_WHLs_foldchange_Sp48hpf.txt

#ATAC-seq concordant peakset annotated with genes.
annotatePeaks.pl Sp_putative_CRMs_48hpf_fixed.bed ~/Data/Spur_3.1.LinearScaffold.fa -gtf ~/Data/Transcriptome.gtf -size given > peakgenes_Sp_putative_CRMs_48hpf_fixed.bed

#Adding target genes to human TF and their target sea urchin pCRM.
paste <(cat Sp_48hpf_gut_hTF_pCRM.txt) <(cut -f 2 Sp_48hpf_gut_hTF_pCRM.txt | while read line; do grep -w $line peakgenes_Sp_putative_CRMs_fixed.txt || echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"; done | cut -f 11 | cut -d "." -f 1,2) > Sp_48hpf_gut_hTF_pCRM_tWHL.txt

#Adding sea urchin TF WHL IDs to the file (adding which human homolog is sea urchin homolog).
paste <(cut -f 1 Sp_48hpf_gut_hTF_pCRM_tWHL.txt | while read line; do grep -w ^$line TF_WHL.txt || echo -e "$line\t$line"; done) <(cat Sp_48hpf_gut_hTF_pCRM_tWHL.txt) > Sp_48hpf_gut_hTF_tfWHL_hTF_pCRM_gWHL.txt

cut -f 4 Sp_putative_CRMs_48hpf_fixed.bed | grep -w -f - Sp_48hpf_gut_hTF_tfWHL_hTF_pCRM_gWHL.txt | sort -V -k 4,4 | uniq > Sp_48hpf_gut_hTF_tfWHL_hTF_pCRM_gWHL_2dpf_gut_pCRMs.txt

#Keeping only info from sea urchin, removing human TF information
cut -f 2,5 Sp_48hpf_gut_hTF_tfWHL_hTF_pCRM_gWHL_2dpf_gut_pCRMs.txt > Sp_48hpf_gut_tfWHL_gWHL_2dpf_gut_pCRMs.txt

#Filtering to keep only Hindgut (1) effector genes and their targets.
awk 'NR==FNR{a[$0];next}$1 in a && $2 in a' sp48_hindgut_cells_whls.txt Sp_48hpf_gut_tfWHL_gWHL_2dpf_gut_pCRMs.txt > Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_2dpf_gut_pCRMs.txt

#Filtering to keep only TF factor as both effector and target genes.
awk 'NR==FNR{a[$0];next}$1 in a && $2 in a' tf_WHL_ids_hand_checked.txt Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_2dpf_gut_pCRMs.txt > Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt

sort Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt | uniq > sorted_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt

#Convert WHL IDs to names.
paste <(cut -f 1 sorted_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt | while read line; do grep -w -m1 $line names.txt || echo -e "$line\t$line"; done | cut -f 2) <(cut -f 2 sorted_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt| while read line; do grep -w -m1 $line names.txt || echo -e "$line\t$line"; done | cut -f 2) > named_sorted_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs_for_GRN.txt

#Filtering to keep only genes co-expressed with Sp-Pdx1.
awk 'NR==FNR{a[$0];next}$1 in a && $2 in a' pdx_positive_whls.txt Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt > pdx_pos_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt

sort pdx_pos_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt | uniq > sorted_pdx_pos_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt

#Convert WHL IDs to names.
paste <(cut -f 1 sorted_pdx_pos_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt | while read line; do grep -w -m1 $line names.txt || echo -e "$line\t$line"; done | cut -f 2) <(cut -f 2 sorted_pdx_pos_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs.txt| while read line; do grep -w -m1 $line names.txt || echo -e "$line\t$line"; done | cut -f 2) > named_sorted_pdx_pos_Sp_48hpf_gut_tfWHL_gWHL_hindgut_cells_TFs_only_2dpf_gut_pCRMs_for_GRN.txt

#Get positions of TF footprints in hindgut cells for visualizations.
paste Sp48hpf_gut_kept_TFs.txt Sp_48hpf_gut_hTF_tfWHL_hTF_pCRM_gWHL.txt | cut -f 1-3,13 | sort -k1,1 -k2,2n | bedtools intersect -wa -a - -b ../Sp_putative_CRMs_48hpf.bed | sort -k1,1 -k2,2n | uniq | grep -w -f ../sp48_hindgut_cells_whls.txt - | sort -k1,1 -k2,2n | uniq > Sp48_hindgut_gut_peaks_tfWHL.bed

paste <(cut -f 1-3 Sp48_hindgut_gut_peaks_tfWHL.bed) <(cut -f 4 Sp48_hindgut_gut_peaks_tfWHL.bed | while read line; do grep -w -m1 $line ../names.txt || echo -e "$line\t$line"; done | cut -f 2) | sort -k1,1 -k2,2n | uniq > named_Sp48_hindgut_gut_peaks_tfWHL.bed
