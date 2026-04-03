# DEEP-seq: scalable, deterministic encapsulation and transcriptomic enrichment for rare cell profiling
# DEEP-seq description
Single-cell transcriptomics has revolutionized our understanding of rare cell populations, uncovering profound heterogeneity within these subsets. However, previous studies on rare cells have been constrained by platform limitations, often focusing on the static properties of the population itself while failing to track dynamic physiological changes or resolve associated lineage trajectories. To overcome this challenge, we introduce DEEP-seq (Deterministic Encapsulation and Enrichment for Profiling), a platform that integrates real-time fluorescence sorting with active microvalve logic to achieve high-throughput, deterministic cell-bead pairing. DEEP-seq demonstrates a superior recovery rate of >86% and high throughput (>28,000 cells/hour) across a broad dynamic range, while actively excluding non-viable cells. 

<img width="2971" height="835" alt="图片1" src="https://github.com/user-attachments/assets/80803105-3919-438c-a2cc-e6413f9b47c1" />

## Step 1: Raw data processing and alignment
Use STARsolo with default parameters to process the raw FASTQ data, performing cell barcode identification and UMI counting.

```
sh 0_Preprocess/run_star.sh
```

## Step 2: Downstream analysis pipeline
The downstream analysis scripts from the DEEP-seq paper have been uploaded to the code directories organized by figures.



