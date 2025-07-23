# orthohmm_filter_tree
Process the results output by [orthohmm](https://github.com/JLSteenwyk/orthohmm) software.

Filter orthogroups and then infer the phylogenetic tree.

Simple analysis process, configure the required environment by yourself.

----

#### step1_prepare.sh 
- Get less variation in the number of homologous genes and more species with homologous genes.
- Mainly use **deal_orthohmm_gene_count.py** script to calculate orthohmm_gene_count.txt and then extract cds sequence.

----

#### step2_make_all_command.sh
- Make all command to get best unique homologous gene.
- Improve cds sequences whose length is not a multiple of three.
- Of course you can choose not to filter.

----

#### step3_OMM_MACSE.sh
- [OMM_MACSE](https://www.agap-ge2pop.org/macse/pipeline-documentation/): Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons

----

#### step4_FastTree.sh
- Constructing gene trees.

----

Finally, inferring species trees from gene trees by [astral-pro3](https://github.com/chaoszhang/ASTER) or [tree-qmc](https://github.com/molloy-lab/TREE-QMC).

Of course, you can use 4d site sequences for analysis, or connect all multiple sequence alignment files into supergene for analysis (merge_phy.py).

example:
```
tree out/OG06164
|-- OG06164.best
|-- OG06164.best.fa
|-- OG06164_best_sequences.fasta
|-- OG06164_config_backup.conf
|-- OG06164_file_manifest.txt
|-- OG06164_pipeline_report.txt
|-- OG06164_processing.log
|-- OG06164_quality_best_sequences.fasta
|-- OG06164_quality_removed_sequences.fasta
|-- OG06164_removed_sequences.fasta
`-- omm_macse
    |-- OG06164_final_align_AA.aln
    |-- OG06164_final_align_AA.aln.phy
    |-- OG06164_final_align_NT.aln
    |-- OG06164_final_align_NT.aln.phy
    |-- OG06164_final_align_NT.aln.phy.gaps.4d.phy
    |-- OG06164_final_align_NT.aln.phy.nogaps.4d
    |-- OG06164_final_align_NT.aln.phy.nogaps.4d.phy
    |-- OG06164_final_mask_align_AA.aln
    |-- OG06164_final_mask_align_NT.aln
    |-- OG06164_final_unmask_align_AA.aln
    |-- OG06164_final_unmask_align_NT.aln
    |-- OG06164_final_unmask_align_NT.aln.phy.nogaps.4d
    |-- OG06164_maskFull_detail.fasta
    |-- OG06164_maskFull_stat.csv
    |-- OG06164_maskHomolog_stat.csv
    |-- all.phy.gz
    |-- all.tree.gz
    |-- all.tree.nwk
    `-- readme_output.txt
```
