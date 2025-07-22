# orthohmm_filter_tree
Process the results output by [orthohmm](https://github.com/JLSteenwyk/orthohmm) software.

Filter orthogroups and then infer the phylogenetic tree.

Simple analysis process, configure the required environment by yourself.

----

#### step1_prepare.sh 
- Get less variation in the number of homologous genes and more species with homologous genes.
- Mainly use deal_orthohmm_gene_count.py script to calculate orthohmm_gene_count.txt and then extract cds sequence.

----

#### step2_make_all_command.sh
- Make all command to get best unique homologous gene.
- Improve cds sequences whose length is not a multiple of three.
- Of course you can choose not to filter.

----

#### step3_OMM_MACSE.sh
- [OMM_MACSE](https://www.agap-ge2pop.org/macse/pipeline-documentation/): Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons
