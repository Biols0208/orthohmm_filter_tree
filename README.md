# orthohmm_filter_tree
Process the results output by orthohmm software.

Filter orthogroups and then infer the phylogenetic tree.

Simple analysis process, configure the required environment by yourself.

----

#### step1_prepare.sh 
- Get less variation in the number of homologous genes and more species with homologous genes.
- Mainly use step1_deal_orthohmm_gene_count.py script to calculate orthohmm_gene_count.txt and then extract cds sequence.

