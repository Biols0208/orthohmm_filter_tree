python step1_deal_orthohmm_gene_count.py orthohmm_gene_count.txt orthohmm_gene_count.txt.out1

##  StdDev<1 && Species_Present_Percent>60
awk '$69<1 && $75>36' orthohmm_gene_count.txt.out1 > orthohmm_gene_count.txt.out1.deal1
awk '{print $1}' orthohmm_gene_count.txt.out1.deal1 | sed 's/://g' > orthohmm_gene_count.txt.out1.deal1.group
mkdir cds

cat orthohmm_gene_count.txt.out1.deal1.group | parallel --jobs 48 '
    grep -f <(echo {}) orthohmm_orthogroups.txt \
    | sed "s/ /\n/g" | sed "1d" \
    | seqtk subseq all.deal.cds - > cds/{}.cds
'

#for i in $(cat orthohmm_gene_count.txt.out1.deal1.group); do echo ${i} | grep -f - orthohmm_orthogroups.txt | sed 's/ /\n/g' | sed '1d' | seqtk subseq all.deal.cds - > cds/${i}.cds; done
