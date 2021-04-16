# Code for sequence processing usrinf [USEARCH11](https://drive5.com/usearch/)


## [Demultiplex forward reads](https://www.drive5.com/usearch/manual/cmd_fastx_demux.html)
In working directory

```
mkdir forward_reads

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_demux Undetermined_S0_L001_R1_001.fastq -index Undetermined_S0_L001_I1_001.fastq -barcodes combined_barcode.fa -fastqout forward_reads/fwd_reads_demux.fq

#05:21 76Mb    100.0% Demuxed 5330416 / 18913331 (28.2%)

```

## [USEARCH phix removal](https://www.drive5.com/usearch/manual/cmd_filter_phix.html)

```
mkdir no_phix_USEARCH

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -filter_phix forward_reads/fwd_reads_demux.fq -output no_phix_USEARCH/fwd_reads_demux_phix_filtered.fq -alnout no_phix_USEARCH/fwd_reads_demux_phix_hits.txt
#00:58 882Mb   100.0% Filtering for phix, 114082 hits (2.1%)
```

## [USEARCH primer removal](https://www.drive5.com/usearch/manual/cmd_fastx_trim_primer.html)

make a fasta file with primers
```
nano forward_reads/Taylor_primers.fa

>939F  
TTGACGGGGGCCCGCACAAG
>5.8S-Fun 
AACTTTYRRCAAYGGATCWCT

mkdir USEARCHv11

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_trim_primer no_phix_USEARCH/fwd_reads_demux_phix_filtered.fq -db forward_reads/Taylor_primers.fa -strand both  -maxdiffs 5 -width 8 -fastqout USEARCHv11/fwd_reads_demux_phix_filtered_trimmed_reads.fq -tabbedout USEARCHv11/fwd_reads_demux_phix_filtered_trimmed_reads_primer_remove.txt
```



## [Quality check of the sequences](https://drive5.com/usearch/manual/cmd_fastq_eestats2.html)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 no_phix_USEARCH/fwd_reads_demux_phix_filtered.fq -output USEARCHv11/fwd_reads_demux_phix_filtered_eestats2.txt

```
```
5216334 reads, max len 301, avg 301.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4785001( 91.7%)    5191827( 99.5%)    5216290(100.0%)
   100    3655270( 70.1%)    4250100( 81.5%)    4818059( 92.4%)
   150    2975285( 57.0%)    3630645( 69.6%)    4268179( 81.8%)
   200    2051749( 39.3%)    2846562( 54.6%)    3652789( 70.0%)
   250     925939( 17.8%)    1709749( 32.8%)    2724598( 52.2%)
   300     207108(  4.0%)     592435( 11.4%)    1416616( 27.2%)
```
   
## Filtering and truncating the merged sequences by MaxEE and set length using [fastq_filter](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter no_phix_USEARCH/fwd_reads_demux_phix_filtered.fq -fastq_maxee 1 -fastq_trunclen 150 -fastaout USEARCHv11/fwd_reads_demux_phix_filtered_fil.fa
```
```
00:31 628Mb   100.0% Filtering, 69.6% passed
   5216334  Reads (5.2M)
         0  Discarded reads length < 150
   1585689  Discarded reads with expected errs > 1.00
   3630645  Filtered reads (3.6M, 69.6%)

```

## Filtering to unique sequences with [fastx_uniques](https://drive5.com/usearch/manual/cmd_fastx_uniques.html)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques USEARCHv11/fwd_reads_demux_phix_filtered_fil.fa  -fastaout USEARCHv11/uniques_fwd_reads_demux_phix_filtered_fil.fa -sizeout

#00:09 2.0Gb  3630645 seqs, 648384 uniques, 519937 singletons (80.2%)
#00:09 2.0Gb  Min size 1, median 1, max 85030, avg 5.60
#00:14 2.0Gb   100.0% Writing USEARCHv11/uniques_fwd_reads_demux_phix_filtered_f                                                                              il.fa
```



## Clustering into 0.97 OTUs using UPARSE and [cluster_otus](https://drive5.com/usearch/manual/cmd_cluster_otus.html)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_otus USEARCHv11/uniques_fwd_reads_demux_phix_filtered_fil.fa -otus USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa -uparseout USEARCHv11/fwd_reads_demux_phix_filtered_fil_otus_uparse.txt -relabel OTU

#00:15 65Mb    100.0% 4833 OTUs, 1135 chimeras
```

## Mapping reads back to OTUs at a 97% similarity using [otutab](http://www.drive5.com/usearch/manual/cmd_otutab.html)

```
cd /mnt/home/belldere/MUD_SequenceData/

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab forward_reads/fwd_reads_demux.fq -otus USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa -mapout USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU_map.txt -otutabout USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU_table.txt -biomout USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU_jsn.biom -notmatchedfq USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU_unmapped.fq -dbmatched USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU_with_sizes.fa -sizeout


#02:41:08 55Mb    100.0% Searching fwd_reads_demux.fq, 79.5% matched
#4238866 / 5330416 mapped to OTUs (79.5%)
#02:41:08 55Mb   Writing forward_reads/fwd_reads_demux_nophix_cut_fil_OTU_table.txt
#02:41:08 55Mb   Writing forward_reads/fwd_reads_demux_nophix_cut_fil_OTU_table.txt ...done.
#02:41:08 55Mb   Writing forward_reads/fwd_reads_demux_nophix_cut_fil_OTU_jsn.biom
#02:41:08 55Mb   Writing forward_reads/fwd_reads_demux_nophix_cut_fil_OTU_jsn.biom ...done.
```

## Classifying taxa against [SILVAv123](https://www.drive5.com/usearch/manual/sintax_downloads.html) and [UNITev8_02.02.2019](https://unite.ut.ee/repository.php) database using [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)

```
#Fungal
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa -db /mnt/research/EvansLab/Databases/Silva132_release/silva_16s_v123.fa -tabbedout USEARCHv11/16s_fwd_reads_demux_phix_filtered_fil_OTU_tax_v123.sintax -strand both 

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax USEARCHv11/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa -db /mnt/research/EvansLab/Databases/UNITE_8.0/utax_reference_dataset_all_02.02.2019.udb -tabbedout USEARCHv11/ITS_fwd_reads_demux_phix_filtered_fil_OTU_tax_v8.sintax -strand both 
```
