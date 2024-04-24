#!/bin/bash

path='/40t_1/caix/Bol_analysis/re-annotations/augustus/train/JZS'

perl   ${path}/geneModels2AugusutsTrainingInput  ${path}/BOL.seq.20110802.chr_check.gff.added.new_20130704  ${path}/BOL.seq.lst.new.chr20110802_check.fa.added_20130704  --cpu 50
