#!/bin/bash

wget -O liftOver http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O - | gunzip > hg19ToHg38.over.chain
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz -O - | gunzip > hg18ToHg38.over.chain
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg38.over.chain.gz -O - | gunzip > hg17ToHg38.over.chain

