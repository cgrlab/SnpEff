#!/bin/sh

# VCF="1kg/ALL.wgs.phase1.projectConsensus.snps.sites.vcf"
VCF="$HOME/snpEff/db/GRCh37/1kg/1kg.vcf"

ANN=`dirname $VCF`/`basename $VCF .vcf`.ann.vcf
REF=GRCh37.75

SNPEFF="java -Xmx8G -jar snpEff.jar"
SNPSIFT="java -Xmx4G -jar SnpSift.jar"

#---
# Annotations
#---

time $SNPEFF ann -v -stats $VCF.html $REF $VCF > $ANN

echo Annotate using dbSnp, dbNSFP and ClinVar
time $SNPSIFT annotate -v $ANN
java -Xmx1g -jar SnpSift.jar \
    annotate \
    -v \
    protocols/db/clinvar_00-latest.vcf \
    protocols/ex1.eff.cc.vcf \



