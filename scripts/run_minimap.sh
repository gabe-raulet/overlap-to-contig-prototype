#!/bin/bash

if (($# != 3))
then
    echo "usage: $0 <reads.fa> <reference.fa> <mapping.bam>"
    exit
fi

READS=$1
REF=$2
MAPPING=$3

minimap2 -x map-hifi -t12 -a $REF $READS | samtools view -b | samtools sort > $MAPPING
samtools index $MAPPING
