#!/bin/bash

READS=$1
REF=$2
MAPPING=$3

minimap2 -x map-hifi -t12 -a $REF $READS | samtools view -b | samtools sort > $MAPPING
samtools index $MAPPING
