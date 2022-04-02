#!/bin/bash

if (($# != 2))
then
    echo "usage: $0 <reads.fa> <assembly.fa>"
    exit
fi

READS=$1
ASM=$2

hifiasm -o _tmp.asm -t12 $READS
gfatools gfa2fa _tmp.asm.bp.p_utg.gfa | seqtk seq -A -l0 > $ASM
rm _tmp.asm.{bp,ec,ovlp}.*
