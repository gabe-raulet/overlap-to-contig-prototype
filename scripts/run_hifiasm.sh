#!/bin/bash

INPUT=$1
OUTPUT=$2

hifiasm -o _tmp.asm -t12 $INPUT
gfatools gfa2fa _tmp.asm.bp.p_utg.gfa > $OUTPUT
rm _tmp.asm.{bp,ec,ovlp}.*
