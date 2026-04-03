#!/bin/bash

/usr/bin/time -v PreProcessR -i /home/Users/rl152/Strainify/4-strain-mock-exp/genomes 2> build.log
/usr/bin/time -v StrainR -1 /home/Users/rl152/Strainify/4-strain-mock-exp/fastq/SRR13355226_r1.fq -2 /home/Users/rl152/Strainify/4-strain-mock-exp/fastq/SRR13355226_r2.fq -r StrainR2DB/ 2> SRR13355226.log