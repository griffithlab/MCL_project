#!/bin/bash

#This shell script will bsub the mosdepth_caluclate.sh.
bsub -oo /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/run_scripts/logs/mosdepth_average.log -M 650000000 -R 'select[mem>65000] span[hosts=1] rusage[mem=65000]' -q research-hpc -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:20)' /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/run_scripts/mosdepth_average.sh 
