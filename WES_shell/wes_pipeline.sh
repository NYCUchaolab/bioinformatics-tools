#!/bin/bash
source "../config.sh"

python "${tools_dir}/others/send_msg.py" "start bqsr"
./run_bqsr.sh
python "${tools_dir}/others/send_msg.py" "start muse"
./run_muse.sh
python "${tools_dir}/others/send_msg.py" "start varscan"
./run_varscan.sh 
python "${tools_dir}/others/send_msg.py" "start mutect"
./run_mutect.sh
python "${tools_dir}/others/send_msg.py" "start pindel"
./run_pindel.sh
# ./run_vcf2maf.sh
