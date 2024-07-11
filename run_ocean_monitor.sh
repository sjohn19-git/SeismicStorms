#!/bin/bash
#source /opt/conda/mini/etc/profile.d/conda.sh
source ~/.bashrc
conda activate data-vis
python3 /home/sjohn/projects/data-visuals/SeismicStorms/SeismicStorms.py 2> /home/sjohn/projects/data-visuals/SeismicStorms/error_log.txt
python3 /home/sjohn/projects/data-visuals/SeismicStorms/video_generator.py

