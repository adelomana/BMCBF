#!/bin/bash

rsync -azPhv --bwlimit=40000 -e "ssh -i /Users/adrian/.ssh/elja" /Users/adrian/scratch/lilit_2026.03.10 adrian@elja.hi.is:/hpcdata/Mimir/adrian/research/071_ushguli/data/.