#!/bin/bash

#rsync -av --bwlimit=40000 -e "ssh -i //Users/adrian/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/adrian/entrance2nextcloud/021_arnarholt_data/Sheepseq_sample1 /Users/adrian/Nextcloud/door2necio5/021/.
rsync -av --bwlimit=40000 -e "ssh -i /Users/adrian/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/adrian/entrance2nextcloud/021_arnarholt_data/Sheepseq_sample2 /Users/adrian/Nextcloud/door2necio5/021/.
