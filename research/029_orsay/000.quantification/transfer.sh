rsync -av --bwlimit=40000 -e "ssh -i /Users/adrian/.ssh/elja" adrian@elja.hi.is:/users/home/adrian/research/029_orsay .

rsync -av --bwlimit=40000 -e "ssh -i /Users/adrian/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/adrian/research/029_orsay/results .
