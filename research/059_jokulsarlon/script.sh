#!/bin/bash
    
#SBATCH --job-name=transfer
#SBATCH --partition=mimir
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=2
#SBATCH --hint=nomultithread  
#SBATCH --output=messages/messages.out.txt   
#SBATCH --error=messages/messages.err.txt 

date
pwd

#
# sync files
#
echo ""
echo "about to sync files"
date

# lets run a test. make sure you have appropriate destination path
#rsync -azPhv /hpcdata/Mimir/adrian/door2next/dav/28-89.ndpi /proj/hpcdeeplearningtissue/hrb44/data/dav/.

# this is the full transfer
#rsync -azPhv /hpcdata/Mimir/adrian/door2next/dav/* /proj/hpcdeeplearningtissue/hrb44/data/dav/.

date

#
# last message
#
echo "all done."
date