#!/bin/bash

rsync -av --bwlimit=40000 -e "ssh -i /Users/adrian/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/shared/snaevar/Hildur/F24A430000789_CANfgyzR/soapnuke/clean/HRH88* /Users/adrian/research/bmcbf/011.askja/doc/geo/.
