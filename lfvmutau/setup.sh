#!/bin/bash

# Get the data
export datasrc=/hdfs/store/user/nkellams/
export jobid=MiniAODSIMv2-Spring15-25ns_LFV_October13
export jobid8TeV=$jobid
export afile=`find $datasrc/$jobid | grep root | head -n 1`

## Build the cython wrappers
rake "make_wrapper[$afile, mt/final/Ntuple, MuTauTree]"

ls *pyx | sed "s|pyx|so|" | xargs rake 

rake "meta:getinputs[$jobid, $datasrc,mt/metaInfo, mt/eventCount]"
rake "meta:getmeta[inputs/$jobid, mt/metaInfo, 13, mt/metaInfo/summedWeights]"


