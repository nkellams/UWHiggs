#!/bin/bash

# Get the data
export datasrc=/hdfs/store/user/nkellams/
export jobid=mmt_13TeV
export jobid8TeV=$jobid
export afile=`find $datasrc/$jobid | grep root | head -n 1`

## Build the cython wrappers
rake "make_wrapper[$afile, mmt/final/Ntuple, MuMuTauTree]"

ls *pyx | sed "s|pyx|so|" | xargs rake 

rake "meta:getinputs[$jobid, $datasrc,mmt/metaInfo, mmt/eventCount]"
rake "meta:getmeta[inputs/$jobid, mmt/metaInfo, 13, mmt/metaInfo/summedWeights]"


