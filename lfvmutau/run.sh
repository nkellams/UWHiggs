#!/bin/bash

set -o nounset
set -o errexit

export jobid=MiniAODSIMv2-Spring15-25ns_LFV_October13

rake analyzeMuTauMC
rake analyzeMuTauData

