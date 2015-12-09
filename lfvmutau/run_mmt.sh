#!/bin/bash

set -o nounset
set -o errexit

export jobid=mmt_13TeV

rake analyzeMuMuTau_FakeRate

