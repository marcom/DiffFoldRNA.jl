#!/bin/bash

set -o errexit
set -o nounset

# cd to dir where this script is
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# unit tests
julia --project=. -e 'import Pkg; Pkg.test()'

# validation scripts
julia --project=validation check-seqstruct-partition-All1Model-vienna-vs-nussinov.jl
