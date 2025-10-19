#!/bin/bash

# Quick test that everything is runnable

set -o errexit
set -o nounset

# cd to dir where this script is
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# unit tests
julia --project=. -e 'import Pkg; Pkg.test()'

# validation scripts
julia --project=validation validation/run-all-scripts.jl

# benchmark scripts
julia --project=benchmark benchmark/old-bench-seqstruct-partition.jl
