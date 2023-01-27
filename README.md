# DiffFoldRNA.jl

Differentiable RNA partition function calculation for sequence
probability distributions and sequence design.

## Installation

```julia
using Pkg
pkg"add https://github.com/marcom/FoldRNA.jl"
pkg"add https://github.com/marcom/DiffFoldRNA.jl"
```

## Usage

```julia
using DiffFoldRNA
```

### Sequence design

```julia
target = "(((...)))"
design_ptarget(target, ViennaModel(); optim_options=(; iterations=10))
```

### Sequence-structure partition function

```julia
seqstruct_partition(ones(10,4)/4, ViennaModel(); hpmin=3)
```

## Citation

TODO
