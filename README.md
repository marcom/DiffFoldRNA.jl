# DiffFoldRNA.jl

[![Build Status](https://github.com/marcom/DiffFoldRNA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcom/DiffFoldRNA.jl/actions/workflows/CI.yml?query=branch%3Amain)

Differentiable RNA partition function calculation for sequence
probability distributions and sequence design.  This implementation is
currently CPU-only and can only deal with smallish sequences in a
reasonable time.

There is also a [jax implementation](https://github.com/rkruegs123/jax-rnafold)
that can run on a GPU and can deal with somewhat longer sequences.

## Installation

```julia
using Pkg
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

Matthies, Krueger, Torda, Ward. Differentiable Partition Function
Calculation for RNA, bioRxiv, 2023.
https://www.biorxiv.org/content/10.1101/2023.01.30.526001v2.abstract
