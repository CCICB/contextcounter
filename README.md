# contextcounter

> \[!WARNING\]  
> This package is in early development and not ready for use

**`contextcounter`** is a command-line tool for counting the opportunities to mutate different sequence contexts offered by genome, exome, or gene panel reference sequences.

These context counts are essential for **renormalizing** mutational signatures derived from whole-genome data, enabling their accurate application to sequencing data from exomes or custom panels. This functionality is used by tools including the [`sigstats`](https://github.com/CCICB/sigstats) R package.

---

## 🔧 Features

- Counts **di-**, **tri-**, and **pentanucleotide** contexts
- Skips user-specified contigs (so you can exclude mitochondrial or sex chromosomes)
- Outputs context count tables per type
- Streamless integration with downstream signature tools (sigverse)

---

## 🚀 Installation

Install directly from GitHub using Cargo:

```bash
cargo install --git https://github.com/CCICB/contextcounter
```


## Quick Start

```
contextcounter genome.fasta
```