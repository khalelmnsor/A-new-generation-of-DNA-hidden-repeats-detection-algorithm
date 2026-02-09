# A New Generation of DNA Hidden Repeats Detection Algorithm  
### Segmentation-Based Statistical and Consensus-Driven Analysis of Genomic Periodicity

## Overview

This repository presents a **segmentation-based algorithm for detecting hidden and imperfect DNA repeats** in genomic sequences.  
Unlike classical tandem repeat detectors that rely on exact matching or fixed periodicity, this method is designed to:

- Detect **weak, mutated, and hidden repeats**
- Operate at **genome scale with bounded memory**
- Combine **consensus-based similarity** and **statistical significance testing**
- Validate detected patterns against **composition-preserving null models**

The algorithm was developed as part of a capstone research project and is intended for large-scale genomic analysis, including applications to **isochore research and genome organization**.

---

## Key Contributions

- Segment-based genome processing enabling scalable analysis
- Mismatch-tolerant representative word detection
- Dual detection pipelines:
  - AVG-based similarity optimization
  - P-value–based statistical optimization
- Column-wise binomial significance testing
- Null model validation using locally composition-preserving random DNA
- Extensive quantitative evaluation and visualization framework

---

## Algorithm Pipeline

### 1. Input

- DNA sequence `S`
- Fixed segment length `L`
- Word length range `[K_min, K_max]`
- Maximum allowed mismatches `M`
- Statistical significance threshold `α`
- Null model window size `W`

---

### 2. Sequence Segmentation
S = S₁ | S₂ | S₃ | ... | Sₙ

Each segment is analyzed independently, allowing genome-scale processing with constant memory usage.

---

### 3. Candidate Word Generation

For each segment and for each candidate word length `K`:

- The segment is partitioned into non-overlapping `K`-mers
- Candidate words are extracted **only from observed k-mers**
- This avoids exponential enumeration and ensures biological relevance

---

### 4. Representative Word Selection (Mismatch-Tolerant)

For each candidate word `c`, its **support** is computed:
 support(c) = number of k-mers with Hamming distance ≤ M from c
The representative word is selected as: w* = argmax_c support(c)
v
This allows detection of repeats even under substitutions and mutations.

---

### 5. Similarity Scoring (AVG Pipeline)

For a given representative word `w*`, the **average positional match score** is computed:

- Each segment k-mer is aligned to `w*`
- Optional cyclic rotations are allowed
- AVG score is defined as the percentage of matching positions

The best `K` is selected by maximizing the AVG score.

---

### 6. Statistical Significance Testing (P-value Pipeline)

For each position (column) in the segment:

- Nucleotide counts `{A, C, G, T}` are computed
- Let:
  - `n` = number of k-mers
  - `k` = count of the most frequent nucleotide

The binomial tail probability is computed:


The DNA sequence is divided into **non-overlapping segments of fixed length `L`**:
 P(X ≥ k) = Σ_{i=k}^n (n choose i) (1/4)^i (3/4)^(n-i)

**Note:**  
No multiplicative factor of 4 is applied, since the most frequent nucleotide is already observed.

Column P-values are combined using **Fisher’s method**: χ² = -2 Σ log(P_i)

The best `K` is selected by minimizing the combined P-value.

---

### 7. Segment Annotation

Each segment stores:

- Representative word
- Representative word count
- Two most frequent k-mers
- AVG score
- Combined P-value
- Selected `K`
- Segment classification (`strong`, `weak`, `noise`)

---

### 8. Null Model Generation

To validate statistical significance, a **composition-preserving random DNA sequence** is generated:

- Local nucleotide frequencies estimated using a sliding window `W`
- Pseudocount smoothing applied
- Random sequence generated position-by-position

The **same pipeline** is applied to the null DNA.

---

### 9. Strength Scoring via Real vs Null Comparison

For each segment: strengthScore = P_null / P_real

Classification:
- `strong` : ratio ≥ 10
- `weak`   : 1 < ratio < 10
- `noise`  : ratio ≤ 1

---

## Experimental Evaluation

The repository includes a full Python analysis framework that produces:

- CDFs of statistical significance
- Histograms of −log10(P-values)
- AVG score distributions
- Representative word length distributions
- Representative vs frequent word agreement
- Significance threshold statistics

All plots are generated for:
- AVG vs P-value pipelines
- Real vs null sequences

---

## Repository Structure

.
├── src/ # C++ implementation
│ ├── segmentation.cpp
│ ├── representative.cpp
│ ├── statistics.cpp
│ └── pipeline.cpp
│
├── analysis/ # Python analysis scripts
│ ├── analyze_results.py
│
├── data/ # Input DNA and CSV outputs
│
├── figures/ # Generated plots
│
├── README.md

---

## Why This Algorithm Matters

- Detects hidden periodicity missed by classical repeat finders
- Robust to mutation noise
- Statistically grounded
- Scalable to full genomes
- Suitable for isochore-level genomic analysis

---

## Future Directions

- Integration with GC-content segmentation
- Adaptive segment lengths
- Extension to amino acid sequences
- Genome-wide isochore-repeat correlation analysis

---

## Authors

**Khalil Mansour || Email: Khalel.Mnsor@e.braude.ac.il**  
**Fatmeh Zoabi || Email: FatmehZo3bi10@gmail.com**  
Supervised by **Dr. Zakharia Frenkel**

---

## License

This project is intended for academic and research use.  
Please cite appropriately if used in publications.
