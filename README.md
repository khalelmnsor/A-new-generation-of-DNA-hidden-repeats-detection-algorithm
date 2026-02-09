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

The DNA sequence is divided into **non-overlapping segments of fixed length `L`**:

