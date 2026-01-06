# A New Generation of DNA Hidden Repeats Detection Algorithm  
## With Applications to Isochore Research

---

## 📌 Overview
This repository presents the implementation of a **novel segmentation-based algorithm** for detecting **hidden DNA repeats** using **dominant k-mer patterns** and **statistical significance testing**.

The method is designed to identify **latent repetitive structures** that are not detectable by classical repeat-finding algorithms, and to explore their relationship with **GC-content variability** and **isochore organization**.

This project was developed as part of an **academic research / capstone project** in the field of **Computational Biology and Bioinformatics**.

---

## 🧬 Motivation
Classical repeat detection algorithms primarily focus on:
- Exact or approximate tandem repeats  
- Explicit sequence similarity  

However, many biologically meaningful repeats are:
- **Statistical rather than explicit**
- **Distributed across genomic segments**
- **Hidden within heterogeneous genomic regions**

This work addresses these limitations by introducing a **dominant k-mer–based segmentation framework** combined with **rigorous statistical testing**.

---

## 🧠 Core Ideas
- Segment the DNA sequence into fixed-length windows  
- Identify **dominant k-mers** at each positional index  
- Compute **position-specific statistical significance**  
- Merge adjacent statistically similar segments  
- Analyze correlations with **GC-content** and **isochore structures**

---

## ⚙️ Algorithm Summary

### 1. Segmentation
- Divide the DNA sequence into segments of length `L`
- Each segment is decomposed into **non-overlapping k-mers**

### 2. Representative Word Construction
- For each segment, compute a **representative k-mer**
- The dominant nucleotide is selected at each k-mer position

### 3. Statistical Testing
- Compute **position-wise p-values** using a **binomial model**
- Combine p-values using **Fisher’s method**
- Assign statistical significance to each segment

### 4. Segment Merging
- Adjacent segments with compatible dominant patterns and significance levels are merged
- This step reduces noise and over-segmentation

### 5. Biological Interpretation
- Compare detected segments with **GC-content profiles**
- Investigate alignment with known or inferred **isochore regions**

---

## 📊 Key Features
- Detection of **hidden (implicit) DNA repeats**
- Robust against noise and local mutations
- Strong statistical foundation (binomial tests + Fisher’s method)
- Suitable for **large genomic sequences**
- Designed for **isochore-level analysis**

---

## 🛠️ Implementation Details
- **Language:** C++  
- **Paradigm:** Modular, algorithm-oriented design  
- **Input:** DNA sequences (FASTA or raw string)  
- **Output:** Segmentation results, dominant patterns, statistical scores  

---

## 📁 Project Structure

---

## 🎓 Academic Context
This project builds upon and extends ideas from:
- Statistical sequence segmentation
- Dominant word representations
- Hidden repeat detection
- Isochore theory in genome organization

It is intended for:
- Bioinformatics researchers
- Computational biology students
- Genomic sequence analysis studies

---

## 🚀 Future Work
- Support for overlapping k-mers  
- Adaptive segment length selection  
- Integration with genome browsers  
- Validation on full eukaryotic genomes  
- Performance optimization and parallelization  

---

## 👤 Authors
**Khalil Mansour**  
📧 Khalel.Mnsor@e.braude.ac.il  

**Fatmeh Zoabi**  
📧 fatmehzo3bi10@gmail.com  

---

## 📜 License
This project is released under the **MIT License**, allowing academic and research use with attribution.
