#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <array>
#include <chrono>
#include <random>
#include <unordered_map>
#include <fstream>
#include <iomanip>

using namespace std;
using CountTable = vector<map<char, int>>;

/*
 * Global parameters
 * L defines the fixed segment length used throughout the pipeline.
 * Segment-based processing enables genome-scale analysis with bounded memory.
 */
const int L = 420; // Segment length
const int KMAX = 6;
const int KMIN = 3;
const int MISTAKE = 1;
/*
 * Segment structure
 * -----------------
 * This structure represents a single DNA segment and all associated
 * statistics computed by the detection pipelines.
 *
 * The design intentionally separates:
 *  - observed sequence data
 *  - representative (consensus) information
 *  - statistical significance
 *  - empirical validation against null models
 */
struct Segment
{
    // Raw DNA sequence of the segment
    string sequence;

    // Canonical (rotation-invariant) representative word
    string canonicalRepresentativeWord;

    // Representative word selected by the pipeline (may be non-canonical)
    string representativeWord;

    // Most frequent k-mer observed in the segment
    string frequentWord1;

    // Second most frequent k-mer (used for ambiguity and robustness analysis)
    string frequentWord2;

    // Absolute occurrence counts
    int frequentCount1 = 0;
    int frequentCount2 = 0;
    int representativeCount = 0;

    // Position-wise nucleotide count table (A,C,G,T per position)
    vector<array<int, 4>> positioncounts;

    // Combined statistical significance (Fisher-combined p-value)
    double combinedPValue;

    // Genomic location metadata
    int startIndex;
    int length;

    // Average positional match percentage with representative word
    double avg = 0.0;

    // Noise annotation (used during post-processing / merging)
    string noise = "no";

    // Empirical significance score based on real vs null comparison
    double strengthScore = 0.0;

    // Final qualitative classification of the segment
    string strength = "unknown"; // noise / weak / strong
};

/*
 * Pipeline entry points
 * ---------------------
 * Each function implements a different selection criterion:
 *  - SegmentSequence1: AVG-based representative selection
 *  - SegmentSequence2: P-value–based representative selection
 *  - SegmentSequence3: frequency-driven baseline (control)
 */
vector<Segment> SegmentSequence1(const string &dna);
vector<Segment> SegmentSequence2(const string &dna);
vector<Segment> SegmentSequence3(const string &dna);

/*
 * Representative word construction
 * --------------------------------
 * Multiple strategies are provided:
 *  - exact consensus (position-wise majority)
 *  - mismatch-aware consensus (robust to mutations)
 */
string ComputeRepresentativeWord1(const string &segment, int K);
string ComputeRepresentativeWord_Mismatch(const string &segment, int K, int M);

/*
 * Mismatch-aware matching
 * -----------------------
 * Counts how many k-mers match a candidate word within a mismatch budget M.
 * This enables detection of weak or mutated periodic signals.
 */
int CountMatchesWithMismatches(const string &segment, const string &word, int M);

/*
 * Frequency analysis
 * ------------------
 * Extracts the two most frequent k-mers in a segment.
 * Used for validation, ambiguity detection, and agreement analysis.
 */
void ComputeTwoMostFrequentWords(
    const string &segment,
    int K,
    string &word1, int &count1,
    string &word2, int &count2);

/*
 * Counts exact occurrences of a representative word.
 * Used to compare consensus strength against raw frequency.
 */
int CountRepresentativeWord(const string &segment, const string &repWord);

/*
 * AVG-based similarity metrics
 * ----------------------------
 * Measure average positional agreement between the segment
 * and the representative word.
 */
double ComputeAvgMatch1(const string &segment, const string &repWord);
double ComputeAvgMatch2(const string &segment, const string &repWord);

/*
 * Statistical modeling
 * --------------------
 * Position-wise nucleotide counts and binomial tests
 * are used to evaluate deviation from random background.
 */
vector<array<int, 4>> ComputePositionCounts(const string &segment, int k);
vector<double> ComputePValuesFromCounts(
    const vector<array<int, 4>> &table,
    int numKmers,
    int k);

/*
 * Binomial and Fisher statistics
 * ------------------------------
 * Provide statistically interpretable significance values.
 */
double BinomialPValue(int n, int k, double p);
double BinomialCoefficient(int n, int k);
double CombinePValuesFisher(const vector<double> &pValues);
double ChiSquareSurvivalEvenDF(double X, int df);

/*
 * Post-processing utilities
 * -------------------------
 * Used to merge coherent segments and reduce fragmentation.
 */
vector<Segment> MergeSameWordSegments(const vector<Segment> &segments);
vector<Segment> MergeNoiseSegments(const vector<Segment> &segments, int k);

/*
 * Canonical rotation utilities
 * ----------------------------
 * Enforce rotation invariance for periodic repeats.
 */
string CanonicalRotation(const string &w);
int CanonicalRotationOffset(const string &w);
vector<array<int, 4>> RotateCountsLeft(
    const vector<array<int, 4>> &table,
    int shift);

/*
 * I/O and formatting
 * ------------------
 * FASTA reading and human-readable formatting utilities.
 */
string ReadFNA(const string &filename);
string FormatSequenceByK(const string &seq, int K);
vector<string> SplitByLength(const string &dna);

/*
 * Null model generation
 * ---------------------
 * Generates random DNA sequences that preserve local
 * nucleotide composition, enabling fair statistical comparison.
 */
vector<array<double, 4>> ComputeLocalDistributions(
    const string &dna,
    int W,
    double pseudo);

string GenerateRandomDNA(const vector<array<double, 4>> &probs);

/*
 * Pipeline orchestration
 * ----------------------
 * Allows running identical logic on real and null data.
 */
vector<Segment> RunPipeline(const string &dna, int mode);
vector<Segment> AnnotateSegmentsWithNull1(
    vector<Segment> &realSegs,
    const string &dna,
    int W,
    int mode);

/*
 * Export utilities
 * ----------------
 * Simple CSV export for downstream statistical analysis
 * and visualization (Python / R).
 */
void ExportSegmentsCSV_Simple(
    const vector<Segment> &segs,
    const string &filename,
    const string &pipeline,
    const string &type); // "real" or "null"
void AddCounts(vector<array<int, 4>> &a, const vector<array<int, 4>> &b);

int main()
{
    // Load genomic DNA sequence from FASTA (.fna) file
    // Header lines are ignored; sequence is concatenated into a single string
    string dna = ReadFNA("C:/Users/user/Desktop/project/Data/gene.fna");

    cout << "Loaded sequence length: " << dna.length() << endl;

    // ============================================================
    // 1. SEGMENTATION + AVG PIPELINE (REAL DNA)
    // ============================================================
    // The DNA is segmented into fixed-length windows (L),
    // and each segment is analyzed using the AVG-based pipeline
    vector<Segment> segments = SegmentSequence1(dna);
    cout << "class the segment done done: " << endl;

    // ============================================================
    // 2. NULL MODEL ANNOTATION (PARALLEL CONTROL)
    // ============================================================
    // Generate a composition-preserving random DNA sequence
    // and apply the same pipeline to obtain null segments.
    // Each real segment is compared to its null counterpart
    vector<Segment> nullSegs = AnnotateSegmentsWithNull1(
        segments,
        dna,
        L / 10, // Local window size for composition preservation
        1       // Pipeline mode: AVG
    );
    cout << "randme Segment done: " << endl;

    // ============================================================
    // Optional post-processing steps (disabled here)
    // ============================================================
    // Merge consecutive segments with identical representative words
    // segments = MergeSameWordSegments(segments);

    // Merge weak or noisy segments into neighboring strong segments
    // segments = MergeNoiseSegments(segments);

    // ============================================================
    // 3. CONSOLE OUTPUT (AVG PIPELINE)
    // ============================================================
    cout << "\nDetected Hidden Repeat Segments (AVG PIPELINE):\n";

    // Print segment-by-segment comparison between real and null DNA
    for (size_t i = 0; i < segments.size() && i < nullSegs.size(); ++i)
    {
        Segment &s = segments[i];
        Segment &nulls = nullSegs[i];

        cout << "====================================================\n";
        cout << "[SEGMENT " << i << "]\n";

        // Genomic coordinates
        cout << "Location   : start=" << s.startIndex
             << ", length=" << s.length << "\n";

        // Representative word selected by the pipeline (real DNA)
        cout << "Representative word (real): "
             << s.representativeWord
             << "  count=" << s.representativeCount << "\n";

        // Most frequent exact k-mers (real DNA)
        cout << "Frequent word #1 (real): "
             << s.frequentWord1
             << "  count=" << s.frequentCount1 << "\n";

        cout << "Frequent word #2 (real): "
             << s.frequentWord2
             << "  count=" << s.frequentCount2 << "\n";

        // Representative word in the null (randomized) DNA
        cout << "Representative word (null): "
             << nulls.representativeWord
             << "  count=" << nulls.representativeCount << "\n";

        // Most frequent k-mers in the null segment
        cout << "Frequent word #1 (null): "
             << nulls.frequentWord1
             << "  count=" << nulls.frequentCount1 << "\n";

        cout << "Frequent word #2 (null): "
             << nulls.frequentWord2
             << "  count=" << nulls.frequentCount2 << "\n";

        // Statistical comparison between real and null segments
        cout << "P-values   : real=" << s.combinedPValue
             << " | null(parallel)=" << nulls.combinedPValue;

        // Empirical strength classification based on null/real ratio
        cout << "Strength   : score=" << s.strengthScore
             << " | class=" << s.strength
             << " | merge-noise=" << s.noise << "\n";

        // Average positional match percentage
        cout << "AVG match  : " << s.avg << "%\n";
        cout << "AVG null match  : " << nulls.avg << "%\n";

        // Formatted sequence view (grouped by representative word length)
        cout << "Sequence (real):\n"
             << FormatSequenceByK(
                    s.sequence,
                    s.canonicalRepresentativeWord.size())
             << "\n";

        cout << "Sequence (null):\n"
             << FormatSequenceByK(
                    nulls.sequence,
                    nulls.canonicalRepresentativeWord.size())
             << "\n";
    }

    // ============================================================
    // 4. EXPORT RESULTS FOR OFFLINE ANALYSIS (AVG PIPELINE)
    // ============================================================
    ExportSegmentsCSV_Simple(
        segments,
        "segments_real_avg.csv",
        "avg",
        "real");
    cout << "ANALYSIS11 done: " << endl;

    ExportSegmentsCSV_Simple(
        nullSegs,
        "segments_null_avg.csv",
        "avg",
        "null");
    cout << "ANALYSIS12 done: " << endl;

    // ============================================================
    // 5. P-VALUE PIPELINE (REAL + NULL)
    // ============================================================
    // Run the alternative pipeline where representative words
    // are selected by minimizing the combined statistical P-value
    vector<Segment> seg_p_real = SegmentSequence2(dna);

    vector<Segment> seg_p_null = AnnotateSegmentsWithNull1(
        seg_p_real,
        dna,
        L / 10,
        2 // Pipeline mode: P-VALUE
    );

    // Export P-value pipeline results
    ExportSegmentsCSV_Simple(
        seg_p_real,
        "segments_real_pvalue.csv",
        "pvalue",
        "real");
    cout << "ANALYSIS21 done: " << endl;

    ExportSegmentsCSV_Simple(
        seg_p_null,
        "segments_null_pvalue.csv",
        "pvalue",
        "null");
    cout << "ANALYSIS22 done: " << endl;

    // ============================================================
    // 6. VISUALIZATION
    // ============================================================
    // Launch Python-based visualization scripts
    cout << "\nRunning Python visualization...\n";
    system("py plot.py");

    return 0;
}

// ===============================
// AVG PIPELINE
// ===============================
// This pipeline selects the representative word by maximizing
// the average positional match (AVG) between the segment and
// a mismatch-tolerant consensus k-mer.
vector<Segment> SegmentSequence1(const string &dna)
{
    // Container for all detected segments
    vector<Segment> segments;

    // Sliding buffer used to accumulate a fixed-length segment
    string buffer;
    buffer.reserve(L);

    int bufferStart = 0;

    // Scan the DNA sequence character by character
    for (size_t i = 0; i < dna.size(); ++i)
    {
        // Normalize nucleotide to uppercase
        char c = toupper(dna[i]);

        // Skip non-ACGT characters (e.g., N, gaps, headers)
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            continue;

        // Append valid nucleotide to the current segment buffer
        buffer.push_back(c);

        // Process the segment when it reaches length L
        // (or a sufficiently long tail at the end of the sequence)
        if (buffer.size() == L || (i == dna.size() - 1 && buffer.size() >= L / 2))
        {
            Segment seg;
            seg.sequence = buffer;

            string w;
            string bestW;
            double bestAvg = -1.0;
            int bestK = -1;

            // Test multiple candidate periods (K values)
            // and select the one maximizing the AVG match score
            for (int k = KMIN; k <= KMAX; k++)
            {
                // Compute a mismatch-aware representative word (M = 1)
                w = ComputeRepresentativeWord_Mismatch(buffer, k, MISTAKE);

                // Measure average positional agreement with the segment
                double avg = ComputeAvgMatch1(buffer, w);

                // Keep the K and word that maximize AVG
                if (avg > bestAvg)
                {
                    bestAvg = avg;
                    bestW = w;
                    bestK = k;
                }
            }

            // Extract the two most frequent exact k-mers
            // (used for agreement and robustness analysis)
            ComputeTwoMostFrequentWords(
                buffer,
                bestK,
                seg.frequentWord1,
                seg.frequentCount1,
                seg.frequentWord2,
                seg.frequentCount2);

            // Store the selected representative word
            seg.representativeWord = bestW;

            // Count exact occurrences of the representative word
            seg.representativeCount =
                CountRepresentativeWord(buffer, seg.representativeWord);

            // Convert representative word to its canonical (rotation-invariant) form
            seg.canonicalRepresentativeWord = CanonicalRotation(bestW);

            // Store AVG score as the primary selection criterion
            seg.avg = bestAvg;

            // Align position-wise statistics with the canonical rotation
            int shift = CanonicalRotationOffset(bestW);

            int numKmers = seg.sequence.size() / bestK;

            // Compute nucleotide counts per position
            seg.positioncounts = ComputePositionCounts(seg.sequence, bestK);

            // Rotate counts to canonical orientation
            seg.positioncounts = RotateCountsLeft(seg.positioncounts, shift);

            // Compute position-wise binomial p-values
            vector<double> positionPValues =
                ComputePValuesFromCounts(seg.positioncounts, numKmers, bestK);

            // Combine positional p-values using Fisher's method
            seg.combinedPValue = CombinePValuesFisher(positionPValues);

            // Record genomic coordinates
            seg.startIndex = i - buffer.size() + 1;
            seg.length = buffer.size();

            // Store the completed segment
            segments.push_back(seg);

            // Reset buffer for the next segment
            buffer.clear();
        }
    }

    return segments;
}

// ===============================
// P-VALUE PIPELINE
// ===============================
// This pipeline selects the representative word and period (K)
// by minimizing a statistically combined P-value derived from
// position-wise nucleotide distributions.
vector<Segment> SegmentSequence2(const string &dna)
{
    vector<Segment> segments;
    string buffer;
    buffer.reserve(L);

    // Iterate through the DNA sequence
    for (size_t i = 0; i < dna.size(); ++i)
    {
        char c = toupper(dna[i]);

        // Ignore non-standard nucleotides
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            continue;

        buffer.push_back(c);

        // Process segment once full (or final partial segment)
        if (buffer.size() == L || (i == dna.size() - 1 && buffer.size() >= L / 2))
        {
            Segment seg;
            seg.sequence = buffer;

            // Select K by minimizing the combined P-value
            double bestP = 1.0;
            int bestK = -1;
            string bestWord;
            vector<array<int, 4>> bestCounts;

            // Evaluate candidate periods
            for (int k = KMIN; k <= KMAX; ++k)
            {
                int numKmers = buffer.size() / k;
                if (numKmers == 0)
                    continue;

                // Compute mismatch-aware representative word
                string w = ComputeRepresentativeWord_Mismatch(buffer, k, MISTAKE);

                // Enforce rotation invariance
                int shift = CanonicalRotationOffset(w);
                string canonW = CanonicalRotation(w);

                // Compute position-wise nucleotide counts
                auto counts = ComputePositionCounts(buffer, k);
                counts = RotateCountsLeft(counts, shift);

                // Compute position-wise binomial p-values
                auto pvals = ComputePValuesFromCounts(counts, numKmers, k);

                // Combine positional p-values using Fisher's method
                double combinedP = CombinePValuesFisher(pvals);

                // Retain the K and word with the strongest statistical signal
                if (combinedP < bestP)
                {
                    bestP = combinedP;
                    bestK = k;
                    bestWord = canonW;
                    seg.representativeWord = w;
                    bestCounts = counts;
                }
            }

            // Compute most frequent exact k-mers for comparison
            ComputeTwoMostFrequentWords(
                buffer,
                bestK,
                seg.frequentWord1,
                seg.frequentCount1,
                seg.frequentWord2,
                seg.frequentCount2);

            // Store canonical representative and statistics
            seg.canonicalRepresentativeWord = bestWord;
            seg.positioncounts = bestCounts;
            seg.combinedPValue = bestP;

            // AVG is reported as a descriptive metric (not selection criterion)
            seg.avg = ComputeAvgMatch1(buffer, seg.representativeWord);

            // Count exact representative occurrences
            seg.representativeCount =
                CountRepresentativeWord(buffer, seg.representativeWord);

            // Record genomic coordinates
            seg.length = buffer.size();
            seg.startIndex = i - buffer.size() + 1;

            segments.push_back(seg);
            buffer.clear();
        }
    }

    return segments;
}

// ============================================================================
// NOTE:
// This pipeline selects the period length K by maximizing the absolute
// occurrence count of the most frequent exact k-mer in the segment.
//
// This criterion is intentionally simple but biased toward smaller K values
// (e.g., K=3), since shorter k-mers naturally occur more frequently.
// A normalized alternative would be c1 / numKmers, but here we keep
// the absolute count to serve as a baseline reference pipeline.
// ============================================================================
vector<Segment> SegmentSequence3(const string &dna)
{
    vector<Segment> segments;
    string buffer;
    buffer.reserve(L);

    // Sequential scan over the DNA sequence
    for (size_t i = 0; i < dna.size(); ++i)
    {
        // Convert to uppercase and ignore non-ACGT characters
        char c = toupper(dna[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            continue;

        buffer.push_back(c);

        // Process a segment once it reaches length L (or final partial segment)
        if (buffer.size() == L || (i == dna.size() - 1 && buffer.size() >= L / 2))
        {
            Segment seg;
            seg.sequence = buffer;

            int bestK = -1;
            int bestCount1 = -1;
            string bestWordRaw; // Most frequent exact k-mer (non-canonical)

            // Test candidate period lengths
            for (int k = KMIN; k < KMAX; ++k)
            {
                int numKmers = (int)buffer.size() / k;
                if (numKmers <= 0)
                    continue;

                string w1, w2;
                int c1 = 0, c2 = 0;

                // Identify the two most frequent exact k-mers for this K
                ComputeTwoMostFrequentWords(buffer, k, w1, c1, w2, c2);

                // Select K by maximizing the absolute frequency of the top k-mer
                if (c1 > bestCount1)
                {
                    bestCount1 = c1;
                    bestK = k;
                    bestWordRaw = w1;
                }
            }

            // Skip segment if no valid K was found
            if (bestK == -1)
            {
                buffer.clear();
                continue;
            }

            // Recompute frequent words using the selected K
            // (ensures w2 and c2 are consistent with bestK)
            string w1, w2;
            int c1 = 0, c2 = 0;
            ComputeTwoMostFrequentWords(buffer, bestK, w1, c1, w2, c2);

            seg.frequentWord1 = w1;
            seg.frequentCount1 = c1;
            seg.frequentWord2 = w2;
            seg.frequentCount2 = c2;

            // Representative word is computed position-wise
            // using majority voting across aligned k-mers
            seg.representativeWord =
                ComputeRepresentativeWord1(seg.sequence, bestK);

            // Canonical form ensures rotation-invariant representation
            seg.canonicalRepresentativeWord =
                CanonicalRotation(bestWordRaw);

            // Exact occurrence count of the representative word
            seg.representativeCount =
                CountRepresentativeWord(buffer, seg.representativeWord);

            // Average positional match percentage (diagnostic metric)
            seg.avg =
                ComputeAvgMatch1(buffer, seg.representativeWord);

            // Align position-specific counts to canonical rotation
            int shift = CanonicalRotationOffset(seg.representativeWord);
            int numKmers = (int)buffer.size() / bestK;

            seg.positioncounts =
                ComputePositionCounts(buffer, bestK);
            seg.positioncounts =
                RotateCountsLeft(seg.positioncounts, shift);

            // Statistical significance (not used for selection in this pipeline)
            auto pvals =
                ComputePValuesFromCounts(seg.positioncounts, numKmers, bestK);
            seg.combinedPValue =
                CombinePValuesFisher(pvals);

            // Genomic coordinates
            seg.length = (int)buffer.size();
            seg.startIndex = (int)(i - buffer.size() + 1);

            segments.push_back(seg);
            buffer.clear();
        }
    }

    return segments;
}

// ============================================================================
// ComputeRepresentativeWord1
// ----------------------------------------------------------------------------
// Constructs a representative k-mer by majority voting at each position.
// For each position j in [0, k-1], the nucleotide with the highest frequency
// among aligned k-mers is selected.
//
// This method is fast (O(k * numKmers)) and robust to random noise,
// but does not explicitly model mismatches between k-mers.
// ============================================================================
string ComputeRepresentativeWord1(const string &segment, int k)
{
    string word;
    int numKmers = segment.size() / k;

    // Iterate over positions within the k-mer
    for (int pos = 0; pos < k; ++pos)
    {
        int countA = 0, countC = 0, countG = 0, countT = 0;

        // Count nucleotides at the given position across all k-mers
        for (int i = 0; i < numKmers; ++i)
        {
            char c = segment[i * k + pos];
            switch (c)
            {
            case 'A':
                countA++;
                break;
            case 'C':
                countC++;
                break;
            case 'G':
                countG++;
                break;
            case 'T':
                countT++;
                break;
            }
        }

        // Select the nucleotide with the maximum count
        char maxNuc = 'A';
        int maxCount = countA;

        if (countC > maxCount)
        {
            maxCount = countC;
            maxNuc = 'C';
        }
        if (countG > maxCount)
        {
            maxCount = countG;
            maxNuc = 'G';
        }
        if (countT > maxCount)
        {
            maxCount = countT;
            maxNuc = 'T';
        }

        word += maxNuc;
    }

    return word;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputePositionCounts
// ------------------------------------------------------------------------------------------
// Builds a position-specific nucleotide count table for a given segment and period length k.
//
// The segment is divided into aligned k-mers. For each position j ∈ [0, k−1],
// we count the number of occurrences of A, C, G, and T across all aligned k-mers.
//
// The output is a k × 4 table:
//   table[j][0] = count of 'A' at position j
//   table[j][1] = count of 'C' at position j
//   table[j][2] = count of 'G' at position j
//   table[j][3] = count of 'T' at position j
//
// This table is the core statistical representation used for significance testing.
// It captures positional bias indicative of hidden periodic structure.
// ------------------------------------------------------------------------------------------
vector<array<int, 4>> ComputePositionCounts(const string &segment, int k)
{
    vector<array<int, 4>> table(k, {0, 0, 0, 0});
    int numKmers = segment.size() / k;

    // Iterate over all aligned k-mers
    for (int i = 0; i < numKmers; ++i)
    {
        int base = i * k;

        // Count nucleotides position-wise
        for (int pos = 0; pos < k; ++pos)
        {
            char c = segment[base + pos];
            if (c == 'A')
                table[pos][0]++;
            else if (c == 'C')
                table[pos][1]++;
            else if (c == 'G')
                table[pos][2]++;
            else if (c == 'T')
                table[pos][3]++;
        }
    }
    return table;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputePValuesFromCounts
// ------------------------------------------------------------------------------------------
// Computes a positional P-value for each position in the k-mer using a binomial model.
//
// For each position j:
//   - We take the maximum nucleotide count among {A,C,G,T}
//   - Under the null hypothesis of random DNA (p = 0.25),
//     we compute the probability of observing at least this many matches
//     using a binomial tail probability.
//
// The result is a vector of k P-values, one per position.
// These P-values are later combined using Fisher's method.
// ------------------------------------------------------------------------------------------
vector<double> ComputePValuesFromCounts(
    const vector<array<int, 4>> &table,
    int numKmers, int k)
{
    vector<double> pValues(k);

    for (int pos = 0; pos < k; ++pos)
    {
        // Maximum nucleotide count at this position
        int maxCount = max({table[pos][0],
                            table[pos][1],
                            table[pos][2],
                            table[pos][3]});

        // Binomial tail probability under uniform background
        pValues[pos] = BinomialPValue(numKmers, maxCount, 0.25);
    }
    return pValues;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BinomialPValue
// ------------------------------------------------------------------------------------------
// Computes the probability of observing k or more successes out of n trials,
// assuming a binomial distribution with success probability p.
//
// Here, a "success" corresponds to matching a specific nucleotide at a given position.
// After computing the tail probability B for a single nucleotide,
// we apply a multiple-hypothesis correction over the four nucleotides
// by computing:
//
//     P = 1 − (1 − B)^4
//
// This yields a conservative position-wise P-value.
// ------------------------------------------------------------------------------------------
double BinomialPValue(int n, int k, double p)
{
    double tail = 0.0;

    // Compute binomial tail: P(X ≥ k)
    for (int i = k; i <= n; ++i)
        tail += BinomialCoefficient(n, i) *
                pow(p, i) *
                pow(1.0 - p, n - i);

    // Tail probability for a single nucleotide
    double B = tail;

    // Bonferroni-style correction over four nucleotides
    double P = 1.0 - pow(1.0 - B, 4.0);

    return P;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BinomialCoefficient
// ------------------------------------------------------------------------------------------
// Computes the binomial coefficient "n choose k" using a numerically stable
// multiplicative formulation.
//
// This implementation avoids factorials and reduces overflow risk.
// ------------------------------------------------------------------------------------------
double BinomialCoefficient(int n, int k)
{
    if (k < 0 || k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;

    // Use symmetry to reduce computation
    if (k > n / 2)
        k = n - k;

    double res = 1.0;
    for (int i = 1; i <= k; ++i)
        res = res * (n - i + 1) / i;

    return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CombinePValuesFisher
// ------------------------------------------------------------------------------------------
// Combines multiple independent P-values using Fisher's method.
//
// The test statistic is:
//     X = −2 Σ log(p_i)
// which follows a Chi-square distribution with 2k degrees of freedom,
// where k is the number of combined P-values.
//
// We compute the survival function P(ChiSq ≥ X), which serves as a single
// segment-level significance score.
//
// Note: This implementation returns the combined P-value directly,
// not the raw Chi-square statistic.
// ------------------------------------------------------------------------------------------
double CombinePValuesFisher(const vector<double> &pValues)
{
    int k = pValues.size();
    double X = 0.0;

    for (double p : pValues)
    {
        // Numerical safeguard against log(0)
        if (p > 1e-300)
            X += -2.0 * log(p);
        else
            X += -2.0 * log(1e-300);
    }

    // Compute survival probability for Chi-square with 2k degrees of freedom
    double p_combined = ChiSquareSurvivalEvenDF(X, 2 * k);
    return p_combined;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ChiSquareSurvivalEvenDF
// ------------------------------------------------------------------------------------------
// Computes the survival function P(ChiSquare(df) ≥ X) for an even number of degrees of freedom.
//
// For df = 2m (m ≥ 1), the Chi-square survival function admits a closed-form expression:
//
//   P(ChiSq ≥ X) = exp(−X/2) * Σ_{j=0}^{m−1} ( (X/2)^j / j! )
//
// This implementation avoids numerical integration and special functions,
// making it fast and stable for repeated use inside the pipeline.
//
// The function is used as the final step of Fisher's method to convert
// the aggregated test statistic into a combined P-value.
// ------------------------------------------------------------------------------------------
double ChiSquareSurvivalEvenDF(double X, int df)
{
    // Degrees of freedom must be positive and even
    if (df <= 0 || (df % 2) != 0)
        return 0.0; // invalid input

    int m = df / 2;
    double t = 0.5 * X;

    double term = 1.0; // t^0 / 0!
    double sum = 1.0;

    // Compute Σ_{j=0}^{m−1} (t^j / j!)
    for (int j = 1; j < m; ++j)
    {
        term *= t / j; // iterative update of t^j / j!
        sum += term;
    }

    // Survival probability
    double p = std::exp(-t) * sum;
    return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MergeSameWordSegments
// ------------------------------------------------------------------------------------------
// Merges consecutive segments that share the same canonical representative word
// and the same period length.
//
// The merge operation:
//   - Concatenates sequences
//   - Aggregates position-wise nucleotide counts
//   - Recomputes the combined P-value for the merged segment
//   - Updates the average match score as a length-weighted mean
//
// This step reduces artificial fragmentation caused by fixed window segmentation
// and produces longer, biologically coherent repeat regions.
// ------------------------------------------------------------------------------------------
vector<Segment> MergeSameWordSegments(const vector<Segment> &segments)
{
    vector<Segment> merged;
    if (segments.empty())
        return merged;

    Segment current = segments[0];

    for (size_t i = 1; i < segments.size(); ++i)
    {
        int K1 = segments[i].canonicalRepresentativeWord.size();
        int K2 = current.canonicalRepresentativeWord.size();

        // Merge only if canonical word and period length are identical
        if (segments[i].canonicalRepresentativeWord == current.canonicalRepresentativeWord && K1 == K2)
        {
            int oldLen = current.length;
            current.length += segments[i].length;

            // Concatenate sequences
            current.sequence += segments[i].sequence;

            // Aggregate positional counts
            AddCounts(current.positioncounts, segments[i].positioncounts);

            // Recompute statistical significance for the merged segment
            int numKmers = current.length / K1;
            vector<double> positionPValues =
                ComputePValuesFromCounts(current.positioncounts, numKmers, K1);
            current.combinedPValue = CombinePValuesFisher(positionPValues);

            // Update average match score (length-weighted)
            current.avg =
                (current.avg * oldLen + segments[i].avg * segments[i].length) /
                (double)current.length;
        }
        else
        {
            merged.push_back(current);
            current = segments[i];
        }
    }

    merged.push_back(current);
    return merged;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HammingDistance
// ------------------------------------------------------------------------------------------
// Computes the Hamming distance between two strings of equal length.
//
// If the lengths differ, the distance is defined as the maximum length,
// ensuring a strict mismatch penalty.
//
// This function is used in mismatch-aware repeat detection to quantify
// similarity between k-mers.
// ------------------------------------------------------------------------------------------
int HammingDistance(const string &a, const string &b)
{
    if (a.size() != b.size())
        return max(a.size(), b.size());

    int d = 0;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
            d++;

    return d;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CanonicalRotation
// ------------------------------------------------------------------------------------------
// Computes the canonical (lexicographically minimal) cyclic rotation of a word.
//
// Periodic repeats are invariant under cyclic shifts.
// By mapping each representative word to its canonical rotation,
// the algorithm enforces rotation-invariant comparison and merging.
//
// Example:
//   "GAT" → {"GAT","ATG","TGA"} → canonical = "ATG"
// ------------------------------------------------------------------------------------------
string CanonicalRotation(const string &w)
{
    string best = w;
    int k = w.size();

    for (int i = 1; i < k; ++i)
    {
        string rotated = w.substr(i) + w.substr(0, i);
        if (rotated < best)
            best = rotated;
    }

    return best;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ReadFNA
// ------------------------------------------------------------------------------------------
// Reads a DNA sequence from a FASTA (.fna) file.
//
// Lines starting with '>' are treated as headers and ignored.
// All remaining lines are concatenated into a single DNA string.
//
// The function assumes that the input file is valid FASTA format
// and terminates the program if the file cannot be opened.
// ------------------------------------------------------------------------------------------
string ReadFNA(const string &filename)
{
    ifstream file(filename);
    if (!file)
    {
        cerr << "Error: cannot open file " << filename << endl;
        exit(1);
    }

    string line;
    string dna;

    while (getline(file, line))
    {
        if (line.empty())
            continue;

        if (line[0] == '>')
            continue; // FASTA header

        dna += line;
    }

    file.close();
    return dna;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputeAvgMatch1
// ------------------------------------------------------------------------------------------
// Computes the average positional match percentage between the segment and
// the representative word, assuming a fixed alignment (no rotation).
//
// The segment is partitioned into non-overlapping k-mers of length K.
// For each k-mer, positions are compared directly against the representative word.
//
// Output:
//   Percentage of matching positions over all k-mers and all positions.
// ------------------------------------------------------------------------------------------
double ComputeAvgMatch1(const string &segment, const string &repWord)
{
    int K = repWord.size();
    int numKmers = segment.size() / K;

    // No valid k-mers → no meaningful match
    if (numKmers == 0)
        return 0.0;

    // Total number of aligned positions
    double totalPositions = numKmers * K;
    double matchCount = 0;

    // Compare each k-mer to the representative word position-wise
    for (int i = 0; i < numKmers; ++i)
    {
        for (int pos = 0; pos < K; ++pos)
        {
            if (segment[i * K + pos] == repWord[pos])
                matchCount++;
        }
    }

    // Return percentage of exact positional matches
    return 100.0 * matchCount / totalPositions;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputeAvgMatch2
// ------------------------------------------------------------------------------------------
// Computes the average positional match percentage allowing cyclic shifts
// of the representative word.
//
// For each k-mer in the segment:
//   - All K cyclic rotations of the representative word are tested
//   - The best local alignment (maximum number of matches) is selected
//
// This measure is rotation-invariant and is suitable for periodic repeats
// whose phase may vary across the segment.
// ------------------------------------------------------------------------------------------
double ComputeAvgMatch2(const string &segment, const string &repWord)
{
    int K = repWord.size();
    int numKmers = segment.size() / K;

    // Invalid input: no k-mers or empty word
    if (numKmers == 0 || K == 0)
        return 0.0;

    double totalPositions = numKmers * K;
    double totalMatch = 0.0;

    // Process each k-mer independently
    for (int i = 0; i < numKmers; ++i)
    {
        string km = segment.substr(i * K, K);
        int bestLocalMatch = 0;

        // Try all cyclic rotations of the representative word
        for (int shift = 0; shift < K; ++shift)
        {
            int localMatch = 0;

            for (int pos = 0; pos < K; ++pos)
            {
                char w = repWord[(pos + shift) % K];
                if (km[pos] == w)
                    localMatch++;
            }

            // Keep the best alignment over all rotations
            bestLocalMatch = max(bestLocalMatch, localMatch);
        }

        totalMatch += bestLocalMatch;
    }

    // Return rotation-invariant average match percentage
    return 100.0 * totalMatch / totalPositions;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RotateCountsLeft
// ------------------------------------------------------------------------------------------
// Performs a cyclic left rotation of position-wise nucleotide count tables.
//
// This operation is required after canonical rotation of the representative word
// in order to maintain alignment between the word and its positional statistics.
// ------------------------------------------------------------------------------------------
vector<array<int, 4>> RotateCountsLeft(const vector<array<int, 4>> &table, int shift)
{
    int k = (int)table.size();
    if (k == 0)
        return table;

    shift %= k;
    if (shift < 0)
        shift += k;

    vector<array<int, 4>> out(k);
    for (int pos = 0; pos < k; ++pos)
        out[pos] = table[(pos + shift) % k];

    return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CanonicalRotationOffset
// ------------------------------------------------------------------------------------------
// Returns the left rotation offset that transforms the word into its
// canonical (lexicographically minimal) rotation.
//
// The offset is later applied to position-wise statistics to ensure
// rotation-invariant analysis.
// ------------------------------------------------------------------------------------------
int CanonicalRotationOffset(const string &w)
{
    string best = w;
    int bestShift = 0;
    int k = (int)w.size();

    for (int s = 1; s < k; ++s)
    {
        string rotated = w.substr(s) + w.substr(0, s);
        if (rotated < best)
        {
            best = rotated;
            bestShift = s;
        }
    }

    return bestShift;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AddCounts
// ------------------------------------------------------------------------------------------
// Element-wise addition of two position-wise nucleotide count tables.
//
// This function is used when merging adjacent segments in order to
// accumulate positional statistics before recomputing significance.
// ------------------------------------------------------------------------------------------
void AddCounts(vector<array<int, 4>> &a, const vector<array<int, 4>> &b)
{
    for (int pos = 0; pos < (int)a.size(); ++pos)
        for (int j = 0; j < 4; ++j)
            a[pos][j] += b[pos][j];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputeLocalDistributions
// ------------------------------------------------------------------------------------------
// Computes local nucleotide probability distributions using a sliding window.
//
// For each position i:
//   - A window of size ±W is considered
//   - Nucleotide frequencies are estimated with optional pseudo-count smoothing
//   - The resulting probabilities are used for composition-preserving null models
//
// Output:
//   A per-position probability vector [pA, pC, pG, pT].
// ------------------------------------------------------------------------------------------
vector<array<double, 4>> ComputeLocalDistributions(
    const string &dna, int W, double pseudo)
{
    int n = dna.size();

    // probs[i] = local nucleotide probability distribution at position i
    // Index order: [A, C, G, T]
    vector<array<double, 4>> probs(n, {0, 0, 0, 0});

    // Sliding window nucleotide counts
    array<int, 4> cnt{0, 0, 0, 0};

    // Map nucleotide to index
    auto idx = [](char c)
    {
        c = toupper(c);
        if (c == 'A')
            return 0;
        if (c == 'C')
            return 1;
        if (c == 'G')
            return 2;
        if (c == 'T')
            return 3;
        return -1;
    };

    // Initialize sliding window boundaries
    int left = 0;
    int right = min(n - 1, W);

    // Initialize counts for the first window [0, W]
    for (int i = left; i <= right; ++i)
        if (idx(dna[i]) >= 0)
            cnt[idx(dna[i])]++;

    // Iterate over all positions in the sequence
    for (int i = 0; i < n; ++i)
    {
        // Compute smoothed local probabilities
        // P(nucleotide | local window centered at i)
        //
        // Pseudocount prevents zero probabilities and stabilizes estimation
        double sum = cnt[0] + cnt[1] + cnt[2] + cnt[3] + 4 * pseudo;
        for (int j = 0; j < 4; j++)
            probs[i][j] = (cnt[j] + pseudo) / sum;

        // Update window boundaries for the next position
        int newL = max(0, i + 1 - W);
        int newR = min(n - 1, i + 1 + W);

        // Remove nucleotides exiting the window on the left
        while (left < newL)
            if (idx(dna[left]) >= 0)
                cnt[idx(dna[left++])]--;

        // Add nucleotides entering the window on the right
        while (right < newR)
            if (idx(dna[++right]) >= 0)
                cnt[idx(dna[right])]++;
    }

    return probs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenerateRandomDNA
// ------------------------------------------------------------------------------------------
// Generates a random DNA sequence given per-position nucleotide probabilities.
//
// Each position i is sampled independently according to probs[i] = {pA, pC, pG, pT}.
// This function is typically used to generate null-model sequences that preserve
// local nucleotide composition estimated from real DNA.
//
// Input:
//   probs[i][j] = probability of nucleotide j at position i
//
// Output:
//   A random DNA string sampled from the given distributions
// ------------------------------------------------------------------------------------------
string GenerateRandomDNA(
    const vector<array<double, 4>> &probs)
{
    random_device rd;
    mt19937 rng(rd());
    uniform_real_distribution<double> U(0, 1);

    string s(probs.size(), 'A');
    for (size_t i = 0; i < probs.size(); i++)
    {
        double r = U(rng);
        double c = probs[i][0];

        if (r < c)
            s[i] = 'A';
        else if (r < (c += probs[i][1]))
            s[i] = 'C';
        else if (r < (c += probs[i][2]))
            s[i] = 'G';
        else
            s[i] = 'T';
    }

    return s;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RunPipeline
// ------------------------------------------------------------------------------------------
// Dispatch function that runs one of the segmentation pipelines based on mode.
//
// mode = 1 : AVG-based pipeline
// mode = 2 : P-value–based pipeline
// mode = 3 : Frequency-based pipeline (baseline)
//
// This function guarantees that the same interface is used for both
// real DNA and null-model DNA, ensuring fair comparison.
// ------------------------------------------------------------------------------------------
vector<Segment> RunPipeline(const string &dna, int mode)
{
    vector<Segment> segs;

    if (mode == 1)
    {
        segs = SegmentSequence1(dna);
    }
    else if (mode == 2)
    {
        segs = SegmentSequence2(dna);
    }
    else if (mode == 3)
    {
        segs = SegmentSequence3(dna);
    }
    else
    {
        cerr << "RunPipeline: invalid mode = " << mode << endl;
        return {};
    }

    return segs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenerateRandomDNA_Global
// ------------------------------------------------------------------------------------------
// Generates a composition-preserving random DNA sequence using a sliding window model.
//
// For each position i:
//   - A window of size ±W is maintained over the original DNA
//   - Local nucleotide frequencies are estimated with pseudo-count smoothing
//   - A nucleotide is sampled according to these local probabilities
//
// This produces a null model that preserves local compositional structure
// (e.g., GC-content and isochores) while destroying higher-order repeat patterns.
//
// Input:
//   dna    : original DNA sequence
//   W      : half-window size
//   pseudo : smoothing constant to avoid zero probabilities
//
// Output:
//   A random DNA sequence of the same length as the input
// ------------------------------------------------------------------------------------------
string GenerateRandomDNA_Global(
    const string &dna,
    int W,
    double pseudo = 1)
{
    int n = dna.size();

    // Output random DNA sequence (same length as input)
    string rnd;
    rnd.reserve(n);

    // Sliding window nucleotide counts: [A, C, G, T]
    array<int, 4> cnt{0, 0, 0, 0};

    // Map nucleotide character to index
    auto idx = [](char c)
    {
        if (c == 'A')
            return 0;
        if (c == 'C')
            return 1;
        if (c == 'G')
            return 2;
        if (c == 'T')
            return 3;
        return -1;
    };

    // Initialize counts for the first window [0, W]
    // This seeds the local composition for the beginning of the sequence
    for (int i = 0; i <= W && i < n; i++)
        if (idx(dna[i]) >= 0)
            cnt[idx(dna[i])]++;

    // Random number generator for multinomial sampling
    random_device rd;
    mt19937 rng(rd());
    uniform_real_distribution<double> U(0, 1);

    // Iterate over genome positions
    for (int i = 0; i < n; i++)
    {
        // Compute smoothed local nucleotide probabilities
        // P(X = nucleotide | local window)
        // Pseudocount prevents zero-probability artifacts
        double sum = cnt[0] + cnt[1] + cnt[2] + cnt[3] + 4 * pseudo;

        array<double, 4> p{
            (cnt[0] + pseudo) / sum,
            (cnt[1] + pseudo) / sum,
            (cnt[2] + pseudo) / sum,
            (cnt[3] + pseudo) / sum};

        // Sample a nucleotide according to the local distribution
        double r = U(rng);
        char c;
        if (r < p[0])
            c = 'A';
        else if (r < p[0] + p[1])
            c = 'C';
        else if (r < p[0] + p[1] + p[2])
            c = 'G';
        else
            c = 'T';

        rnd.push_back(c);

        // Slide the window forward:
        // remove nucleotide exiting on the left
        // add nucleotide entering on the right
        int L = i - W;
        int R = i + W + 1;

        if (L >= 0 && idx(dna[L]) >= 0)
            cnt[idx(dna[L])]--;

        if (R < n && idx(dna[R]) >= 0)
            cnt[idx(dna[R])]++;
    }

    return rnd;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AnnotateSegmentsWithNull1
// ------------------------------------------------------------------------------------------
// Annotates real segments with an empirical significance score using a null model.
//
// The null model preserves local nucleotide composition using a sliding window,
// ensuring that detected repeat signals are not explained by GC-content or
// compositional bias alone.
//
// The strengthScore is defined as:
//     strengthScore = p_null / p_real
//
// Interpretation:
//   - strong  : real segment is much more significant than null (ratio ≥ 10)
//   - weak    : real segment slightly more significant than null (1 < ratio < 10)
//   - noise   : no evidence beyond background (ratio ≤ 1)
// ------------------------------------------------------------------------------------------
vector<Segment> AnnotateSegmentsWithNull1(
    vector<Segment> &realSegs,
    const string &dna,
    int W, int mode)
{
    // Step 1: Generate composition-preserving random DNA (local null model)
    string rndDNA = GenerateRandomDNA_Global(dna, W);
    cout << "GenerateRandomDNA done: " << endl;

    // Step 2: Run the exact same detection pipeline on the null sequence
    vector<Segment> nullSegs = RunPipeline(rndDNA, mode);
    cout << "random segment done: " << endl;

    // Step 3: One-to-one comparison between real and null segments
    int N = min(realSegs.size(), nullSegs.size());

    for (int i = 0; i < N; ++i)
    {
        // Combined statistical significance for real and null segments
        double pReal = realSegs[i].combinedPValue;
        double pNull = nullSegs[i].combinedPValue;

        // Empirical significance ratio (stabilized to avoid division by zero)
        double ratio = max(pNull, 1e-300) / max(pReal, 1e-300);

        realSegs[i].strengthScore = ratio;

        // Discrete classification based on empirical evidence
        if (ratio >= 10)
            realSegs[i].strength = "strong";
        else if (ratio > 1)
            realSegs[i].strength = "weak";
        else
            realSegs[i].strength = "noise";
    }

    return nullSegs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SplitByLength
// ------------------------------------------------------------------------------------------
// Splits the DNA sequence into non-overlapping chunks of fixed length L.
// Mainly used for baseline analysis or alternative segmentation strategies.
// ------------------------------------------------------------------------------------------
vector<string> SplitByLength(const string &dna)
{
    vector<string> chunks;
    for (int i = 0; i + L <= (int)dna.size(); i += L)
        chunks.push_back(dna.substr(i, L));
    return chunks;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FormatSequenceByK
// ------------------------------------------------------------------------------------------
// Formats a DNA sequence by inserting a space every K characters.
// This function is used strictly for visualization and debugging.
// ------------------------------------------------------------------------------------------
string FormatSequenceByK(const string &seq, int K)
{
    string out;
    for (size_t i = 0; i < seq.size(); ++i)
    {
        if (i > 0 && i % K == 0)
            out += ' ';
        out += seq[i];
    }
    return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ExportSegmentsCSV_Simple
// ------------------------------------------------------------------------------------------
// Exports segment-level statistics into a CSV file for downstream analysis.
//
// The output is designed to support:
//   - Statistical comparison (real vs null)
//   - Visualization (histograms, scatter plots, CDFs)
//   - Agreement analysis between representative and frequent words
// ------------------------------------------------------------------------------------------
void ExportSegmentsCSV_Simple(
    const vector<Segment> &segs,
    const string &filename,
    const string &pipeline,
    const string &type) // "real" or "null"
{
    // Open output file stream
    ofstream out(filename);
    if (!out)
    {
        cerr << "Cannot write " << filename << endl;
        return;
    }

    // Write CSV header:
    // - pipeline: algorithmic pipeline used (avg / pvalue)
    // - type: real genomic data or composition-preserving null model
    // - segment_id: index of the segment in the sequence
    // - start, length: genomic coordinates
    // - AVG: average positional match percentage
    // - pvalue: combined Fisher p-value
    // - strengthScore: empirical null/real significance ratio
    // - class: qualitative classification (noise / weak / strong)
    // - representativeWord + count: consensus repeat word and its exact support
    // - frequentWord{1,2} + count: most frequent observed k-mers
    out << "pipeline,type,segment_id,start,length,"
           "AVG,pvalue,strengthScore,class,"
           "representativeWord,representativeCount,"
           "frequentWord1,frequentCount1,"
           "frequentWord2,frequentCount2\n";

    // Write one row per segment
    for (int i = 0; i < (int)segs.size(); i++)
    {
        const Segment &s = segs[i];

        out << pipeline << ","              // Algorithmic pipeline identifier
            << type << ","                  // Data type: real or null
            << i << ","                     // Segment index
            << s.startIndex << ","          // Start position in genome
            << s.length << ","              // Segment length
            << s.avg << ","                 // Average positional match (%)
            << s.combinedPValue << ","      // Combined statistical significance
            << s.strengthScore << ","       // Empirical null/real score
            << s.strength << ","            // Discrete classification
            << s.representativeWord << ","  // Representative repeat word
            << s.representativeCount << "," // Exact occurrence count
            << s.frequentWord1 << ","       // Most frequent observed k-mer
            << s.frequentCount1 << ","      // Its count
            << s.frequentWord2 << ","       // Second most frequent k-mer
            << s.frequentCount2             // Its count
            << "\n";
    }

    // Close file
    out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputeTwoMostFrequentWords
// ------------------------------------------------------------------------------------------
// Identifies the two most frequent non-overlapping k-mers in a segment.
// This provides a purely frequency-based baseline for comparison with
// the representative word selected by the algorithm.
// ------------------------------------------------------------------------------------------
void ComputeTwoMostFrequentWords(
    const string &segment,
    int K,
    string &word1, int &count1,
    string &word2, int &count2)
{
    unordered_map<string, int> freq;
    int numKmers = segment.size() / K;

    // Count occurrences of all k-mers
    for (int i = 0; i < numKmers; ++i)
    {
        string km = segment.substr(i * K, K);
        freq[km]++;
    }

    // Initialize outputs
    word1 = "";
    word2 = "";
    count1 = 0;
    count2 = 0;

    // Track the two most frequent k-mers
    for (auto &p : freq)
    {
        if (p.second > count1)
        {
            word2 = word1;
            count2 = count1;

            word1 = p.first;
            count1 = p.second;
        }
        else if (p.second > count2 && p.first != word1)
        {
            word2 = p.first;
            count2 = p.second;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CountRepresentativeWord
// ------------------------------------------------------------------------------------------
// Counts exact occurrences of the representative word in the segment.
// Used for agreement analysis and validation against frequent k-mers.
// ------------------------------------------------------------------------------------------
int CountRepresentativeWord(
    const string &segment,
    const string &repWord)
{
    // Length of the representative k-mer
    int K = repWord.size();

    // Guard against invalid representative words
    if (K == 0)
        return 0;

    // Number of non-overlapping k-mers in the segment
    int numKmers = segment.size() / K;

    // Counter for exact matches
    int count = 0;

    // Scan all k-mers in the segment
    for (int i = 0; i < numKmers; ++i)
    {
        // Compare the current k-mer to the representative word
        // using an exact string comparison (Hamming distance = 0)
        if (segment.compare(i * K, K, repWord) == 0)
            count++;
    }

    // Return the total number of exact matches
    return count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CountMatchesWithMismatches
// ------------------------------------------------------------------------------------------
// Counts how many k-mers in the segment match a candidate word
// with at most M mismatches (Hamming distance ≤ M).
//
// This enables robust detection of mutated or noisy repeats.
// ------------------------------------------------------------------------------------------
int CountMatchesWithMismatches(
    const string &segment,
    const string &word,
    int M)
{
    // Length of the k-mer (period length)
    int K = word.size();

    // Number of non-overlapping k-mers in the segment
    int numKmers = segment.size() / K;

    // Counter for k-mers that satisfy the mismatch constraint
    int count = 0;

    // Iterate over all k-mers in the segment
    for (int i = 0; i < numKmers; ++i)
    {
        int mism = 0; // Number of mismatches for the current k-mer

        // Compare the current k-mer to the candidate word position-wise
        for (int pos = 0; pos < K; ++pos)
        {
            if (segment[i * K + pos] != word[pos])
            {
                mism++;

                // Early termination: no need to continue if mismatches exceed M
                if (mism > M)
                    break;
            }
        }

        // Count this k-mer if it matches within the allowed mismatch budget
        if (mism <= M)
            count++;
    }

    // Return the total number of approximate matches
    return count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ComputeRepresentativeWord_Mismatch
// ------------------------------------------------------------------------------------------
// Selects a representative k-mer by maximizing the number of segment k-mers
// that match it with at most M mismatches.
//
// Candidate words are restricted to k-mers that actually appear in the segment,
// ensuring computational efficiency and biological plausibility.
// ------------------------------------------------------------------------------------------
string ComputeRepresentativeWord_Mismatch(
    const string &segment,
    int K,
    int M)
{
    // Number of non-overlapping k-mers in the segment
    int numKmers = segment.size() / K;

    // Edge case: segment too short to contain even one k-mer
    if (numKmers == 0)
        return "";

    // Collect all observed k-mers in the segment.
    // Candidate representative words are restricted to observed k-mers
    // to avoid introducing artificial or biologically implausible patterns.
    vector<string> kmers;
    for (int i = 0; i < numKmers; ++i)
        kmers.push_back(segment.substr(i * K, K));

    // Track the best representative word and its support count
    string bestWord;
    int bestCount = -1;

    // For each candidate k-mer, count how many segment k-mers
    // match it with at most M mismatches
    for (const string &cand : kmers)
    {
        int cnt = CountMatchesWithMismatches(segment, cand, M);

        // Select the k-mer that maximizes the number of approximate matches
        if (cnt > bestCount)
        {
            bestCount = cnt;
            bestWord = cand;
        }
    }

    // Return the representative word with maximal approximate support
    return bestWord;
}
