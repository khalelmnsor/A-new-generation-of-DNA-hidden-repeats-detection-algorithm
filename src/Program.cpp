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

using namespace std;
using CountTable = vector<map<char, int>>;

// Parameters
const int L = 200;        // Segment length
const double tau1 = 0.25; // High p-value threshold
const double tau2 = 0.05; // Low p-value threshold
// const double alpha = 0.05;

struct Segment
{
    string sequence;
    string canonicalRepresentativeWord;
    string representativeWord;
    vector<array<int, 4>> positioncounts;
    double combinedPValue;
    double pNull;
    double nullParallelP;
    int startIndex;
    int length;
    double avg = 0.0;
    string noise = "no";
    double strengthScore = 0.0;  // למשל -log10(empiricalP)
    string strength = "unknown"; // noise / weak / strong
};

// Function declarations
vector<Segment> SegmentSequence1(const string &dna);
vector<Segment> SegmentSequence2(const string &dna);

string ComputeRepresentativeWord1(const string &segment, int K);
string ComputeRepresentativeWord2(const string &segment, int k);
string ComputeRepresentativeWord3(const string &segment, int k);

double ComputeAvgMatch1(const string &segment, const string &repWord);
double ComputeAvgMatch2(const string &segment, const string &repWord);

vector<array<int, 4>> ComputePositionCounts(const string &segment, int k);
vector<double> ComputePValuesFromCounts(const vector<array<int, 4>> &table, int numKmers, int k);
vector<double> ComputePositionPValues(const string &segment, int k);

double BinomialPValue(int n, int k, double p);
double BinomialCoefficient(int n, int k);
double CombinePValuesFisher(const vector<double> &pValues);
vector<Segment> MergeSameWordSegments(const vector<Segment> &segments);
vector<Segment> MergeNoiseSegments(const vector<Segment> &segments);
double ChiSquareSurvivalEvenDF(double X, int df);                         //
string FilterACGT(const string &s);                                       //
vector<Segment> SplitIncoherentSegments(const vector<Segment> &segments); //
string CanonicalRotation(const string &w);
string ReadFNA(const string &filename);
string NormalizeAndCompressN(const string &s);

vector<array<int, 4>> RotateCountsLeft(const vector<array<int, 4>> &table, int shift);
int CanonicalRotationOffset(const string &w);
void AddCounts(vector<array<int, 4>> &a, const vector<array<int, 4>> &b);

vector<array<double, 4>> ComputeLocalDistributions(const string &dna, int W, double pseudo = 0);
string GenerateRandomDNA(const vector<array<double, 4>> &probs);
vector<Segment> RunPipeline(const string &dna, int mode);
vector<Segment> AnnotateSegmentsWithNull1(vector<Segment> &realSegs, const string &dna, int W, int mode);
void AnnotateSegmentsWithNull2(vector<Segment> &realSegs, const string &dna, int W);
string FormatSequenceByK(const string &seq, int K);
vector<string> SplitByLength(const string &dna);

int main()
{
    // FIX: Removed '+' signs to allow implicit string literal concatenation
    /*string dna = "GTGGCGGTGTTG"  // strong repeat GTGs
                 "GCGATTGGAATG"  // weak noise ATG
                 "GTGATGGTGGAG"; // strong repeat GTG again*/

    /* string dna = "GTGGGTGGGTGG"
                      "ACGTTCAGATCG"
                      "GTGGGTGGGTGG";*/

    // string dna = "CATACCCCACATCAACATGGAGATCACACACCTCAACACGGAGATCACACACATCAACAT"; // eror ttts
    // string dna = "ACGT GTCA GATCTGACGTTCGAGTCTAGCTAGTCGATCGATGCTAGTCGATCGTAGCTAG"; // eror4*
    string dna = ReadFNA("C:/Users/user/Desktop/project/Data/gene.fna");

    // string dna = "atgtgagat";
    //   dna = NormalizeAndCompressN(dna);

    cout
        << "Loaded sequence length: " << dna.length() << endl;
    // cout << "DNA Sequence Preview: " << dna << endl;s

    // 1. Segment sequence
    vector<Segment> segments = SegmentSequence1(dna);
    //////////////////
    vector<Segment> nullSegs = AnnotateSegmentsWithNull1(segments, dna, L / 10, 1);

    // 2. Merge same-word segments
    // segments = MergeSameWordSegments(segments); ////////////////////////////////

    // 3. Merge weak/noisy segments
    // segments = MergeNoiseSegments(segments);/////////////////////////////////////////
    // cout << "S: " << segments[1].positionPValues[0];
    // segments = SplitIncoherentSegments(segments);
    // Output results
    cout << "\nDetected Hidden Repeat Segments:\n";
    for (size_t i = 0; i < segments.size() && i < nullSegs.size(); ++i)
    {
        Segment &s = segments[i];
        Segment &nulls = nullSegs[i];
        cout << "====================================================\n";
        cout << "[SEGMENT " << i << "]\n";

        cout << "Location   : start=" << s.startIndex
             << ", length=" << s.length << "\n";

        cout << "Word(real) : " << s.representativeWord
             << "  | canonical: " << s.canonicalRepresentativeWord << "\n";

        cout << "Word(null) : " << nulls.representativeWord
             << "  | canonical: " << nulls.canonicalRepresentativeWord << "\n";

        cout << "P-values   : real=" << s.combinedPValue
             << " | null(parallel)=" << nulls.combinedPValue
             << " | null(K-matched)=" << s.pNull << "\n";

        cout << "Strength   : score=" << s.strengthScore
             << " | class=" << s.strength
             << " | merge-noise=" << s.noise << "\n";

        cout << "AVG match  : " << s.avg << "%\n";

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

    return 0;
}
// לפי AVG
vector<Segment> SegmentSequence1(const string &dna)
{
    vector<Segment> segments;
    string buffer;
    buffer.reserve(L);
    int bufferStart = 0;
    for (size_t i = 0; i < dna.size(); ++i)
    {
        char c = toupper(dna[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            continue;

        if (buffer.empty())
            buffer.reserve(L);

        buffer.push_back(c);

        if (buffer.size() == L || (i == dna.size() - 1 && buffer.size() >= L / 2))
        {
            Segment seg;
            seg.sequence = buffer;
            string w;
            string bestW;
            double bestAvg = -1;
            int bestK = -1;
            for (int k = 3; k < 9; k++)
            {
                w = ComputeRepresentativeWord1(buffer, k);
                double avg = ComputeAvgMatch1(buffer, w);
                if (avg > bestAvg)
                {
                    bestAvg = avg;
                    bestW = w;
                    bestK = k;
                }
            }
            seg.representativeWord = bestW;
            seg.canonicalRepresentativeWord = CanonicalRotation(bestW);

            seg.avg = bestAvg;
            int shift = CanonicalRotationOffset(bestW);
            seg.canonicalRepresentativeWord = CanonicalRotation(bestW);

            int numKmers = seg.sequence.size() / bestK;
            seg.positioncounts = ComputePositionCounts(seg.sequence, bestK);
            seg.positioncounts = RotateCountsLeft(seg.positioncounts, shift);
            vector<double> positionPValues = ComputePValuesFromCounts(seg.positioncounts, numKmers, bestK);
            seg.combinedPValue = CombinePValuesFisher(positionPValues);

            seg.startIndex = i - buffer.size() + 1;
            seg.length = buffer.size();

            segments.push_back(seg);
            buffer.clear();
        }
    }

    return segments;
}
// לפי PVALUE
vector<Segment> SegmentSequence2(const string &dna)
{
    vector<Segment> segments;
    string buffer;
    buffer.reserve(L);

    for (size_t i = 0; i < dna.size(); ++i)
    {
        char c = toupper(dna[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            continue;

        buffer.push_back(c);

        if (buffer.size() == L || (i == dna.size() - 1 && buffer.size() >= L / 2))
        {
            Segment seg;
            seg.sequence = buffer;

            // ===== בחירת K לפי P-VALUE =====
            double bestP = 1.0;
            int bestK = -1;
            string bestWord;
            vector<array<int, 4>> bestCounts;

            for (int k = 3; k <= 6; ++k)
            {
                int numKmers = buffer.size() / k;
                if (numKmers == 0)
                    continue;

                // מילה מייצגת
                string w = ComputeRepresentativeWord1(buffer, k);

                // Canonical rotation
                int shift = CanonicalRotationOffset(w);
                string canonW = CanonicalRotation(w);

                // counts
                auto counts = ComputePositionCounts(buffer, k);
                counts = RotateCountsLeft(counts, shift);

                // P-values
                auto pvals = ComputePValuesFromCounts(counts, numKmers, k);
                double combinedP = CombinePValuesFisher(pvals);

                if (combinedP < bestP)
                {
                    bestP = combinedP;
                    bestK = k;
                    bestWord = canonW;
                    seg.representativeWord = w;
                    bestCounts = counts;
                }
            }

            // ===== מילוי הסגמנט =====
            seg.canonicalRepresentativeWord = bestWord;
            seg.positioncounts = bestCounts;
            seg.combinedPValue = bestP;
            seg.avg = ComputeAvgMatch1(buffer, seg.representativeWord);

            seg.length = buffer.size();
            seg.startIndex = i - buffer.size() + 1;

            segments.push_back(seg);
            buffer.clear();
        }
    }

    return segments;
}

string FilterACGT(const string &s)
{
    string clean;
    clean.reserve(s.size());

    for (char c : s)
    {
        c = toupper(c);

        if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
            clean.push_back(c);
    }
    return clean;
}

// Compute representative k-mer word (robust & faster)
string ComputeRepresentativeWord1(const string &segment, int k)
{
    string word;
    int numKmers = segment.size() / k;

    for (int pos = 0; pos < k; ++pos)
    {
        int countA = 0, countC = 0, countG = 0, countT = 0;

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

string ComputeRepresentativeWord3(const string &segment, int k)
{
    string word;
    int numKmers = segment.size() / k;

    for (int pos = 0; pos < k; ++pos)
    {
        int countA = 0, countC = 0, countG = 0, countT = 0;

        for (int i = 0; i < numKmers; ++i)
        {
            int base = i * k;

            // חילוץ ה-k-mer
            string km = segment.substr(base, k);

            // חישוב ההזזה הקנונית
            int shift = CanonicalRotationOffset(km);

            // האות בפוזיציה אחרי סיבוב
            char c = km[(pos + shift) % k];

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

string ComputeRepresentativeWord2(const string &segment, int K)
{
    unordered_map<string, int> freq;

    for (int i = 0; i + K <= (int)segment.size(); ++i)
    {
        string km = segment.substr(i, K);
        string canon = CanonicalRotation(km);
        freq[canon]++;
    }

    string best;
    int bestCnt = -1;

    for (auto &p : freq)
    {
        if (p.second > bestCnt)
        {
            bestCnt = p.second;
            best = p.first;
        }
    }

    return best;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<array<int, 4>> ComputePositionCounts(const string &segment, int k)
{
    vector<array<int, 4>> table(k, {0, 0, 0, 0});
    int numKmers = segment.size() / k;

    for (int i = 0; i < numKmers; ++i)
    {
        int base = i * k;
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

vector<double> ComputePValuesFromCounts(
    const vector<array<int, 4>> &table,
    int numKmers, int k)
{
    vector<double> pValues(k);

    for (int pos = 0; pos < k; ++pos)
    {
        int maxCount = max({table[pos][0],
                            table[pos][1],
                            table[pos][2],
                            table[pos][3]});

        pValues[pos] = BinomialPValue(numKmers, maxCount, 0.25);
    }
    return pValues;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute p-values for each position
vector<double> ComputePositionPValues(const string &segment, int k)
{
    vector<double> pValues(k);
    int numKmers = segment.size() / k;

    for (int pos = 0; pos < k; ++pos)
    {
        map<char, int> counts{{'A', 0}, {'C', 0}, {'G', 0}, {'T', 0}};

        for (int i = 0; i < numKmers; ++i)
        {
            int idx = i * k + pos;
            if (idx >= segment.size())
                continue;
            counts[segment[idx]]++;
        }
        int maxCount = 0;
        for (auto &kv : counts)
            maxCount = max(maxCount, kv.second);
        pValues[pos] = BinomialPValue(numKmers, maxCount, 0.25);
    }
    return pValues;
}

double BinomialPValue(int n, int k, double p)
{
    double tail = 0.0;
    for (int i = k; i <= n; ++i)
        tail += BinomialCoefficient(n, i) * pow(p, i) * pow(1.0 - p, n - i);

    // B = P(X >= k) for one nucleotide
    double B = tail;

    // Corrected combined P-value over 4 nucleotides:
    double P = 1.0 - pow(1.0 - B, 4.0);

    return P;
}

// Binomial coefficient nCk
double BinomialCoefficient(int n, int k)
{
    if (k < 0 || k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
    if (k > n / 2)
        k = n - k;

    double res = 1.0;
    for (int i = 1; i <= k; ++i)
        res = res * (n - i + 1) / i;
    return res;
}

// Combine p-values using Fisher's method
double CombinePValuesFisher(const vector<double> &pValues)
{

    int k = pValues.size();
    double X = 0.0;
    for (double p : pValues)
    {
        // Avoid log(0) if p is extremely small
        if (p > 1e-300)
            X += -2.0 * log(p);
        else
            X += -2.0 * log(1e-300);
    }
    // Note: This simplified return acts as a score, not a direct Chi-sq p-value

    double p_combined = ChiSquareSurvivalEvenDF(X, 2 * k);
    return p_combined;
}

// Computes P(ChiSq(df) >= X) for even df using the closed-form series.
double ChiSquareSurvivalEvenDF(double X, int df)
{
    // df must be positive and even
    if (df <= 0 || (df % 2) != 0)
        return 0.0; // or handle error

    int m = df / 2;
    double t = 0.5 * X;

    double term = 1.0; // (t^0 / 0!)
    double sum = 1.0;

    for (int j = 1; j < m; ++j)
    {
        term *= t / j; // (t^j / j!) = (t^(j-1)/(j-1)!) * t/j
        sum += term;
    }

    double p = std::exp(-t) * sum; // e^{-t} * Σ_{j=0}^{m-1} t^j/j!
    return p;
}

// Merge consecutive segments with same word
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

        if (segments[i].canonicalRepresentativeWord == current.canonicalRepresentativeWord && K1 == K2)
        {
            // קודם להגדיל אורך
            int oldLen = current.length;
            current.length += segments[i].length;

            // לחבר מחרוזת אם אתה רוצה (לא חובה לסטטיסטיקה, אבל בסדר)
            current.sequence += segments[i].sequence;

            // לחבר counts
            AddCounts(current.positioncounts, segments[i].positioncounts);

            // עכשיו numKmers הנכון לפי האורך החדש
            int numKmers = current.length / K1;

            vector<double> positionPValues =
                ComputePValuesFromCounts(current.positioncounts, numKmers, K1);
            current.combinedPValue = CombinePValuesFisher(positionPValues);

            // avg
            current.avg = (current.avg * oldLen + segments[i].avg * segments[i].length) / (double)current.length;
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

// Parameters
/*const int L = 12;         // Segment length
const int K = 3;          // k-mer size
const double tau1 = 1;    // High p-value threshold
const double tau2 = 0.01; // Low p-value threshold
const double alpha = 0.05;*/

// Merge noisy segments
/*vector<Segment> MergeNoiseSegments(const vector<Segment> &segments)
{
    vector<Segment> result;
    if (segments.empty())
        return result;

    // We need a copy of segments because we might modify flow,
    // but here we are just iterating.
    // Since we are "skipping" indices, a while loop is appropriate.

    result.push_back(segments[0]); // Start with the first

    size_t i = 1;
    while (i < segments.size())
    {
        // We look at the *last pushed* segment (left) and the current one (middle)
        // We need to check if there is a 'right' segment to bridge to.
        if (i < segments.size() - 1)
        {
            Segment &left = result.back();
            const Segment &middle = segments[i];
            const Segment &right = segments[i + 1];

            bool leftStrong = left.combinedPValue < tau2;
            bool rightStrong = right.combinedPValue < tau2;
            bool middleWeak = middle.combinedPValue > tau1; // High p-value = weak signal

            // Check if middle is noise between two strong signals
            // Note: Original logic checked words weren't equal to middle,
            // implying middle is an interruption.
            if (middle.representativeWord != left.representativeWord &&
                middle.representativeWord != right.representativeWord &&
                left.representativeWord == right.representativeWord && /************
                leftStrong && rightStrong && middleWeak)
            {
                // Merge Middle and Right into Left
                left.sequence += middle.sequence + right.sequence;
                left.length = left.sequence.size();

                int numKmers = left.sequence.size() / K;
                left.positioncounts = ComputePositionCounts(left.sequence);
                vector<double> positionPValues = ComputePValuesFromCounts(left.positioncounts, numKmers);
                left.combinedPValue = CombinePValuesFisher(positionPValues);
                left.noise = "yes";

                i += 2; // Skip middle and right
                continue;
            }
        }
        // Otherwise just add the current segment
        result.push_back(segments[i]);
        i++;
    }
    return result;
}*/

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

/*vector<Segment> SplitIncoherentSegments(const vector<Segment> &segments)
{
    vector<Segment> result;

    const int W = K * 5; // חלון קטן
    // const int WORD_DIFF = 0; // סף שינוי במילה

    for (const Segment &seg : segments)
    {
        if ((int)seg.sequence.size() < 2 * W)
        {
            result.push_back(seg);
            continue;
        }

        int lastCut = 0;

        string prevWindow = seg.sequence.substr(0, W);
        string RepresentativeWordprev = ComputeRepresentativeWord(prevWindow);

        string prevWord = CanonicalRotation(RepresentativeWordprev);
        double prevP = CombinePValuesFisher(
            ComputePositionPValues(prevWindow));

        for (int i = W; i + W <= (int)seg.sequence.size(); i += W)
        {
            string currWindow = seg.sequence.substr(i, W);
            string RepresentativeWordcurr = ComputeRepresentativeWord(currWindow);

            string currWord = CanonicalRotation(RepresentativeWordcurr);
            double currP = CombinePValuesFisher(
                ComputePositionPValues(currWindow));

            bool bothStrong = (prevP < tau2 && currP < tau2);

            if (bothStrong)
            {
                Segment left;
                left.sequence = seg.sequence.substr(lastCut, i - lastCut);
                left.startIndex = seg.startIndex + lastCut;
                left.length = left.sequence.size();
                string w = ComputeRepresentativeWord(left.sequence);
                left.representativeWord = CanonicalRotation(w);
                left.avg = ComputeAvgMatch(left.sequence, w);

                int numKmers = left.sequence.size() / K;
                left.positioncounts = ComputePositionCounts(left.sequence);
                vector<double> positionPValues = ComputePValuesFromCounts(left.positioncounts, numKmers);
                left.combinedPValue = CombinePValuesFisher(positionPValues);

                result.push_back(left);
                lastCut = i;
            }

            prevWord = currWord;
            prevP = currP;
        }

        // הזנב
        if (lastCut < (int)seg.sequence.size())
        {
            Segment tail;
            tail.sequence = seg.sequence.substr(lastCut);
            tail.startIndex = seg.startIndex + lastCut;
            tail.length = tail.sequence.size();
            string w = ComputeRepresentativeWord(tail.sequence);
            tail.representativeWord = CanonicalRotation(w);
            tail.avg = ComputeAvgMatch(tail.sequence, w);

            int numKmers = tail.sequence.size() / K;
            tail.positioncounts = ComputePositionCounts(tail.sequence);
            vector<double> positionPValues = ComputePValuesFromCounts(tail.positioncounts, numKmers);
            tail.combinedPValue = CombinePValuesFisher(positionPValues);

            result.push_back(tail);
        }
    }

    return result;
}*/

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

#include <fstream>

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

string NormalizeAndCompressN(const string &s)
{
    string out;
    out.reserve(s.size());

    bool lastWasN = false;

    for (char c : s)
    {
        c = toupper(c);

        if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
        {
            out.push_back(c);
            lastWasN = false;
        }
        else if (c == 'N')
        {
            if (!lastWasN)
            {
                out.push_back('N'); // גבול אחד
                lastWasN = true;
            }
        }
        // כל תו אחר – מדולג
    }

    return out;
}

double ComputeAvgMatch1(const string &segment, const string &repWord)
{
    int K = repWord.size();
    int numKmers = segment.size() / K;

    if (numKmers == 0)
        return 0.0;

    double totalPositions = numKmers * K;
    double matchCount = 0;

    for (int i = 0; i < numKmers; ++i)
    {
        for (int pos = 0; pos < K; ++pos)
        {
            if (segment[i * K + pos] == repWord[pos])
                matchCount++;
        }
    }

    return 100.0 * matchCount / totalPositions;
}

double ComputeAvgMatch2(const string &segment, const string &repWord)
{
    int K = repWord.size();
    int numKmers = segment.size() / K;

    if (numKmers == 0 || K == 0)
        return 0.0;

    double totalPositions = numKmers * K;
    double totalMatch = 0.0;

    // עבור כל k-mer בסיגמנט
    for (int i = 0; i < numKmers; ++i)
    {
        string km = segment.substr(i * K, K);

        int bestLocalMatch = 0;

        // עבור כל סיבוב של המילה
        for (int shift = 0; shift < K; ++shift)
        {
            int localMatch = 0;

            for (int pos = 0; pos < K; ++pos)
            {
                char w = repWord[(pos + shift) % K];
                if (km[pos] == w)
                    localMatch++;
            }

            bestLocalMatch = max(bestLocalMatch, localMatch);
        }

        totalMatch += bestLocalMatch;
    }

    return 100.0 * totalMatch / totalPositions;
}

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
        out[pos] = table[(pos + shift) % k]; // סיבוב שמאלה

    return out;
}

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
    return bestShift; // כמה הזזנו שמאלה
}

void AddCounts(vector<array<int, 4>> &a, const vector<array<int, 4>> &b)
{
    for (int pos = 0; pos < (int)a.size(); ++pos)
        for (int j = 0; j < 4; ++j)
            a[pos][j] += b[pos][j];
}
//------------------------>>>>>>>>>>>>>>>
vector<array<double, 4>> ComputeLocalDistributions(
    const string &dna, int W, double pseudo)
{
    int n = dna.size();
    vector<array<double, 4>> probs(n, {0, 0, 0, 0});
    array<int, 4> cnt{0, 0, 0, 0};

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

    int left = 0, right = min(n - 1, W);

    for (int i = left; i <= right; ++i)
        if (idx(dna[i]) >= 0)
            cnt[idx(dna[i])]++;

    for (int i = 0; i < n; ++i)
    {
        double sum = cnt[0] + cnt[1] + cnt[2] + cnt[3] + 4 * pseudo;
        for (int j = 0; j < 4; j++)
            probs[i][j] = (cnt[j] + pseudo) / sum;

        int newL = max(0, i + 1 - W);
        int newR = min(n - 1, i + 1 + W);

        while (left < newL)
            if (idx(dna[left]) >= 0)
                cnt[idx(dna[left++])]--;
        while (right < newR)
            if (idx(dna[++right]) >= 0)
                cnt[idx(dna[right])]++;
    }
    return probs;
}

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
    else
    {
        cerr << "RunPipeline: invalid mode = " << mode << endl;
        return {};
    }

    return segs;
}

vector<Segment> AnnotateSegmentsWithNull1(
    vector<Segment> &realSegs,
    const string &dna,
    int W, int mode)
{
    // 1️⃣ התפלגות מקומית
    auto probs = ComputeLocalDistributions(dna, W);

    // 2️⃣ יצירת DNA רנדומלי מקביל
    string rndDNA = GenerateRandomDNA(probs);

    // 3️⃣ הרצת אותו pipeline בדיוק
    vector<Segment> nullSegs = RunPipeline(rndDNA, mode);

    // 4️⃣ השוואה אחד-לאחד
    int N = min(realSegs.size(), nullSegs.size());

    for (int i = 0; i < N; ++i)
    {
        double pReal = realSegs[i].combinedPValue;
        double pNull = nullSegs[i].combinedPValue;
        realSegs[i].pNull = pNull;
        // יחס פשוט וברור
        double ratio = max(pNull, 1e-300) / max(pReal, 1e-300);

        realSegs[i].strengthScore = ratio;

        if (ratio > 10)
            realSegs[i].strength = "strong";
        else if (ratio > 1)
            realSegs[i].strength = "weak";
        else
            realSegs[i].strength = "noise";
    }
    return nullSegs;
}

double ComputeSegmentPValue(const string &segment, int K)
{
    if ((int)segment.size() < K)
        return 1.0;

    // 3. position counts
    int numKmers = segment.size() / K;
    auto counts = ComputePositionCounts(segment, K);

    // 4. position p-values
    auto posP = ComputePValuesFromCounts(counts, numKmers, K);

    // 5. combined p-value
    return CombinePValuesFisher(posP);
}

vector<string> SplitByLength(const string &dna)
{
    vector<string> chunks;
    for (int i = 0; i + L <= (int)dna.size(); i += L)
        chunks.push_back(dna.substr(i, L));
    return chunks;
}

void AnnotateSegmentsWithNull2(
    vector<Segment> &realSegs,
    const string &dna,
    int W)
{
    // 1️⃣ התפלגות מקומית
    auto probs = ComputeLocalDistributions(dna, W);

    // 2️⃣ יצירת DNA רנדומלי
    string rndDNA = GenerateRandomDNA(probs);

    // 3️⃣ חיתוך רנדומלי לפי L בלבד
    auto rndChunks = SplitByLength(rndDNA);

    int N = min((int)realSegs.size(), (int)rndChunks.size());

    for (int i = 0; i < N; ++i)
    {
        Segment &s = realSegs[i];

        int K = s.canonicalRepresentativeWord.size(); // K של האמיתי

        double pReal = s.combinedPValue;

        // p-value של הסיגמנט המקביל ברנדומלי
        double pNull = ComputeSegmentPValue(rndChunks[i], K);
        s.pNull = pNull;
        double ratio = max(pNull, 1e-300) / max(pReal, 1e-300);

        s.strengthScore = ratio;

        if (pReal < pNull * 0.1)
            s.strength = "strong";
        else if (pReal < pNull * 0.65)
            s.strength = "good";
        else
            s.strength = "noise";
    }
}

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
