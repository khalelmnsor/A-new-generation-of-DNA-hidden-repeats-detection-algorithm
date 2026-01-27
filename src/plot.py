import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator
from scipy.stats import ks_2samp

# ============================================================
# Environment setup
# ============================================================
# Set working directory to the script location to ensure
# reproducible relative file paths (CSV input / figure output)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# Output directory for all generated figures
OUT_DIR = "figures"
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# Helper functions
# ============================================================

def cdf(values):
    """
    Compute the empirical cumulative distribution function (CDF).

    Parameters
    ----------
    values : array-like
        Input values (e.g., -log10(P-values)).

    Returns
    -------
    x : ndarray
        Sorted values.
    y : ndarray
        Empirical CDF values in [0,1].
    """
    x = np.sort(values)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def beautify(ax, xM=None, xm=None, yM=None, ym=None):
    """
    Apply consistent grid and tick formatting to matplotlib axes.
    This is a visualization utility only (no analytical impact).
    """
    if xM: ax.xaxis.set_major_locator(MultipleLocator(xM))
    if xm: ax.xaxis.set_minor_locator(MultipleLocator(xm))
    if yM: ax.yaxis.set_major_locator(MultipleLocator(yM))
    if ym: ax.yaxis.set_minor_locator(MultipleLocator(ym))
    ax.grid(True, which="major", alpha=0.4)
    ax.grid(True, which="minor", alpha=0.15)


def logP(df):
    """
    Compute -log10(P-value) while discarding invalid zero entries.
    Used for statistical significance visualization.
    """
    return -np.log10(df[df["pvalue"] > 0]["pvalue"])


def agreement_ratio(df):
    """
    Fraction of segments where the representative word
    exactly matches the most frequent k-mer.

    Interpretation:
    Measures consensus quality between representative-word
    inference and raw frequency dominance.
    """
    return (df["representativeWord"] == df["frequentWord1"]).mean()


def agreement_ratio_by_count(df):
    """
    Fraction of segments where the representative word
    and the most frequent k-mer have identical occurrence counts.

    Interpretation:
    Stronger agreement criterion than string equality,
    focusing on quantitative dominance.
    """
    return (df["representativeCount"] == df["frequentCount1"]).mean()

# ============================================================
# Load segmentation results
# ============================================================
# Each CSV corresponds to a pipeline (AVG / P-value)
# and to real genomic data vs composition-preserving null model
avg_real = pd.read_csv("segments_real_avg.csv")
avg_null = pd.read_csv("segments_null_avg.csv")
pv_real  = pd.read_csv("segments_real_pvalue.csv")
pv_null  = pd.read_csv("segments_null_pvalue.csv")

# Remove invalid segments with undefined statistical scores
avg_real = avg_real[avg_real["pvalue"] > 0].copy()
avg_null = avg_null[avg_null["pvalue"] > 0].copy()
pv_real  = pv_real[pv_real["pvalue"] > 0].copy()
pv_null  = pv_null[pv_null["pvalue"] > 0].copy()

# Precompute -log10(P-value) for downstream analysis
avg_real["logP"] = -np.log10(avg_real["pvalue"])
pv_real["logP"]  = -np.log10(pv_real["pvalue"])
avg_null["logP"] = -np.log10(avg_null["pvalue"])
pv_null["logP"]  = -np.log10(pv_null["pvalue"])

# Segment strength labels (used later for stratified analysis)
classes = ["strong", "weak", "noise"]
colors  = {"strong": "red", "weak": "orange", "noise": "blue"}

# ============================================================
# Representative–Frequent agreement analysis (string equality)
# ============================================================
labels = ["AVG-REAL", "AVG-NULL", "PV-REAL", "PV-NULL"]

values = [
    agreement_ratio(avg_real),
    agreement_ratio(avg_null),
    agreement_ratio(pv_real),
    agreement_ratio(pv_null),
]

plt.figure(figsize=(7, 5))
plt.bar(labels, values, edgecolor="black")
plt.ylabel("Fraction of segments\n(repWord == frequentWord)")
plt.title("Representative–Frequent agreement")
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig01_rep_freq_agreement.png", dpi=300)

# Complementary view: false consensus rate
plt.figure(figsize=(7, 5))
plt.bar(labels, [1 - v for v in values],
        edgecolor="black", color="orange")
plt.ylabel("Fraction of segments\n(repWord ≠ frequentWord)")
plt.title("False consensus rate")
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig02_false_consensus.png", dpi=300)

# ============================================================
# Representative–Frequent agreement (by count equality)
# ============================================================
values = [
    agreement_ratio_by_count(avg_real),
    agreement_ratio_by_count(avg_null),
    agreement_ratio_by_count(pv_real),
    agreement_ratio_by_count(pv_null),
]


# ============================================================
# Agreement analysis (by COUNT equality)
# ============================================================
# This analysis evaluates whether the representative word
# captures the dominant repeat quantitatively, not only
# lexically. A match occurs when the representative word
# and the most frequent k-mer appear the same number of times
# within the segment.

plt.figure(figsize=(7, 5))
plt.bar(labels, values, edgecolor="black")

plt.ylabel(
    "Fraction of segments\n"
    "(repCount == frequentCount)"
)

plt.title(
    "Representative–Frequent agreement (by count)"
)

plt.ylim(0, 1)
plt.tight_layout()

plt.savefig(
    f"{OUT_DIR}/fig01_rep_freq_count_agreement.png",
    dpi=300
)

# ============================================================
# False consensus rate (by COUNT)
# ============================================================
# Complementary view: segments where the representative word
# does NOT match the dominant k-mer in terms of occurrence count.
# This highlights cases where a consensus pattern is inferred,
# but does not correspond to the most frequent repeat.

plt.figure(figsize=(7, 5))
plt.bar(
    labels,
    [1 - v for v in values],
    edgecolor="black",
    color="orange"
)

plt.ylabel(
    "Fraction of segments\n"
    "(repCount ≠ frequentCount)"
)

plt.title(
    "False consensus rate (by count)"
)

plt.ylim(0, 1)
plt.tight_layout()

plt.savefig(
    f"{OUT_DIR}/fig02_false_consensus_count.png",
    dpi=300
)

# ============================================================
# AVG distribution as percentage of segments (clean style)
# ============================================================

real = avg_real["AVG"].values
null = avg_null["AVG"].values

# Define bin step for AVG (%)
STEP = 1

# Determine range
x_min = int(np.floor(min(real.min(), null.min())))
x_max = int(np.ceil(max(real.max(), null.max())))

bins = np.arange(x_min, x_max + STEP, STEP)

# Convert counts to percentages
weights_real = np.ones_like(real) / len(real) * 100
weights_null = np.ones_like(null) / len(null) * 100

plt.figure(figsize=(12, 5))

plt.hist(real, bins=bins, weights=weights_real, alpha=0.7, label="Real")
plt.hist(null, bins=bins, weights=weights_null, alpha=0.7, label="Null")

# Labels and title
plt.xlabel("AVG match (%)")
plt.ylabel("Percentage of segments (%)")
plt.title("AVG distribution (AVG pipeline)")

# X-axis formatting
plt.xticks(np.arange(x_min, x_max + 1, 1), fontsize=8)
plt.xlim(x_min, x_max)

# Y-axis formatting (small fixed steps)
y_max = plt.ylim()[1]
y_max_round = int(np.ceil(y_max / 2) * 2)
plt.ylim(0, y_max_round)
plt.yticks(np.arange(0, y_max_round + 0.1, 2), fontsize=8)

# Legend outside
plt.legend(
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=True
)

plt.tight_layout()
plt.savefig(
    f"{OUT_DIR}/fig_avg_distribution_percent_clean.png",
    dpi=300,
    bbox_inches="tight"
)
# ============================================================
# AVG distribution as percentage of segments (P-VALUE pipeline)
# ============================================================

real = pv_real["AVG"].values
null = pv_null["AVG"].values

# Define bin step for AVG (%)
STEP = 1

# Determine range
x_min = int(np.floor(min(real.min(), null.min())))
x_max = int(np.ceil(max(real.max(), null.max())))

bins = np.arange(x_min, x_max + STEP, STEP)

# Convert counts to percentages
weights_real = np.ones_like(real) / len(real) * 100
weights_null = np.ones_like(null) / len(null) * 100

plt.figure(figsize=(12, 5))

plt.hist(real, bins=bins, weights=weights_real,
         alpha=0.7, label="Real")

plt.hist(null, bins=bins, weights=weights_null,
         alpha=0.7, label="Null")

# Labels and title
plt.xlabel("AVG match (%)")
plt.ylabel("Percentage of segments (%)")
plt.title("AVG distribution (P-value pipeline)")

# X-axis formatting
plt.xticks(np.arange(x_min, x_max + 1, 1), fontsize=8)
plt.xlim(x_min, x_max)

# Y-axis formatting (small fixed steps)
y_max = plt.ylim()[1]
y_max_round = int(np.ceil(y_max / 2) * 2)
plt.ylim(0, y_max_round)
plt.yticks(np.arange(0, y_max_round + 0.1, 2), fontsize=8)

# Legend outside
plt.legend(
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=True
)

plt.tight_layout()
plt.savefig(
    f"{OUT_DIR}/fig_avg_distribution_percent_pvalue_pipeline.png",
    dpi=300,
    bbox_inches="tight"
)




def word_length(df):
    """
    Returns the length (K) of the representative word
    for each detected segment.
    """
    return df["representativeWord"].str.len()

# ============================================================
# Representative word length distribution – AVG pipeline
# ============================================================
# This histogram shows the distribution of selected word lengths (K)
# in the AVG-based pipeline. Since K is chosen by maximizing average
# positional agreement, the pipeline may exhibit a bias toward longer
# words that accumulate smoother matches.

plt.figure(figsize=(7, 5))

plt.hist(
    word_length(avg_real),
    bins=np.arange(2, 10) - 0.5,
    alpha=0.7,
    label="Real"
)

plt.hist(
    word_length(avg_null),
    bins=np.arange(2, 10) - 0.5,
    alpha=0.7,
    label="Null"
)

plt.xlabel("Word length (K)")
plt.ylabel("Number of segments")
plt.title("Distribution of representative word length (AVG pipeline)")
plt.legend()
plt.tight_layout()

plt.savefig(
    f"{OUT_DIR}/fig_word_length_avg_pipeline.png",
    dpi=300
)

# ============================================================
# Representative word length distribution – P-value pipeline
# ============================================================
# This histogram shows the distribution of K selected by the
# statistical (P-value–based) pipeline. Here, K is chosen by
# minimizing the combined positional P-value, allowing the
# pipeline to adaptively favor word lengths that yield the
# strongest statistical signal rather than maximal smoothness.

plt.figure(figsize=(7, 5))

plt.hist(
    word_length(pv_real),
    bins=np.arange(2, 10) - 0.5,
    alpha=0.7,
    label="Real"
)

plt.hist(
    word_length(pv_null),
    bins=np.arange(2, 10) - 0.5,
    alpha=0.7,
    label="Null"
)

plt.xlabel("Word length (K)")
plt.ylabel("Number of segments")
plt.title("Distribution of representative word length (P-value pipeline)")
plt.legend()
plt.tight_layout()

plt.savefig(
    f"{OUT_DIR}/fig_word_length_pvalue_pipeline.png",
    dpi=300
)

# ============================================================
# Percentage histogram with clean axes and external legend
# ============================================================

real = avg_real["logP"].values
null = avg_null["logP"].values

CUT = 35

x_max = int(np.ceil(np.max(np.concatenate([real, null]))))

if x_max <= CUT:
    bins = np.arange(0, x_max + 1, 1)
else:
    bins = np.concatenate([
        np.arange(0, CUT + 1, 1),
        [x_max + 1]
    ])


# Convert counts to percentages
weights_real = np.ones_like(real) / len(real) * 100
weights_null = np.ones_like(null) / len(null) * 100

plt.figure(figsize=(12, 5))

plt.hist(real, bins=bins, weights=weights_real, alpha=0.7, label="Real")
plt.hist(null, bins=bins, weights=weights_null, alpha=0.7, label="Null")

# Labels and title
plt.xlabel("-log10(P-value)")
plt.ylabel("Percentage of segments (%)")
plt.title("Statistical significance distribution (AVG pipeline)")

# X-axis: unit steps
plt.xticks(np.arange(0, CUT + 1, 1), fontsize=8)
plt.xlim(0, CUT + 1)

# Y-axis: small fixed steps
y_max = plt.ylim()[1]
y_max_round = int(np.ceil(y_max / 2) * 2)  # round to nearest 2%
plt.ylim(0, y_max_round)
plt.yticks(np.arange(0, y_max_round + 0.1, 2), fontsize=8)

# Overflow bin label
plt.text(
    CUT + 0.5,
    y_max_round * 0.95,
    "≥35",
    ha="center",
    va="top",
    fontsize=10,
    fontweight="bold"
)

# Legend outside the plot
plt.legend(
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=True
)

plt.tight_layout()
plt.savefig(
    f"{OUT_DIR}/fig_logp_distribution_0_35_overflow_percent_clean.png",
    dpi=300,
    bbox_inches="tight"
)
# ============================================================
# P-VALUE pipeline: percentage histogram with overflow bin
# ============================================================

real = pv_real["logP"].values
null = pv_null["logP"].values

CUT = 35

# Determine maximum logP
x_max = int(np.ceil(np.max(np.concatenate([real, null]))))

# Define bins: unit bins up to 35 + overflow bin
if x_max <= CUT:
    bins = np.arange(0, x_max + 1, 1)
else:
    bins = np.concatenate([
        np.arange(0, CUT + 1, 1),
        [x_max + 1]
    ])


# Convert counts to percentages
weights_real = np.ones_like(real) / len(real) * 100
weights_null = np.ones_like(null) / len(null) * 100

plt.figure(figsize=(12, 5))

# Plot histograms
plt.hist(real, bins=bins, weights=weights_real, alpha=0.7, label="Real")
plt.hist(null, bins=bins, weights=weights_null, alpha=0.7, label="Null")

# Labels and title
plt.xlabel("-log10(P-value)")
plt.ylabel("Percentage of segments (%)")
plt.title("Statistical significance distribution (P-value pipeline)")

# X-axis: unit steps
plt.xticks(np.arange(0, CUT + 1, 1), fontsize=8)
plt.xlim(0, CUT + 1)

# Y-axis: fixed small steps
y_max = plt.ylim()[1]
y_max_round = int(np.ceil(y_max / 2) * 2)
plt.ylim(0, y_max_round)
plt.yticks(np.arange(0, y_max_round + 0.1, 2), fontsize=8)

# Overflow bin label
plt.text(
    CUT + 0.5,
    y_max_round * 0.95,
    "≥35",
    ha="center",
    va="top",
    fontsize=10,
    fontweight="bold"
)

# Legend outside the plot
plt.legend(
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=True
)

plt.tight_layout()
plt.savefig(
    f"{OUT_DIR}/fig_logp_distribution_pvalue_pipeline.png",
    dpi=300,
    bbox_inches="tight"
)

plt.figure(figsize=(8,6))

# CDF of -log10(P-value) for real segments only
# This comparison isolates the behavior of the two pipelines
# on real genomic data, without interference from null models.

# AVG pipeline (Real)
x1, y1 = cdf(avg_real["logP"])

# P-value pipeline (Real)
x2, y2 = cdf(pv_real["logP"])

plt.plot(x1, y1, label="AVG pipeline (Real)")
plt.plot(x2, y2, label="P-value pipeline (Real)")

plt.xlabel("-log10(P-value)")
plt.ylabel("CDF")
plt.title("CDF comparison: AVG vs P-value pipeline (Real segments)")
plt.legend()
plt.tight_layout()

plt.savefig(
    f"{OUT_DIR}/fig_cdf_avg_vs_pvalue_real.png",
    dpi=300
)



plt.figure(figsize=(8,6))

# AVG pipeline
plt.plot(*cdf(avg_real["logP"]), label="AVG-Real")
plt.plot(*cdf(avg_null["logP"]), "--", label="AVG-Null")

# P-value pipeline
plt.plot(*cdf(pv_real["logP"]), label="P-value-Real")
plt.plot(*cdf(pv_null["logP"]), "--", label="P-value-Null")

plt.xlabel("-log10(P-value)")
plt.ylabel("CDF")
plt.title("CDF of statistical significance (all pipelines)")
plt.legend()
plt.tight_layout()

plt.savefig(
    f"{OUT_DIR}/fig_cdf_all_pipelines.png",
    dpi=300
)



THRESHOLDS = [5, 10, 20]

def fraction_above(df, t):
    """
    Percentage of segments whose statistical significance
    exceeds a given -log10(P-value) threshold.
    """
    return (df["logP"] >= t).mean() * 100

print("\n=== Percentage of segments above significance thresholds ===")

for t in THRESHOLDS:
    print(f"\nThreshold: -log10(P) >= {t}")

    print(f"AVG-Real   : {fraction_above(avg_real, t):.2f}%")
    print(f"AVG-Null   : {fraction_above(avg_null, t):.2f}%")
    print(f"PVal-Real  : {fraction_above(pv_real, t):.2f}%")
    print(f"PVal-Null  : {fraction_above(pv_null, t):.2f}%")
plt.show()
