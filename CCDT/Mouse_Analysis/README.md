# Scripts Used for Analysis of Mouse Data

## Quick Install Guide

### Prerequisites
Before running the MATLAB scripts, make sure you have:
- **MATLAB R2022b** or later (recommended)  
- The following MATLAB toolboxes installed:  
  - Signal Processing Toolbox  
  - Statistics and Machine Learning Toolbox  
  - (Optional) Parallel Computing Toolbox for faster processing
  
### Installation Steps
These steps are estimated to take approximately 10 minutes.
1. **Clone or Download the Repository**
   ```bash
   git clone https://github.com/mhedlund/RNET_Public/CCDT/Mouse_Analysis/analysis_scripts.git
   cd analysis_scripts
   ```
   Or download the .zip file from GitHub and extract it.  
2. **Add Project to MATLAB Path**  
   Open MATLAB and run:  
   ```bash
   addpath(genpath(pwd)) % pwd prints current working directory
   savepath
   ```  
   This makes all functions and scripts accessible.

---

## Data Organization

### Expected directory layout
Set `basePath` to the root directory containing the raw data:
```
<basePath>/
  └─ <protocol>/
       └─ <mouse_id>/
            ├─ experiment_000/
            │   ├─ experiment_000.mat          % photometry (signals)
            │   └─ experiment_000_bpod.mat     % Bpod/behavior (events & RTs)
            ├─ experiment_001/
            └─ ...
```
> The main script constructs paths like: `fullfile(basePath, protocol, mouse, sprintf('experiment_%03d', i))`.

### Required contents of each experiment
- **Photometry file**: `experiment_###.mat` containing (per fiber):
  - `sig` — fluorescence signal (e.g., 470 nm), size `[nFibers × T]` (or `[T × nFibers]`; transposed internally)
  - `ref` — isosbestic/control (e.g., 405 nm), same shape as `sig`
- **Behavior file**: `experiment_###_bpod.mat` containing trial timing and behavior, including fields used to derive:
  - Trial start times, first‑tone onset times, per‑trial reaction times (RT), and inter‑trial intervals (ITI)
  - Code adjusts Bpod vs photometry clocks and computes aligned times

> **Sampling**: The code sets `fs = 1/20`, treating this as the **time step** (0.05 s). PSTHs span −10 s to +10 s around the first tone (401 samples).

---

## File Descriptions (order of use)

### `processDataForExperiment_zscored.m` (function)
**Purpose.** Preprocess one experiment and return trial‑wise matrices/vectors.

**Key steps:**  
1. **Load files**: `<filename>.mat` (photometry) and `<filename>_bpod.mat` (behavior).  
2. **Motion/isosbestic regression**: For each fiber, fit `ref` → `sig` with a 1st‑order polynomial, scale the control, and subtract (`sig_norm = sig − scaled_ref + offset`).  
3. **Z‑scoring**: Convert to `sig_normZ` (z‑scored), using pre‑tone baseline.  
4. **Peri‑event alignment**: Build `psthMatrix` for trials centered on the **first tone** with a window of −10 s…+10 s.  
5. **Trial classification**: Mark **successful** vs **unsuccessful** trials using RT thresholds and NaNs. Create:  
   - `psthSuccessMatrix` — successes aligned to the first tone  
   - `psthWithoutSecondMatrix` — trials without a valid second tone (i.e. unsuccessful trials)  
6. **Behavioral vectors**: Return `RTSuccessful`, `RTUnsuccessful`, `RTToneSuccessful`, `RTToneUnsuccessful`, `ITISuccessful`, `ITIUnsuccessful`.  
7. **Timebase**: Return `psthTime` (−10 s:0.05 s:+10 s).  

**Signature (inputs/outputs).**
```matlab
[psthMatrix, ITISuccessful, ITIUnsuccessful, ...
 RTToneSuccessful, RTToneUnsuccessful, ...
 RTSuccessful, RTUnsuccessful, ...
 psthSuccessMatrix, psthWithoutSecondMatrix, psthTime] = ...
    processDataForExperiment_zscored(filename);
```
> If `<filename>_bpod.mat` is missing, the function returns empty arrays and skips the experiment.

---

### `MFPanalysis.m` (script)
**Purpose.** Batch all experiments across mice/protocols, aggregate matrices, and generate manuscript figures.

**User‑set parameters (top of file).**  
- `basePath` — root path to your data  
- `protocols` — cell array of protocol folder names for either the test or training data (e.g., `{'test'}`)  
- `mice` — cell array of mouse IDs (e.g., `{'m1','m2','m3','m4','m5','m6'}`)  
- `first_experiment_count`, `last_experiment_count` — numeric range for `experiment_###`  
- `numfibers` — number of photometry input channels (keep at `2`)  

**What the script does.**  
1. Loops over `protocols × mice × experiments` and calls `processDataForExperiment_zscored(experimentPath)`.  
2. **Aggregates** across experiments into combined structures (per mouse and pooled): PSTHs, RT vectors, and ITIs for **successful** and **unsuccessful** trials.  
3. Defines **time windows** used in the manuscript:  
   - **Preparatory**: −4 s to 0 s (`timeWindow1`)  
   - **Anticipatory/response**: 0 s to 2 s (`timeWindow2`)  
4. Builds **adaptive RT bins** for successes to ensure sufficient counts:  
   - Start with 20 ms bins from 0→0.5 s; **merge forward** until each bin has ≥ 100 trials.  
5. Produces the figure panels listed below.  

---

## Reproducing the Figures

> These figures are generated directly by `MFPanalysis.m` after aggregation.

### Figure 4
- **4C — Successful RT distribution.** Histogram of `RTToneSuccessful` (bin width 10 ms).
- **4D — Example trials.** Overlay PSTHs for **2nd fastest** vs **2nd slowest** successful trials (Fiber 1), with vertical markers at 0, 2, 2.5 s.
- **4E — Gradient PSTH (adaptive RT bins).** For each adaptive RT bin (≥ 100 trials), plot mean ± SEM PSTH (per fiber) with a jet gradient.
- **4F — RT vs z-scored Ca2+ signal relationship (0–2 s).** Scatter of **bin‑center RT** vs **mean z-scored Ca2+ signal** in 0–2 s with linear fit, R², and slope p‑value (CI bands).
- **4G — Fastest‑fraction gradient.** PSTHs for cumulative fastest segments (1/10 → 1/2 of trials), each with mean ± SEM.
- **4H — RT vs z-scored Ca2+ signal relationship for fastest trials.** Mean z-scored Ca2+ signal in −4–0 s vs mean RT for the same fastest segments, with quadratic fit and CI.
- **4I — Fast vs slow overlay (successes).** PSTHs for **fast** (0–20 ms) vs **slow** (20–500 ms) trials (mean ± SEM).
- **4K — 0–100 ms third‑split overlay.** Using fastest vs slowest **thirds**, overlay PSTHs but constrain the RT range to 0–100 ms.
- **4J — Δ mean z-scored Ca2+ signal barplot.** Difference (Fast–Slow) in mean z-scored Ca2+ signal for **−4–0 s** and **0–2 s** (with two‑sample t‑tests).
- **4L — Δ mean z-scored Ca2+ signal within 0–100 ms** Difference (fastest third vs slowest third, RT<100 ms) in mean z-scored Ca2+ signal for **−4–0 s** and **0–2 s** (with two‑sample t‑tests).

### Supplemental S7–S9
- **S7 — Additional stratifications.** PSTH overlays and summary stats for alternative windows and groupings
- **S8 — Δ mean z-scored Ca2+ signal between MD soma and MD→dmPFC by RT bin** PSTH overlays and summary stats for −4–0 s and 0–2 s.
- **S9 — Fastest / Slower / Near‑Miss / False‑Start.** PSTHs and barplots (means ± SEM) for four groups:
  - **Fastest** (≤ 20 ms)
  - **Slower** (20–500 ms)
  - **Near‑misses** (`RTToneUnsuccessful` > 1.98 s)
  - **False starts** (`RTToneUnsuccessful` < 0.1 s)
  - Includes pairwise t‑tests between groups for both windows.

> All plots are generated for **both fibers** (`numfibers = 2`); vertical lines mark 0s, 2s, and 2.5 s where applicable, with 0s being the start of the first tone.

---

## How to Run (step‑by‑step)

1. **Open** `MFPanalysis.m` and set:
   ```matlab
   basePath = '.../mfp_analysis/data';
   protocols = {'test'};
   mice = {'m1','m2', 'm3', 'm4', 'm5', 'm6'};
   first_experiment_count = 0;
   last_experiment_count  = 34;
   numfibers = 2;
   ```
2. **Run** `MFPanalysis.m`. Figures for **Figure 4** and **Supplemental S7–S9** will be created in MATLAB.
3. **(Optional) Save figures**: A commented block at the end shows how to save all open figures as `.svg` (uncomment and set `outputDir`).

---

## Outputs

- **Aggregated variables (workspace)**: combined PSTH matrices per fiber for successful/unsuccessful trials; vectors of RT (tone‑aligned) and ITI for each group.
- **Figures (MATLAB windows)**: Panels corresponding to the manuscript (Figure 4 & Supplemental S7–S9).

---

## Notes

- **First‑trial exclusion**: The first trial is removed from `psthMatrix` (warm‑up/edge effects).
