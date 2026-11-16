# Payload-integrity audit of the MOVER v2 arterial-pressure waveforms

Code accompanying the manuscript:

> **Payload integrity as a measurement-validity problem in routinely collected waveform data:
> an archive-wide audit of arterial-pressure usability and declared gain.**
> Lu J, Yang W, Li Y, Li Z, Zheng F. *BMC Medical Research Methodology* (submitted).

The audit asks whether the waveform archive of the
[Medical Informatics Operating Room Vitals and Events Repository (MOVER) v2](https://doi.org/10.24432/C5VS5G)
delivers what it declares. Its central result is that the **declared arterial-pressure gain is
correct** while **most of the stored payloads cannot carry the signal**: 94.6% of de-duplicated
`GE_ART` files are unusable (58.5% period-2 alternation, 36.1% noise-like), whereas the 2.27%
of paired cases whose payloads are readable decode to physiological pressures at the declared
gain.

## No data are included

The repository contains **code only**. MOVER is distributed under a data use agreement that
prohibits redistribution, so no patient-level or derived record-level data are included here.
Derived aggregate tables (the per-case gain fits, the six estimator arms, the archive census)
are available from the corresponding author on reasonable request.

## Contents

| File | Purpose |
|---|---|
| `ac1_census_20261030.py` | Payload classifier: single-pass lag-1 autocorrelation rule (classes S / A / N) applied to every de-duplicated `GE_ART` file of the archive-wide census (344,519 files). |
| `spotcheck_20261030.py` | Byte-level spot check: base64 → int16 on randomly sampled files, printing odd/even sample parity, so the alternation can be inspected without going through the classifier. |
| `pred_i_crosschannel_20261007.py` | Co-affliction test: within the same `IP-` files, are `GE_ART` and `GE_ECG` affected together (κ statistic)? |
| `pred_iii_negcontrol_20261007.py` | Negative control: the mutually exclusive `CB-` channel set (`ECG1` / `PLETH`) on the same corpus. |
| `fig1_positive_control_20261014.py` | Figure 1 — gain concordance, payload classes, physiological window. |
| `fig2_packing_path_20261014.py` | Figure 2 — co-affliction matrix, concordance, negative control. |
| `fig3_recoverability_20261014.py` | Figure 3 — six estimator arms, Bland–Altman for the released decoder, bin distributions. |
| `fig4_design_bound_20261014.py` | Figure 4 — anchor-bin-width design bound and the VitalDB second-corpus comparison. |

## Requirements

Python 3.9+ with `numpy` and `matplotlib`:

```bash
pip install -r requirements.txt
```

The census and the two channel tests were run on a machine with read access to the MOVER
distribution; the figure scripts read frozen aggregate tables and need no raw corpus (Figure 1
reads per-case gain tables produced from the corpus).

## Configuring paths

The scripts were written against one installation and their input paths are set at the top of
each file. Every path can be overridden by an environment variable, so no edit is needed:

| Variable | Meaning | Default (original installation) |
|---|---|---|
| `MOVER_WAVE_ROOT` | Root of the unpacked waveform distribution (`epic_wave_1_v2` …) | `/mnt/share/ljfang/exp/raw/wave_pilot_20260909` |
| `MOVER_WAVE_INDEX` | File index (CSV.GZ) listing every file, its channels and `has_wave` | `…/11_wave_gates/wave_file_index_20260911.csv.gz` |
| `MOVER_ECG_QC` | Per-file ECG quality table used by the channel test | `…/20_wp2_ecg/ecg_file_qc_full.csv` |
| `MOVER_ABPCONTROL` | Per-case gain tables behind Figure 1 | `…/run-20260913-abpcontrol` |
| `MOVER_PROJECT_ROOT` | Root containing the frozen analysis tables | parent of the repository |
| `MOVER_FIG_OUT` | Output directory for Figure 4 | `<root>/scratch/figs_20261014` |
| `MOVER_OUT` | Output directory for the census | `…/run-20261030-gate` |

Example:

```bash
export MOVER_WAVE_ROOT=/data/mover/wave_pilot
export MOVER_WAVE_INDEX=/data/mover/wave_file_index.csv.gz
python3 ac1_census_20261030.py --workers 16
```

The census de-duplicates on `(pid, suffix, fn_ts)` before classifying, because the distributed
index lists each file once per volume it appears in (34.4% of nominal hits are cross-volume
duplicates).

## Classifier

Each stored packet is 180 int16 codes (base64 in the XML). For one packet:

```python
ac1 = lag-1 autocorrelation of the 180 codes
```

A file's class is the median of its packets' `ac1`:

| Class | Rule | Reading |
|---|---|---|
| **A** | `ac1 < −0.5` | period-2 alternation (storage/packing fault) |
| **N** | `\|ac1\| ≤ 0.5` | noise-like |
| **S** | `ac1 > +0.5` | readable, smoothly sampled pulsatile series |

The classifier uses stored codes only. It does not read the numeric pressure channel, so it can
be applied to the whole archive without a reference measurement; the positive control in the
manuscript is a separate analysis that does use the numeric channel.

## Provenance

These are the exact versions used for the reported results. Classifier thresholds, sampling
seeds and pre-registered gate definitions are documented in the Supplementary Material
(Additional file 1, S1–S9) of the manuscript.

## License

MIT — see `LICENSE`. If you need a different license, contact the corresponding author.
