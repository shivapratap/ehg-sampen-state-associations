# Data

Nothing in this repository is raw patient data. There are two kinds of input.

## 1. Raw recordings (not included): PhysioNet

**Icelandic 16-electrode electrohysterogram database** (Alexandersson et al., *Scientific Data* 2015), open access.

- Landing page and download: <https://physionet.org/content/ehgdb/> (DOI 10.13026/c2159j; the files section shows the current version, 1.0.0 at the time of writing). Check the landing page for the licence and the exact citation to use.
- Download everything in one folder, for example:
  ```bash
  wget -r -N -c -np -nd -P data/raw https://physionet.org/files/ehgdb/1.0.0/
  ```
  (`-nd` puts all files directly in `data/raw/`, which is the layout this code expects: the WFDB `.dat/.hea/.atr` files together with `RECORDS`, `ANNOTATORS` and `SHA256SUMS.txt`.)
- Size: a few gigabytes. `data/raw/` is git-ignored; only a `.gitkeep` is tracked.
- Point the code at it with `export EHG_RAW_DIR=/path/to/folder` (default: `data/raw/`).
- `python src/data.py` verifies the files against the SHA-256 list that PhysioNet distributes (this repository mints no checksums of its own).

You only need the raw data for the stages marked **raw** in `python scripts/run_all.py --list`.

## 2. Derived tables (included): `data/derived/ehg_derived_inputs.zip`

Everything the statistical stages start from, produced from the raw data by the raw stages (window manifest, annotation tables, the 34-feature tables for EHG9–EHG12, QC tables, the event-unit/model tables, and the GI feature tables). ~40 MB; contents and SHA-256 of every file in `data/derived/MANIFEST.sha256`.

```bash
python scripts/restore_derived.py      # unpacks into work/ and verifies checksums
```

These tables contain window-level feature values and pseudonymous recording identifiers such as `ice002_p_1of3` (as used in the public database). They contain no signal samples. The raw signals cannot be rebuilt from them.

Not in the zip: the 60 audited waveform windows used by the gain audit (`results/manuscript/audits/waveforms/`, 1 MB) are committed under `results/` instead, because the gain audit and the independent grid check run from them without raw data.
