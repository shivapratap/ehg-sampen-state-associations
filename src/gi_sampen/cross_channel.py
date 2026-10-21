#!/usr/bin/env python3
"""Stage F1 -- cross-channel replication of GI-SD and GI-Total on EHG10-12.

Only the two features that survived the resolution ablation are carried
forward. Kurtosis and skewness were dropped there (resolution-sensitive,
sign-unstable) and are not extracted here -- there is no reason to spend
compute replicating a feature that already failed on EHG9.

EHG9 is not re-extracted: its GI-SD / GI-Total values already exist in
work/gi_features_60s_EHG9.csv from Stage A2, at the same resolution
(decimals=3) used here. This stage produces the matching EHG10/11/12 files.

Each channel gets a short invariance spot-check (30 windows x 2 scales)
before its features are trusted -- cheap (~15s), and it catches a
channel-specific bug (e.g. a bad electrode read) that a global correctness
test on EHG9 alone would not.
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from . import paths as _common
from . import gi_profile

CHANNELS = ("EHG10", "EHG11", "EHG12")
SPOT_CHECK_N = 30
SPOT_CHECK_SCALES = (0.5, 2.0)
SPOT_CHECK_ATOL = 1e-10
SPOT_CHECK_RTOL = 1e-8


def spot_check_invariance(channel: str, frame: pd.DataFrame) -> bool:
    sample = frame.sample(n=min(SPOT_CHECK_N, len(frame)), random_state=1) \
        .sort_values(["recording_id", "window_start_s"])
    ok = True
    for row, signal in _common.iter_windows(sample, progress=False, channel=channel):
        base = gi_profile.gi_summaries_rounded(signal)
        for scale in SPOT_CHECK_SCALES:
            scaled = gi_profile.gi_summaries_rounded(scale * signal)
            for key in ("gi_sd", "gi_total"):
                if not np.isclose(scaled[key], base[key],
                                  atol=SPOT_CHECK_ATOL, rtol=SPOT_CHECK_RTOL):
                    print(f"  INVARIANCE SPOT-CHECK FAILED: {channel} "
                          f"{row.window_id} {key} scale={scale} "
                          f"base={base[key]} scaled={scaled[key]}")
                    ok = False
    return ok


def extract_channel(channel: str, quiet: bool) -> pd.DataFrame:
    frame = _common.load_analysis_windows(channel=channel)

    print(f"[{channel}] invariance spot-check ({SPOT_CHECK_N} windows, "
          f"scales {SPOT_CHECK_SCALES})...", flush=True)
    if not spot_check_invariance(channel, frame):
        raise RuntimeError(
            f"{channel}: invariance spot-check failed. Do not trust this "
            "channel's features -- investigate before proceeding.")
    print(f"[{channel}] spot-check passed.", flush=True)

    started = time.perf_counter()
    rows = []
    for row, signal in _common.iter_windows(frame, progress=not quiet, channel=channel):
        features = gi_profile.gi_summaries_rounded(signal)
        rows.append({
            "window_id": row.window_id, "recording_id": row.recording_id,
            "patient_id": row.patient_id, "window_class": row.window_class,
            "channel": channel,
            "gi_sd": features["gi_sd"], "gi_total": features["gi_total"],
        })
    elapsed = time.perf_counter() - started
    print(f"[{channel}] {len(rows)} windows in {elapsed/60:.1f} min", flush=True)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channels", nargs="+", default=list(CHANNELS),
                        choices=list(CHANNELS))
    parser.add_argument("--probe", action="store_true",
                        help="time 100 windows on the first channel and stop")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    if args.probe:
        channel = args.channels[0]
        frame = _common.load_analysis_windows(limit=100, channel=channel)
        started = time.perf_counter()
        for _row, signal in _common.iter_windows(frame, progress=not args.quiet,
                                                  channel=channel):
            gi_profile.gi_summaries_rounded(signal)
        elapsed = time.perf_counter() - started
        total = len(_common.load_analysis_windows(channel=channel))
        per_window = elapsed / 100
        print(f"PROBE ({channel}): 100 windows in {elapsed:.1f}s "
              f"({per_window*1000:.0f} ms/window)")
        print(f"  projected per channel: {per_window*total/60:.1f} min")
        print(f"  projected for {len(args.channels)} channels: "
              f"{per_window*total*len(args.channels)/60:.1f} min")
        return 0

    all_rows = [extract_channel(channel, args.quiet) for channel in args.channels]
    combined = pd.concat(all_rows, ignore_index=True)
    _common.write_csv(combined,
                      _common.WORK_DIR / "cross_channel_gi_features.csv",
                      label="cross-channel GI features")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
