#!/usr/bin/env python3
"""
Check for gaps and inconsistencies in ADV timestamp data.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Load the parquet file
file_path = Path("data/processed/lander/gems_adv.parquet")
print(f"Reading {file_path}...")

try:
    df = pd.read_parquet(file_path)
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
    exit(1)

print(f"\nDataset shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}\n")

# Identify timestamp column
timestamp_cols = [col for col in df.columns if 'time' in col.lower() or 'date' in col.lower()]
print(f"Potential timestamp columns: {timestamp_cols}")

if not timestamp_cols:
    print("Error: No timestamp column found")
    exit(1)

ts_col = timestamp_cols[0]
print(f"\nUsing column: '{ts_col}'")
print(f"Data type: {df[ts_col].dtype}")

# Convert to datetime if needed
if not pd.api.types.is_datetime64_any_dtype(df[ts_col]):
    print(f"Converting {ts_col} to datetime...")
    df[ts_col] = pd.to_datetime(df[ts_col])

# Filter to date range
start_date = pd.Timestamp("2025-07-16", tz="UTC")
end_date = pd.Timestamp("2025-07-17", tz="UTC")
df = df[(df[ts_col] >= start_date) & (df[ts_col] < end_date)].reset_index(drop=True)
print(f"Filtered to date range: {start_date.date()} to {end_date.date()}")
print(f"Filtered dataset shape: {df.shape}\n")

# Analysis
print("\n" + "="*60)
print("TIMESTAMP ANALYSIS")
print("="*60)

# Basic statistics
print(f"\nTime range: {df[ts_col].min()} to {df[ts_col].max()}")
print(f"Total duration: {df[ts_col].max() - df[ts_col].min()}")
print(f"Total records: {len(df)}")

# Check for NaN values
nan_count = df[ts_col].isna().sum()
print(f"\nMissing values: {nan_count}")

# Check for duplicates
dup_count = df[ts_col].duplicated().sum()
print(f"Duplicate timestamps: {dup_count}")
if dup_count > 0:
    print(f"  Duplicate values:\n{df[df[ts_col].duplicated(keep=False)].sort_values(ts_col)}")

# Check if monotonically increasing
is_monotonic = df[ts_col].is_monotonic_increasing
print(f"Monotonically increasing: {is_monotonic}")

if not is_monotonic:
    # Find non-monotonic points
    non_mono = (df[ts_col].diff() <= pd.Timedelta(0)).sum() - 1  # -1 for the first NaN
    print(f"  Non-monotonic transitions: {non_mono}")

# Calculate time differences
time_diffs = df[ts_col].diff()
print(f"\nTime differences (between consecutive records):")
print(f"  Min: {time_diffs.min()}")
print(f"  Max: {time_diffs.max()}")
print(f"  Mean: {time_diffs.mean()}")
print(f"  Median: {time_diffs.median()}")
print(f"  Std Dev: {time_diffs.std()}")

# Find gaps (assuming expected frequency is the median interval)
expected_freq = time_diffs.median()
print(f"\nExpected interval (median): {expected_freq}")

gaps = time_diffs[time_diffs > expected_freq * 1.5]
if len(gaps) > 0:
    print(f"\nGaps detected (>{expected_freq * 1.5}):")
    print(f"  Number of gaps: {len(gaps)}")
    print(f"  Largest gap: {gaps.max()}")
    
    # Check for gaps > 1s
    large_gaps = time_diffs[time_diffs > pd.Timedelta(seconds=1)]
    print(f"  Gaps > 1s: {len(large_gaps)}")
    if len(large_gaps) > 0:
        print(f"  Largest gap > 1s: {large_gaps.max()}")
        print(f"\n  All gaps > 1s (start -> end: duration):")
        for idx in large_gaps.index:
            gap_size = time_diffs[idx]
            start_time = df[ts_col].iloc[idx - 1]
            end_time = df[ts_col].iloc[idx]
            print(f"    {start_time} -> {end_time}: {gap_size}")

    print(f"\n  Gap details (first 10 gaps > {expected_freq * 1.5}):")
    gap_indices = gaps.index[:10]
    for idx in gap_indices:
        gap_size = time_diffs[idx]
        start_time = df[ts_col].iloc[idx - 1]
        end_time = df[ts_col].iloc[idx]
        print(f"    {start_time} -> {end_time}: {gap_size}")
else:
    print(f"\nNo significant gaps detected!")

# Show gap details with full rows
if len(gaps) > 0:
    # Prioritize showing gaps > 1s if they exist
    large_gaps = time_diffs[time_diffs > pd.Timedelta(seconds=1)]
    if len(large_gaps) > 0:
        print(f"\nDetailed information for gaps > 1s (up to 5):")
        display_indices = large_gaps.index[:5]
    else:
        print(f"\nDetailed gap information (first 5 gaps):")
        display_indices = gaps.index[:5]

    for i, idx in enumerate(display_indices, 1):
        gap_size = time_diffs[idx]
        print(f"\n  Gap {i}: {gap_size}")
        print(f"  Last record before gap (index {idx-1}):")
        print(f"    {df.iloc[idx-1]}")
        print(f"  First record after gap (index {idx}):")
        print(f"    {df.iloc[idx]}")

print("\n" + "="*60)
