#!/usr/bin/env python3
import pandas as pd
import glob
from pathlib import Path

def check_gaps(file_path):
    print(f"\nChecking {file_path}...")
    # Files are space-delimited
    try:
        # Using engine='c' for speed on large files
        df = pd.read_csv(file_path, sep=r'\s+', engine='c')
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    if 'time' not in df.columns:
        print(f"Error: 'time' column not found in {file_path}")
        return

    # time is in decimal hours. 8Hz = 0.125s = 0.125 / 3600 hours
    expected_step = 0.125 / 3600
    # Use a slightly larger threshold to account for floating point precision
    threshold = expected_step * 1.05 

    diffs = df['time'].diff()
    
    # Find where the difference is significantly larger than expected
    gaps = diffs[diffs > threshold]

    if not gaps.empty:
        print(f"Found {len(gaps)} total gaps in {file_path}")
        
        # Gaps > 1s
        long_gaps = gaps[gaps * 3600 > 1.0]
        if not long_gaps.empty:
            print(f"  Found {len(long_gaps)} gaps > 1s:")
            for idx in long_gaps.index[:50]: # Show up to 50 long gaps
                t_prev = df.loc[idx-1, 'time']
                t_curr = df.loc[idx, 'time']
                gap_s = (t_curr - t_prev) * 3600
                print(f"    Long gap at index {idx:10d}: {t_prev:12.8f} -> {t_curr:12.8f} ({gap_s:8.3f} s)")
            
            if len(long_gaps) > 50:
                print(f"    ... and {len(long_gaps) - 50} more gaps > 1s.")
        
        print(f"Max gap: {gaps.max() * 3600:.3f} s")
    else:
        print(f"No gaps found in {file_path} (checked {len(df)} records)")

def main():
    # Check files in the matlab directory
    files = sorted(glob.glob("matlab/lecs_ml_data_*.dat"))
    if not files:
        # Also check in the processed data directory as a fallback
        files = sorted(glob.glob("data/processed/matlab_input/lecs_ml_data_*.dat"))
        
    if not files:
        print("No files found matching lecs_ml_data_*.dat in matlab/ or data/processed/matlab_input/")
        return

    for f in files:
        check_gaps(f)

if __name__ == "__main__":
    main()
