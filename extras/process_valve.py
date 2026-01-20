"""
Script to read and combine CSV files from a directory.
Extracts lines starting with "Rebooted" before processing CSV data.
"""

import pandas as pd
from pathlib import Path
from typing import List, Tuple
from datetime import datetime


def extract_reboot_lines(file_path: Path) -> Tuple[List[str], List[str]]:
    """
    Separate reboot lines from CSV data lines.
    
    Args:
        file_path: Path to the file to process
        
    Returns:
        Tuple of (reboot_lines, csv_lines)
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    reboot_lines = [line.strip() for line in lines if line.startswith('Rebooted')]
    csv_lines = [line for line in lines if not line.startswith('Rebooted')]
    
    return reboot_lines, csv_lines


def process_file(file_path: Path) -> pd.DataFrame:
    """
    Process a single CSV file, extracting reboot info.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        DataFrame with data and metadata columns
    """
    reboot_lines, csv_lines = extract_reboot_lines(file_path)
    
    # Read CSV data from the filtered lines
    from io import StringIO
    csv_data = pd.read_csv(StringIO(''.join(csv_lines)))
    
    # Add metadata columns
    csv_data['source_file'] = file_path.name
    csv_data['reboot_info'] = '; '.join(reboot_lines) if reboot_lines else ''
    
    return csv_data


def get_csv_files(directory: str) -> List[Path]:
    """
    Get list of all CSV files in a directory.
    
    Args:
        directory: Path to the directory
        
    Returns:
        List of Path objects for CSV files
    """
    return list(Path(directory).glob('*.csv'))


def combine_csv_files(csv_dir: str, output_path: str = None) -> pd.DataFrame:
    """
    Combine all CSV files from a directory into a single DataFrame.
    
    Args:
        csv_dir: Directory containing CSV files
        output_path: Optional path to save combined data
        
    Returns:
        Combined DataFrame
    """
    csv_files = get_csv_files(csv_dir)
    
    if not csv_files:
        print(f"No CSV files found in {csv_dir}")
        return pd.DataFrame()
    
    # Process and combine all files
    dataframes = [process_file(file) for file in csv_files]
    combined_data = pd.concat(dataframes, ignore_index=True)
    
    # Print summary
    print(f"Successfully combined {len(csv_files)} CSV files")
    print(f"Total rows: {len(combined_data)}")
    print(f"Total columns: {len(combined_data.columns)}")
    
    # Save if output path provided
    if output_path:
        combined_data.to_csv(output_path, index=False)
        print(f"Saved to {output_path}")
    
    return combined_data


def extract_reboot_times(file_path: Path) -> pd.DataFrame:
    """
    Extract reboot timestamps from a file.

    Args:
        file_path: Path to the file to process

    Returns:
        DataFrame with reboot timestamps and source file
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Extract lines starting with "Rebooted at"
    reboot_lines = [line.strip() for line in lines if line.startswith('Rebooted at')]

    # Parse timestamps from reboot lines
    reboot_times = []
    for line in reboot_lines:
        try:
            timestamp_str = line.split('at ')[1]
            # Remove 'Z' if present and parse the timestamp
            timestamp = datetime.strptime(timestamp_str.replace('Z', ''), "%Y-%m-%dT%H:%M:%S")
            reboot_times.append({'timestamp': timestamp, 'source_file': file_path.name})
        except (IndexError, ValueError):
            print(f"Failed to parse timestamp in line: {line}")

    return pd.DataFrame(reboot_times)


def combine_reboot_times(csv_dir: str) -> pd.DataFrame:
    """
    Combine reboot times from all files in a directory.

    Args:
        csv_dir: Directory containing files

    Returns:
        DataFrame with all reboot timestamps
    """
    csv_files = get_csv_files(csv_dir)

    if not csv_files:
        print(f"No CSV files found in {csv_dir}")
        return pd.DataFrame()

    # Extract and combine reboot times
    reboot_dataframes = [extract_reboot_times(file) for file in csv_files]
    combined_reboot_times = pd.concat(reboot_dataframes, ignore_index=True)

    return combined_reboot_times


if __name__ == "__main__":
    # Configuration
    CSV_DIR = "data/raw/GEMS_VALVE/"
    OUTPUT_PATH = "data/processed/valve_combined.csv"
    REBOOT_OUTPUT_PATH = "data/processed/reboot_times.csv"
    
    # Process files
    df = combine_csv_files(CSV_DIR, OUTPUT_PATH)
    
    # Extract and save reboot times
    reboot_df = combine_reboot_times(CSV_DIR)
    reboot_df.to_csv(REBOOT_OUTPUT_PATH, index=False)
    print(f"Reboot times saved to {REBOOT_OUTPUT_PATH}")
    
    # Display first few rows
    print("\nFirst few rows:")
    print(df.head())
    print("\nData info:")
    print(df.info())
    
    # Display first few reboot times
    print("\nFirst few reboot times:")
    print(reboot_df.head())
