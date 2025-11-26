import argparse
import pandas as pd
import subprocess
import sys
from pathlib import Path
from glob import glob


def calculate_ipsae(
    pae_file_path: str,
    structure_file_path: str,
    pae_cutoff: float = 15.0,
    dist_cutoff: float = 15.0,
) -> float:
    """
    Calculate ipSAE score for a protein-protein interaction.

    Parameters:
    -----------
    pae_file_path : str
        Path to the PAE file (NPZ for Boltz)
    structure_file_path : str
        Path to the structure file (mmCIF)
    pae_cutoff : float
        Cutoff value for PAE in score calculations
    dist_cutoff : float
        Cutoff value for distance in score calculations

    Returns:
    --------
    float
        The ipSAE score (lower is better)
    """
    # Try different python commands and workdirs
    # Docker: python + /workdir, Windows: py + .
    success = False
    for python_cmd in ["python", "py"]:
        for workdir in ["/workdir", "."]:
            ipsae_script = Path(workdir) / "ipsae.py"
            if not ipsae_script.exists():
                continue
            try:
                command = [
                    python_cmd,
                    str(ipsae_script),
                    pae_file_path,
                    structure_file_path,
                    str(pae_cutoff),
                    str(dist_cutoff),
                ]
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode == 0:
                    success = True
                    break
                else:
                    print(f"ipsae.py failed with: {result.stderr}")
            except Exception as e:
                print(f"Error running ipsae.py with command: {command}: {e}")
                continue
        if success:
            break
    
    if not success:
        print(f"Failed to run ipsae.py for {structure_file_path}")
        return float('inf')

    # Read results from the output file
    results_path = structure_file_path.replace('.cif', f'_{int(pae_cutoff)}_{int(dist_cutoff)}.txt')
    
    try:
        df = pd.read_csv(results_path)
        results = {}

        for i, row in df[df.Type == "max"].iterrows():
            chainpair = f"{row['Chn1']}-{row['Chn2']}"
            results[chainpair] = {"max": {col: row[col] for col in df.columns[5:-1]}}
            
            mask = (df['Chn1'] == row['Chn1']) & (df['Chn2'] == row['Chn2']) & (df['Type'] != "max")
            min_vals = df[mask][df.columns[5:-1]].min()
            results[chainpair]["min"] = min_vals.to_dict()

        # Try common chain pairs
        for chain_pair in ['A-B', 'B-A', 'A-C', 'C-A']:
            if chain_pair in results:
                return results[chain_pair]['min']['ipSAE']
        
        # Fall back to first available pair
        if results:
            first_pair = list(results.keys())[0]
            return results[first_pair]['min']['ipSAE']
        
        return float('inf')
    
    except Exception as e:
        print(f"Error reading results for {structure_file_path}: {e}")
        return float('inf')


def write_fasta(designs: list, output_path: str):
    """
    Write designs to a FASTA file.

    Parameters:
    -----------
    designs : list
        List of dicts with 'id', 'sequence', 'ipsae' keys
    output_path : str
        Path to output FASTA file
    """
    with open(output_path, 'w') as f:
        for design in designs:
            header = f">{design['id']}|ipSAE={design['ipsae']:.4f}|rank={design['ipsae_rank']}"
            f.write(f"{header}\n")
            # Write sequence in 80-character lines
            seq = design['sequence']
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def process_final_designs(
    final_designs_dir: str,
    pae_cutoff: float = 15.0,
    dist_cutoff: float = 15.0,
    output_fasta: str = None,
):
    """
    Process final ranked designs, calculate ipSAE, and output ranked FASTA.

    Parameters:
    -----------
    final_designs_dir : str
        Path to the final_ranked_designs directory
    pae_cutoff : float
        Cutoff value for PAE in score calculations
    dist_cutoff : float
        Cutoff value for distance in score calculations
    output_fasta : str
        Path to output FASTA file (default: <final_designs_dir>/designs_ranked_by_ipsae.fasta)
    """
    final_dir = Path(final_designs_dir)
    
    # Find the final designs subdirectory (final_X_designs)
    final_subdir = None
    for subdir in final_dir.iterdir():
        if subdir.is_dir() and subdir.name.startswith("final_") and subdir.name.endswith("_designs"):
            final_subdir = subdir
            break
    
    if final_subdir is None:
        print(f"Error: Could not find final_*_designs subdirectory in {final_dir}")
        sys.exit(1)
    
    print(f"Processing designs from: {final_subdir}")
    
    # Find metrics CSV
    metrics_csvs = list(final_dir.glob("final_designs_metrics_*.csv"))
    if not metrics_csvs:
        # Try all_designs_metrics.csv as fallback
        metrics_csvs = list(final_dir.glob("all_designs_metrics.csv"))
    
    if not metrics_csvs:
        print(f"Error: Could not find metrics CSV in {final_dir}")
        sys.exit(1)
    
    metrics_csv = metrics_csvs[0]
    print(f"Loading metrics from: {metrics_csv}")
    
    df_metrics = pd.read_csv(metrics_csv)
    
    # Build a mapping from filename to sequence
    filename_to_seq = {}
    filename_to_chain_seq = {}
    for _, row in df_metrics.iterrows():
        filename_to_seq[row['file_name']] = row['designed_sequence']
        filename_to_chain_seq[row['file_name']] = row['designed_chain_sequence']
    
    # Find CIF and PAE files
    cif_files = sorted(final_subdir.glob("*.cif"))
    pae_dir = final_subdir / "pae"
    
    if not pae_dir.exists():
        print(f"Error: PAE directory not found at {pae_dir}")
        print("Make sure you've run the pipeline with PAE output enabled.")
        sys.exit(1)
    
    print(f"Found {len(cif_files)} CIF files")
    print(f"Calculating ipSAE scores...")
    
    designs = []
    
    for cif_path in cif_files:
        cif_name = cif_path.name
        
        # Extract original filename (remove rank prefix)
        # Format: rank00_original_name.cif
        parts = cif_name.split('_', 1)
        if len(parts) > 1:
            original_name = parts[1]  # original_name.cif
        else:
            original_name = cif_name
        
        # Find corresponding PAE file
        npz_name = cif_name.replace('.cif', '.npz')
        npz_path = pae_dir / npz_name
        
        if not npz_path.exists():
            print(f"Warning: PAE file not found for {cif_name}, skipping...")
            continue
        
        # Get sequence
        sequence = filename_to_seq.get(original_name, "")
        chain_sequence = filename_to_chain_seq.get(original_name, "")
        
        if not sequence and not chain_sequence:
            print(f"Warning: Sequence not found for {original_name}, skipping...")
            continue
        
        # Calculate ipSAE
        print(f"  Processing {cif_name}...")
        ipsae = calculate_ipsae(
            str(npz_path),
            str(cif_path),
            pae_cutoff,
            dist_cutoff
        )
        
        designs.append({
            'id': cif_path.stem,
            'original_name': original_name.replace('.cif', ''),
            'sequence': chain_sequence if chain_sequence else sequence,
            'designed_sequence': sequence,
            'ipsae': ipsae,
            'cif_path': str(cif_path),
        })
    
    if not designs:
        print("Error: No designs were processed successfully")
        sys.exit(1)
    
    # Sort by ipSAE (lower is better)
    designs.sort(key=lambda x: x['ipsae'])
    
    # Add ipSAE rank
    for i, design in enumerate(designs):
        design['ipsae_rank'] = i + 1
    
    # Determine output path
    if output_fasta is None:
        output_fasta = final_dir / "designs_ranked_by_ipsae.fasta"
    
    # Write FASTA
    write_fasta(designs, str(output_fasta))
    print(f"\nFASTA file written to: {output_fasta}")
    
    # Also write a summary CSV
    summary_csv = Path(str(output_fasta).replace('.fasta', '_summary.csv'))
    summary_df = pd.DataFrame([{
        'ipsae_rank': d['ipsae_rank'],
        'id': d['id'],
        'original_name': d['original_name'],
        'ipsae': d['ipsae'],
        'designed_sequence': d['designed_sequence'],
        'chain_sequence': d['sequence'],
    } for d in designs])
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary CSV written to: {summary_csv}")
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"ipSAE Ranking Summary (lower is better)")
    print(f"{'='*60}")
    print(f"{'Rank':<6} {'ID':<40} {'ipSAE':<12}")
    print(f"{'-'*60}")
    for d in designs[:10]:  # Show top 10
        print(f"{d['ipsae_rank']:<6} {d['id']:<40} {d['ipsae']:<12.4f}")
    if len(designs) > 10:
        print(f"... and {len(designs) - 10} more designs")
    print(f"{'='*60}")
    
    return designs


def main():
    parser = argparse.ArgumentParser(
        description="Calculate ipSAE scores for final ranked designs and output ranked FASTA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--final_designs_dir",
        type=str,
        nargs='?',
        default='test/final_ranked_designs',
        help="Path to final_ranked_designs directory"
    )
    parser.add_argument(
        "--pae_cutoff",
        type=float,
        default=15.0,
        help="PAE cutoff for ipSAE calculation (default: 15.0)"
    )
    parser.add_argument(
        "--dist_cutoff",
        type=float,
        default=15.0,
        help="Distance cutoff for ipSAE calculation (default: 15.0)"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default='test/final_ranked_designs/designs_ranked_by_ipsae.fasta',
        help="Output FASTA file path (default: <final_designs_dir>/designs_ranked_by_ipsae.fasta)"
    )
    
    args = parser.parse_args()
    
    # If no directory specified, try to find one
    if args.final_designs_dir is None:
        # Look for common locations
        candidates = [
            "output/final_ranked_designs",
            "test/final_ranked_designs",
            "final_ranked_designs",
        ]
        for candidate in candidates:
            if Path(candidate).exists():
                args.final_designs_dir = candidate
                print(f"Auto-detected final designs directory: {candidate}")
                break
        
        if args.final_designs_dir is None:
            print("Error: Please specify the final_ranked_designs directory")
            parser.print_help()
            sys.exit(1)
    
    process_final_designs(
        args.final_designs_dir,
        args.pae_cutoff,
        args.dist_cutoff,
        args.output
    )


if __name__ == "__main__":
    """
    # Basic usage - auto-detects common directories
    py calc_ipsae.py

    # Specify directory explicitly
    py calc_ipsae.py output/final_ranked_designs

    # Custom cutoffs
    py calc_ipsae.py output/final_ranked_designs --pae_cutoff 12 --dist_cutoff 12

    # Custom output path
    py calc_ipsae.py output/final_ranked_designs -o my_ranked_designs.fasta
    """
    main()
