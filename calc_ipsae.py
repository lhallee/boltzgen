import pandas as pd
import subprocess


def calculate_ipsae(
    pae_file_path,
    structure_file_path,
    pae_cutoff=15.0,
    dist_cutoff=15.0,
):
    """
    Calculate ipSAE and related scores for protein-protein interactions.

    Parameters:
    -----------
    pae_file_path : str
        Path to the PAE file (JSON for AF2/AF3, NPZ for Boltz1)
    structure_file_path : str
        Path to the structure file (PDB for AF2, mmCIF for AF3/Boltz1)
    pae_cutoff : float
        Cutoff value for PAE in score calculations
    dist_cutoff : float
        Cutoff value for distance in score calculations

    Returns:
    --------
    dict
        Dictionary containing all calculated scores
    """

    try:
        command = [
            "python",
            "ipsae.py",
            pae_file_path,
            structure_file_path,
            str(pae_cutoff),
            str(dist_cutoff),
        ]
        subprocess.run(command)
    except:
        command = [
            "py",
            "ipsae.py",
            pae_file_path,
            structure_file_path,
            str(pae_cutoff),
            str(dist_cutoff),
        ]
        subprocess.run(command)
    
    print(f"Reading results from {structure_file_path.replace('.cif',  f'_{int(pae_cutoff)}_{int(dist_cutoff)}.txt')}")

    df = pd.read_csv(structure_file_path.replace('.cif', f'_{int(pae_cutoff)}_{int(dist_cutoff)}.txt'))
    results = {}


    for i, row in df[df.Type=="max"].iterrows():
        chainpair = f"{row['Chn1']}-{row['Chn2']}"

        results[chainpair] = {
            "max": {
                **{col: row[col] for col in df.columns[5:-1]}
            }
        }
        mask = (df['Chn1'] == row['Chn1']) & (df['Chn2'] == row['Chn2']) & (df['Type'] != "max")
        min_vals = df[mask][df.columns[5:-1]].min()
        results[chainpair]["min"] = min_vals.to_dict()

    return results


if __name__ == "__main__":
    from glob import glob
    result_dir = "test/intermediate_designs_inverse_folded"
    cifs = sorted(glob(f"{result_dir}/*.cif"))
    npzs = sorted(glob(f"{result_dir}/*.npz"))
    print(cifs)
    print(npzs)
    for cif, npz in zip(cifs, npzs):
        calculate_ipsae(npz, cif)