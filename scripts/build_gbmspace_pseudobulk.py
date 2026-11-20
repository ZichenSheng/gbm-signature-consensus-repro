"""Build donor- and site-level pseudobulks from the GBM-Space snRNA H5AD."""

from argparse import ArgumentParser
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def decode_strings(values):
    return np.array([item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values])


def read_categorical(observations, name):
    field = observations[name]
    categories = decode_strings(field["categories"][:])
    codes = field["codes"][:]
    return np.where(codes < 0, "NA", categories[np.clip(codes, 0, None)])


def build_pseudobulk(input_path, output_dir):
    with h5py.File(input_path, "r") as source:
        observations = source["obs"]
        labels = read_categorical(observations, "annotation_granular")
        donors = read_categorical(observations, "donor_id")
        sites = read_categorical(observations, "site_id")
        phases = read_categorical(observations, "phase")
        genes = decode_strings(source["var"]["SYMBOL"][:])

        donor_groups, donor_indices = np.unique(pd.Series(donors + "||" + labels), return_inverse=True)
        site_groups, site_indices = np.unique(pd.Series(sites + "||" + labels), return_inverse=True)
        donor_sums = np.zeros((len(donor_groups), len(genes)), dtype=np.float64)
        site_sums = np.zeros((len(site_groups), len(genes)), dtype=np.float32)

        counts = source["X"]
        row_boundaries = counts["indptr"][:]
        matrix_data = counts["data"]
        column_indices = counts["indices"]
        chunk_size = 25_000
        for start in range(0, len(labels), chunk_size):
            stop = min(start + chunk_size, len(labels))
            left, right = row_boundaries[start], row_boundaries[stop]
            rows = np.arange(stop - start)
            block = csr_matrix(
                (matrix_data[left:right], column_indices[left:right], row_boundaries[start:stop + 1] - left),
                shape=(len(rows), len(genes)),
            )
            donor_map = csr_matrix(
                (np.ones(len(rows)), (donor_indices[start:stop], rows)),
                shape=(len(donor_groups), len(rows)),
            )
            site_map = csr_matrix(
                (np.ones(len(rows), dtype=np.float32), (site_indices[start:stop], rows)),
                shape=(len(site_groups), len(rows)),
            )
            donor_sums += (donor_map @ block).toarray()
            site_sums += (site_map @ block).toarray().astype(np.float32)

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "gbmspace_pseudobulk.npz",
        donor_sums=donor_sums.astype(np.float32),
        site_sums=site_sums,
        donor_groups=donor_groups.astype(str),
        site_groups=site_groups.astype(str),
        donor_n=pd.Series(donor_indices).value_counts().sort_index().to_numpy(),
        site_n=pd.Series(site_indices).value_counts().sort_index().to_numpy(),
        genes=genes,
        donor_cycling_frac=pd.Series(np.isin(phases, ["S", "G2M"])).groupby(donor_indices).mean().to_numpy(),
    )
    pd.crosstab(sites, labels).to_csv(output_dir / "gbmspace_site_label_counts.tsv", sep="\t")
    pd.crosstab(donors, labels).to_csv(output_dir / "gbmspace_donor_label_counts.tsv", sep="\t")
    pd.crosstab([sites, labels], phases).to_csv(output_dir / "gbmspace_site_label_phase_counts.tsv", sep="\t")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="GBM-Space raw-count snRNA H5AD from the study portal")
    parser.add_argument("output_dir", type=Path, help="Empty local output directory")
    args = parser.parse_args()
    if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
        parser.error("output_dir must be absent or an empty directory")
    build_pseudobulk(args.input, args.output_dir)


if __name__ == "__main__":
    main()
