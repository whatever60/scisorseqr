"""
Make a exon x cell matrix from mapping result obtained in the previous step.
The input file is a tsv with these 9 columns:
`Read`: <run_id.read_idx>
`Gene`: <ensemble_gene_id>
`Celltype`: a categorical variable
`Barcode`: a DNA sequence marking its cell origin
`A unknown column`
`A unknown column`
`Exons`: `;%;<chr>_<exon_start>_<exon_end>_<strand>;%;...`
`type`: either `knwon` or `novel`
`a unknown column`

Usage: python3 ExonCount.py <input_file> -g <grouping_factor> -t <threshold> -p <num_threads>
"""

import os
import argparse

import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "info_file", type=str, help="File with results from the last step"
    )
    parser.add_argument("--grouping_factor", "-g", type=str, help="Grouping factor")
    parser.add_argument(
        "--num_threads",
        "-p",
        type=int,
        help="Number of threads for parallel processing",
    )
    parser.add_argument(
        "--threshold", "-t", type=int, help="Threshold for filtering data"
    )
    return parser.parse_args()


def checkSpanningReads(
    uniqExons_g, readSE_g, inclusionCounts_g, grouping_factor
) -> pd.DataFrame:
    exons = uniqExons_g.copy()
    exons[["chr", "s", "e", "strand"]] = exons["Exons"].str.split("_", expand=True)
    reads = readSE_g
    spanningReads = reads.merge(exons, on="Gene", how="left")
    spanningReads = spanningReads.astype({"s": int, "e": int, "start": int, "end": int})
    spanningReads = spanningReads[
        (spanningReads["s"] >= spanningReads["start"])
        & (spanningReads["e"] <= spanningReads["end"])
    ]
    spanningReads = (
        spanningReads[["Exons", "Gene", grouping_factor]]
        .assign(
            Total=spanningReads.groupby(["Exons", "Gene", grouping_factor]).transform(
                "size"
            )
        )
        .drop_duplicates()
    )
    inclusionReads = inclusionCounts_g
    inc_tot = inclusionReads.merge(
        spanningReads, on=["Exons", "Gene", grouping_factor], how="right"
    ).fillna(0)
    inc_tot["PSI"] = inc_tot.inclusion / inc_tot.Total
    inc_tot["exclusion"] = inc_tot.Total - inc_tot.inclusion
    inc_tot = inc_tot[
        ["Exons", "Gene", grouping_factor, "Total", "PSI", "inclusion", "exclusion"]
    ].astype({"inclusion": int, "exclusion": int})
    return inc_tot


def main():
    args = parse_arguments()
    info_file = args.info_file
    threshold = args.threshold
    grouping_factor = args.grouping_factor
    num_threads = args.num_threads

    allInfo = pd.read_csv(info_file, sep="\t")
    if allInfo.shape[1] == 11:
        allInfo = allInfo.iloc[:, [0, 1, 2, 3, 8]]
    else:
        allInfo = allInfo.iloc[:, [0, 1, 2, 3, 6]]

    allInfo.columns = ["Read", "Gene", "Celltype", "Barcode", "Exons"]
    allInfo[["start", "end"]] = (
        allInfo["Exons"].str.split("_").apply(lambda x: (x[1], x[-2])).tolist()
    )
    allInfo_SE = allInfo.groupby("Gene").filter(lambda x: len(x) >= threshold)
    if allInfo_SE.shape[0] == 0:
        print("No gene with enough reads")
        return

    # for each read, get rid of the first and last row (i.e., the first and last )
    internalExons = allInfo_SE.copy()
    internalExons["Exons"] = (
        internalExons["Exons"].str.split(";%;").apply(lambda x: x[1:])
    )
    internalExons = internalExons[internalExons.Exons.apply(lambda x: len(x) >= 3)]
    internalExons["Exons"] = internalExons.Exons.apply(lambda x: x[1:-1])
    internalExons = internalExons.explode("Exons")

    uniqExons = internalExons[["Exons", "Gene"]].drop_duplicates()

    inclusionCounts = internalExons[["Exons", "Gene", grouping_factor]].copy()
    inclusionCounts["inclusion"] = inclusionCounts.groupby(
        ["Exons", "Gene", grouping_factor]
    ).transform("size")
    inclusionCounts.drop_duplicates(inplace=True)

    readSE = internalExons[
        ["Read", "Gene", grouping_factor, "start", "end"]
    ].drop_duplicates()

    uniqExons_gs = uniqExons.groupby("Gene")
    readSE_gs = readSE.groupby("Gene")
    inclusionCounts_gs = inclusionCounts.groupby("Gene")
    res = pd.concat(
        Parallel(n_jobs=num_threads)(
            delayed(checkSpanningReads)(
                uniqExons_gs.get_group(gene),
                readSE_gs.get_group(gene),
                inclusionCounts_gs.get_group(gene),
                grouping_factor,
            )
            for gene in tqdm(uniqExons["Gene"].unique())
        )
    )

    outDir = "ExonQuantOutput"
    os.makedirs(outDir, exist_ok=True)
    res.to_csv(f"{outDir}/InclusionExclusionCounts.tsv.gz", index=False, sep="\t")


if __name__ == "__main__":
    main()
