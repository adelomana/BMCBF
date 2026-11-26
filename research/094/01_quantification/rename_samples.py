#!/usr/bin/env python3
"""Replace the pipeline sample names in a TPM table with the original sample
names from the sequencing facility, sanitised for computational use.

The original names, e.g. "(AEG13) DMSO - 1", become e.g. "AEG13_DMSO_1". The
AEG13 / AEG14 prefix is the batch, so the new names make the batch structure
explicit: AZM appears only in AEG13 and EP317 only in AEG14, while DMSO spans
both batches.

Note that the numbering restarts within each batch: the old DMSO_4 becomes
AEG14_DMSO_1, not AEG14_DMSO_4.
"""

import os
import pandas

PROFILES_DIR = "/Users/adrian/research/bmcbf/094/results/profiles"
INPUT_NAME = "DESeq2_TPM_values.tsv"
OUTPUT_NAME = "DESeq2_TPM_values.newnames.tsv"

# old pipeline name -> sanitised original name
NAME_MAP = {
    "AZM_1": "AEG13_AZM_1",
    "AZM_2": "AEG13_AZM_2",
    "AZM_3": "AEG13_AZM_3",
    "AZM_CPAD_1": "AEG13_AZM_CPAD_1",
    "AZM_CPAD_2": "AEG13_AZM_CPAD_2",
    "AZM_CPAD_3": "AEG13_AZM_CPAD_3",
    "DMSO_1": "AEG13_DMSO_1",
    "DMSO_2": "AEG13_DMSO_2",
    "DMSO_3": "AEG13_DMSO_3",
    "DMSO_CPAD_1": "AEG13_DMSO_CPAD_1",
    "DMSO_CPAD_2": "AEG13_DMSO_CPAD_2",
    "DMSO_CPAD_3": "AEG13_DMSO_CPAD_3",
    "DMSO_4": "AEG14_DMSO_1",
    "DMSO_5": "AEG14_DMSO_2",
    "DMSO_6": "AEG14_DMSO_3",
    "DMSO_CPAD_4": "AEG14_DMSO_CPAD_1",
    "DMSO_CPAD_5": "AEG14_DMSO_CPAD_2",
    "DMSO_CPAD_6": "AEG14_DMSO_CPAD_3",
    "EP317_1": "AEG14_EP317_1",
    "EP317_2": "AEG14_EP317_2",
    "EP317_3": "AEG14_EP317_3",
    "EP317_CPAD_1": "AEG14_EP317_CPAD_1",
    "EP317_CPAD_2": "AEG14_EP317_CPAD_2",
    "EP317_CPAD_3": "AEG14_EP317_CPAD_3",
}


def main():
    input_path = os.path.join(PROFILES_DIR, INPUT_NAME)
    output_path = os.path.join(PROFILES_DIR, OUTPUT_NAME)

    table = pandas.read_csv(input_path, sep="\t", index_col=0)
    print("read {}: {} rows, {} columns".format(input_path, table.shape[0], table.shape[1]))

    # Rename by name, never by position, and fail loudly on any mismatch.
    found = list(table.columns)
    unexpected = [name for name in found if name not in NAME_MAP]
    unused = [name for name in NAME_MAP if name not in found]
    if unexpected or unused:
        raise SystemExit(
            "column names do not match the mapping\n"
            "  columns in the file with no mapping: {}\n"
            "  mapping entries not in the file:     {}".format(unexpected, unused)
        )

    renamed = table.rename(columns=NAME_MAP)
    if len(set(renamed.columns)) != len(renamed.columns):
        raise SystemExit("the mapping produced duplicate column names")

    renamed.to_csv(output_path, sep="\t")
    print("wrote {}: {} rows, {} columns".format(output_path, renamed.shape[0], renamed.shape[1]))
    for old, new in zip(found, renamed.columns):
        print("  {:16} -> {}".format(old, new))


if __name__ == "__main__":
    main()
