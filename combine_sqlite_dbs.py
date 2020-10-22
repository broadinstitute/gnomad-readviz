import argparse
import collections
import gzip
import peewee
from pprint import pprint
import logging
import os
import pysam
import random
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

SQLITE_BATCH_SIZE = 10000


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("combined_bamout_id", help="name for this group of bamouts to use as output filename prefix")
    p.add_argument("sqlite_db_paths", help="text file containing single-sample deidentified sqlite db paths (one path per line)")

    args = p.parse_args()

    sqlite_db_paths = []
    with open(args.sqlite_db_paths) as f:
        for line in f:
            path = line.strip()
            if not os.path.isfile(path):
                p.error(f"path {path} in {args.sqlite_db_paths} doesn't exist")
            sqlite_db_paths.append(path)

    return args.combined_bamout_id, sqlite_db_paths


def run(db_filename, sql):
    command = f"sqlite3 {db_filename} '{sql}'"
    print(command)
    os.system(command)


def main():
    combined_bamout_id, sqlite_db_paths = parse_args()

    for chrom in [str(c) for c in range(1, 23)] + ["X", "Y", "M"]:
        output_db_filename = f"{combined_bamout_id}.chr{chrom}.db"
        if os.path.isfile(output_db_filename):
            os.remove(output_db_filename)

        run(output_db_filename,
            'CREATE TABLE "variants" ('
                '"id" INTEGER NOT NULL PRIMARY KEY, '
                '"chrom" VARCHAR(2) NOT NULL, '
                '"pos" INTEGER NOT NULL, '
                '"ref" TEXT NOT NULL, '
                '"alt" TEXT NOT NULL, '
                '"zygosity" INTEGER NOT NULL, '
                '"qual" INTEGER NOT NULL, '
                '"combined_bamout_id" TEXT, '
                '"read_group_id" INTEGER NOT NULL);')

        column_names_string = "chrom, pos, ref, alt, zygosity, qual, combined_bamout_id, read_group_id"
        for input_db_path in sqlite_db_paths:

            run(output_db_filename,
                f'ATTACH "{input_db_path}" as toMerge; '
                f'BEGIN; '
                f'INSERT INTO variants ({column_names_string}) SELECT {column_names_string} FROM toMerge.variants WHERE chrom="{chrom}"; '
                f'COMMIT; '
                f'DETACH toMerge; '
            )

        run(output_db_filename,
            f'UPDATE variants SET combined_bamout_id="{combined_bamout_id}";')

        #run(output_db_filename,
        #    f'SELECT COUNT(*) from variants;')
        #run(output_db_filename,
        #    'CREATE INDEX variant_index ON "variants" ("chrom", "pos", "ref", "alt", "zygosity", "qual");')


if __name__ == "__main__":
    main()
