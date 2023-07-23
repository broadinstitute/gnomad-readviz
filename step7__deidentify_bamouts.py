import hail as hl   # used for hadoop file utils
import logging
import os
import pandas as pd
import subprocess
from tqdm import tqdm

from step_pipeline import pipeline, Backend, Localize, Delocalize

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("deidentify bamouts")
logger.setLevel(logging.INFO)

DOCKER_IMAGE = "weisburd/gnomad-readviz@sha256:c8f02e79221d643dcdbfb0d73e7e28bb6ba9ee7e6d108382b1cf0f35d3dab86c"

OUTPUT_BUCKET = "gs://gnomad-readviz/v4.0/deidentified_bamout"


def parse_args(batch_pipeline):
    """Parse command line args."""

    p = batch_pipeline.get_config_arg_parser()
    p.add_argument(
        "--output-dir",
        help="Where to write HaplotypeCaller output bams.",
        default=OUTPUT_BUCKET,
    )
    debugging_group = p.add_mutually_exclusive_group()
    debugging_group.add_argument(
        "-s", "--sample-id",
        help="Only process this sample id",
        action="append",
    )
    debugging_group.add_argument(
        "-n",
        help="For testing, process only the first N samples",
        type=int)
    p.add_argument(
        "--offset",
        help="Skip the first this many sampmles",
        default=0,
        type=int)
    p.add_argument(
        "--random",
        help="Randomly select -n samples",
        action="store_true")
    p.add_argument(
        "--tsv-has-header",
        help="The input tsv files have a header line",
        action="store_true")
    p.add_argument(
        "tsv_and_bamout_paths_table",
        help="A text file containing at least these columns: "
             "sample_id, output_bamout_bam, output_bamout_bai, variants_tsv_bgz")
    args = p.parse_args()

    return p, args


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE,
                  config_file_path="~/.step_pipeline_gnomad")

    parser, args = parse_args(bp)

    df = pd.read_table(args.tsv_and_bamout_paths_table)

    missing_columns = {"sample_id", "output_bamout_bam", "output_bamout_bai", "variants_tsv_bgz"} - set(df.columns)
    if missing_columns:
        parser.error(f"{args.tsv_and_bamout_paths_table} is missing these columns: {missing_columns}")

    if args.n and args.random:
        df = df.sample(n=args.n)

    if args.sample_id:
        df = df[df["sample_id"].isin(args.sample_id)]

    if not args.force:
        paths = bp.precache_file_paths(os.path.join(args.output_dir, f"*.*"))
        logger.info(f"Found {len(paths)} exising .bam, .bai, .db files")

    bp.name = f"step7: deidentify bamouts ({min(args.n or 10**9, (len(df) - args.offset))} samples)"

    #with open("deidentify_bamout.py", "rt") as f:
    #    deidentify_bamouts_script = f.read()

    logging.info(f"Processing {len(df)} samples")
    for i, (_, row) in tqdm(enumerate(df.iterrows()), unit=" samples"):
        if args.offset and i < args.offset:
            continue
        if args.n and i >= args.offset + args.n:
            break

        s1 = bp.new_step(
            f"run HaplotypeCaller: {row.sample_id}",
            arg_suffix=f"step1",
            image=DOCKER_IMAGE,
            step_number=1,
            cpu=0.25,
            memory="lowmem",
            localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
            delocalize_by=Delocalize.COPY,
            output_dir=args.output_dir,
            all_outputs_precached=True)

        local_tsv_path = s1.input(row.variants_tsv_bgz)
        local_bamout_path, _ = s1.inputs(row.output_bamout_bam, row.output_bamout_bai)

#        s1.command(f"""cat <<EOF > deidentify_bamout.py
#{deidentify_bamouts_script}
#EOF""")

        tsv_has_header = "--tsv-has-header" if args.tsv_has_header else ""

        s1.command(f"""echo --------------
time python3 /deidentify_bamout.py {tsv_has_header} "{row.sample_id}" "{local_bamout_path}" "{local_tsv_path}"

samtools sort -o "{row.sample_id}.deidentified.sorted.bam" "{row.sample_id}.deidentified.bam"
mv "{row.sample_id}.deidentified.sorted.bam" "{row.sample_id}.deidentified.bam"
samtools index "{row.sample_id}.deidentified.bam"
ls -lh
""")
        s1.output(f"{row.sample_id}.deidentified.db")
        s1.output(f"{row.sample_id}.deidentified.bam")
        s1.output(f"{row.sample_id}.deidentified.bam.bai")

    bp.run()


if __name__ == "__main__":
    main()


