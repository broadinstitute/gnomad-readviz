import argparse
import logging
import hail as hl
import hail.expr.aggregators as agg
from gnomad_hail.utils.slack import try_slack
from gnomad_hail.utils.gnomad_functions import adjusted_sex_ploidy_expr
from gnomad_qc.v3.resources import get_full_mt, meta_ht_path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Readviz_prep")
logger.setLevel(logging.INFO)


def het_hom_hemi_take_expr(mt):
    return hl.struct(S=mt.s, GQ=mt.GQ)


def het_expr(mt):
    return mt.GT.is_het()


def hom_expr(mt):
    return mt.GT.is_diploid() & mt.GT.is_hom_var()


def hemi_expr(mt):
    return hl.or_missing(
        mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar(),
        mt.GT.is_haploid() & (mt.meta.sex == "male") & (mt.GT[0] == 1),
    )


def dp_threshold_expr(mt, dp_threshold):
    return hl.if_else(
        (mt.meta.sex == "male") & (mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar()),
        dp_threshold / 2,
        dp_threshold,
    )


def main(args):

    hl.init(log="/readviz_prep", default_reference="GRCh38")
    meta_ht = hl.read_table(meta_ht_path)
    mt = get_full_mt(split=False, key_by_locus_and_alleles=True)
    crams = hl.import_table(args.cram_paths).key_by("s")

    dp_threshold = args.dp_threshold
    gq_threshold = args.gq_threshold
    num_samples = args.num_samples

    if args.test:
        logger.info("Filtering to chrX PAR1 boundary: chrX:2781477-2781900")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chrX:2781477-2781900")])

    meta_join = meta_ht[mt.s]
    mt = mt.annotate_cols(
        meta=hl.struct(
            sex=meta_join.sex,
            release=meta_join.release,
            cram=crams[mt.s].final_cram_path,
            crai=crams[mt.s].final_crai_path,
        )
    )
    logger.info("Filtering to releasable samples with a defined cram path")
    mt = mt.filter_cols(mt.meta.release & hl.is_defined(mt.meta.cram))
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Adjusting samples' sex ploidy")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus,
            mt.GT,
            mt.meta.sex,
            xy_karyotype_str="male",
            xx_karyotype_str="female",
        )
    )
    mt = mt.select_entries("GT", "GQ", "DP")

    logger.info("Filtering to entries meeting GQ and DP thresholds")
    mt = mt.filter_entries(
        (mt.GQ >= gq_threshold)
        & (mt.DP >= dp_threshold_expr(mt, dp_threshold))
        & (mt.GT.is_non_ref())
    )
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

    logger.info(
        f"Taking up to {num_samples} samples per site where samples are het, hom_var, or hemi"
    )
    mt = mt.annotate_rows(
        samples_w_het_var=hl.agg.filter(
            het_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), num_samples, ordering=-mt.GQ),
        ),
        samples_w_hom_var=hl.agg.filter(
            hom_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), num_samples, ordering=-mt.GQ),
        ),
        samples_w_hemi_var=hl.agg.filter(
            hemi_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), num_samples, ordering=-mt.GQ),
        ),
    )

    ht = mt.rows()
    ht = ht.select(ht.samples_w_het_var, ht.samples_w_hom_var, ht.samples_w_hemi_var)
    ht.write(args.output_ht_path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", help="Test on chrX", action="store_true")
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite if object already exists", action="store_true"
    )
    parser.add_argument(
        "--dp-threshold",
        type=int,
        help="Lower depth threshold for sample list",
        default=10,
    )
    parser.add_argument(
        "--gq-threshold",
        type=int,
        help="Lower genotype quality threshold for sample list",
        default=20,
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        help="Number of samples to take from each genotype category at each site",
        default=10,
    )
    parser.add_argument(
        "--cram-paths",
        help="Path to file containing sample, cram, and crai information with headers: s, final_cram_path, final_crai_path",
        required=True,
    )
    parser.add_argument(
        "--output-ht-path", help="Path for output hail table", required=True
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
