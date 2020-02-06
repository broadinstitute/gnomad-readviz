import argparse
import logging
import hail as hl
import hail.expr.aggregators as agg
from gnomad_hail.utils.slack import try_slack
from gnomad_hail.utils.gnomad_functions import get_adj_expr, adjusted_sex_ploidy_expr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Readviz_prep")
logger.setLevel(logging.INFO)


def get_expr_for_het_hom_hemi_take(mt, var_type: str):
    return hl.struct(
        S=mt.s,
        sex=mt.meta.sex,
        cram=mt.meta.cram,
        crai=mt.meta.crai,
        GT=mt.GT,
        GQ=mt.GQ,
        DP=mt.DP,
        AD=mt.AD,
        het_or_hom_or_hemi=var_type,
    )


def get_expr_for_het(mt, gq_threshold, dp_threshold):
    return (
        mt.GT.is_het()
        & (mt.GQ >= gq_threshold)
        & (mt.DP >= dp_threshold)
        & hl.is_defined(mt.meta.cram)
    )


def get_expr_for_hom(mt, gq_threshold, dp_threshold):
    return (
        mt.GT.is_diploid()
        & mt.GT.is_hom_var()
        & (mt.GQ >= gq_threshold)
        & (mt.DP >= dp_threshold)
        & hl.is_defined(mt.meta.cram)
    )


def get_expr_for_hemi(mt, gq_threshold, dp_threshold):
    return (
        mt.GT.is_haploid()
        & (mt.meta.sex == "male")
        & (mt.GT[0] == 1)
        & (mt.GQ >= gq_threshold)
        & (mt.DP >= dp_threshold)
        & hl.is_defined(mt.meta.cram)
    )


def main(args):

    hl.init(log="/readviz_prep")
    meta_ht = hl.read_table(args.meta_ht_path)
    mt = hl.read_matrix_table(args.variant_mt_path).key_rows_by("locus", "alleles")
    crams = hl.import_table(args.cram_paths).key_by("s")

    dp_threshold = args.dp_threshold
    gq_threshold = args.gq_threshold
    num_samples = args.num_samples

    meta_join = meta_ht[mt.s]
    mt = mt.annotate_cols(
        meta=hl.struct(
            project_id=meta_join.project_id,
            sex=meta_join.sex,
            release=meta_join.release,
            sample_filters=meta_join.sample_filters,
            cram=crams[mt.s].final_cram_path,
            crai=crams[mt.s].final_crai_path,
        )
    )
    mt = mt.filter_cols(mt.meta.release)
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Annotating genotypes with adj...")

    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.meta.sex),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    logger.info("Selecting only entry fields needed for densify and output ...")
    mt = mt.select_entries("GT", "GQ", "DP", "AD", "END", "adj")

    mt = hl.experimental.densify(mt)

    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

    logger.info(
        f"Taking up to {num_samples} samples per variant where samples are het, hom_var, or hemi"
    )
    mt = mt.annotate_rows(
        samples_w_het_var=hl.agg.filter(
            get_expr_for_het(mt, gq_threshold, dp_threshold),
            hl.agg.take(
                get_expr_for_het_hom_hemi_take(mt, "het"), num_samples, ordering=-mt.GQ,
            ),
        ),
        samples_w_hom_var=hl.agg.filter(
            get_expr_for_hom(mt, gq_threshold, dp_threshold),
            hl.agg.take(
                get_expr_for_het_hom_hemi_take(mt, "hom"), num_samples, ordering=-mt.GQ
            ),
        ),
        samples_w_hemi_var=hl.agg.filter(
            get_expr_for_hemi(mt, gq_threshold, dp_threshold),
            hl.agg.take(
                get_expr_for_het_hom_hemi_take(mt, "hemi"), num_samples, ordering=-mt.GQ
            ),
        ),
    )

    ht = mt.rows()
    ht = ht.select(
        samples=ht.samples_w_het_var.extend(
            ht.samples_w_hom_var.extend(ht.samples_w_hemi_var)
        )
    )
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
        "--variant-mt-path",
        help="Path to sparse variant matrix table containing variants",
        required=True,
    )
    parser.add_argument(
        "--meta-ht-path", help="Path to sample metadata hail table", required=True
    )
    parser.add_argument(
        "--output-ht-path", help="Path for output hail table", required=True
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
