import argparse
import logging
import hail as hl
import re
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.filtering import filter_to_adj
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

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


def main(args):

    hl.init(log="/select_samples", default_reference="GRCh38", idempotent=True)
    meta_ht = hl.import_table(args.sample_metadata_tsv, force_bgz=True)
    meta_ht = meta_ht.key_by("s")
    meta_ht = meta_ht.repartition(1000)
    meta_ht = meta_ht.checkpoint(re.sub(".tsv(.b?gz)?", "") + ".ht", overwrite=True, _read_if_exists=True)

    vds = get_gnomad_v4_vds(split=True, release_only=True)

    # see https://github.com/broadinstitute/ukbb_qc/pull/227/files
    if args.test:
        logger.info("Filtering to chrX PAR1 boundary: chrX:2781477-2781900")
        vds = hl.vds.filter_intervals(vds, [hl.parse_locus_interval("chrX:2781477-2781900")])

    mt = hl.vds.to_dense_mt(vds)

    meta_join = meta_ht[mt.s]
    mt = mt.annotate_cols(
        meta=hl.struct(
            sex=meta_join.sex_karyotype,
            cram=meta_join.cram_path,
            crai=meta_join.crai_path,
        )
    )

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
    mt = mt.select_entries("GT", "GQ", "DP", "AD")

    logger.info("Filtering to entries meeting GQ, DP and other 'adj' thresholds")
    mt = filter_to_adj(mt)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

    logger.info(
        f"Taking up to {args.num_samples} samples per site where samples are het, hom_var, or hemi"
    )

    def sample_ordering_expr(mt):
        """It can be problematic for downstream steps when several samples have many times more variants selected
        than in other samples. To avoid this, and distribute variants more evenly across samples,
        add a random number as the secondary sort order. This way, when many samples have an identically high GQ
        (as often happens for common variants), the same few samples don't get selected repeatedly for all common
        variants.
        """

        return -mt.GQ, hl.rand_unif(0, 1, seed=1)

    mt = mt.annotate_rows(
        samples_w_het_var=hl.agg.filter(
            het_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), args.num_samples, ordering=sample_ordering_expr(mt)),
        ),
        samples_w_hom_var=hl.agg.filter(
            hom_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), args.num_samples, ordering=sample_ordering_expr(mt)),
        ),
        samples_w_hemi_var=hl.agg.filter(
            hemi_expr(mt),
            hl.agg.take(het_hom_hemi_take_expr(mt), args.num_samples, ordering=sample_ordering_expr(mt)),
        ),
    )

    ht = mt.rows()
    ht = ht.select(ht.samples_w_het_var, ht.samples_w_hom_var, ht.samples_w_hemi_var)
    ht.write(args.output_ht_path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Test on chrX", action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite if object already exists", action="store_true",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        help="Number of samples to take from each genotype category at each site",
        default=10,
    )
    parser.add_argument(
        "--sample-metadata-tsv",
        help="Path of the gnomAD sample metadata TSV with columns: s, cram_path, crai_path, sex_karyotype",
        default="gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.metadata.tsv.gz",
    )
    parser.add_argument(
        "--output-ht-path",
        help="Path for output hail table",
        default="gs://gnomad-readviz/v4.0/gnomad.exomes.v4.0.readviz_crams.ht",
        #default="gs://gnomad/readviz/genomes_v3/gnomad_v3_readviz_crams.ht",
    )
    args = parser.parse_args()

    main(args)


"""
Output schema looks like:

----------------------------------------
Global fields:
    None
----------------------------------------
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'samples_w_het_var': array<struct {
        S: str,
        GQ: int32
    }>
    'samples_w_hom_var': array<struct {
        S: str,
        GQ: int32
    }>
    'samples_w_hemi_var': array<struct {
        S: str,
        GQ: int32
    }>
----------------------------------------
Key: ['locus', 'alleles']
----------------------------------------

"""
