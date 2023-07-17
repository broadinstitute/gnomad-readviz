import argparse
import logging
import hail as hl
import re
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.annotations import get_adj_expr
from gnomad_qc.v4.resources.basics import gnomad_v4_genotypes
from gnomad_qc.v4.resources.meta import meta

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
        mt.GT.is_haploid() & (mt.meta.sex_karyotype == "XY") & (mt.GT[0] == 1),
    )


def main(args):

    hl.init(log="/select_samples", default_reference="GRCh38", idempotent=True, tmp_dir=args.temp_bucket)
    meta_ht = hl.import_table(args.sample_metadata_tsv, force_bgz=True)
    meta_ht = meta_ht.key_by("s")
    meta_ht = meta_ht.filter(hl.is_defined(meta_ht.cram_path) & hl.is_defined(meta_ht.crai_path), keep=True)
    meta_ht = meta_ht.repartition(1000)
    meta_ht = meta_ht.checkpoint(
        re.sub(".tsv(.b?gz)?", "", args.sample_metadata_tsv) + ".ht", overwrite=True, _read_if_exists=True)

    vds = gnomad_v4_genotypes.vds()

    # see https://github.com/broadinstitute/ukbb_qc/pull/227/files
    if args.test:
        logger.info("Filtering to chrX PAR1 boundary: chrX:2781477-2781900")
        vds = hl.vds.filter_intervals(vds, [hl.parse_locus_interval("chrX:2781477-2781900")])

    v4_qc_meta_ht = meta.ht()

    mt = vds.variant_data
    #mt = vds.variant_data._filter_partitions([41229])

    mt = mt.filter_cols(v4_qc_meta_ht[mt.s].release)

    meta_join = meta_ht[mt.s]
    mt = mt.annotate_cols(
        meta=hl.struct(
            sex_karyotype=meta_join.sex_karyotype,
            cram=meta_join.cram_path,
            crai=meta_join.crai_path,
        )
    )

    logger.info("Adjusting samples' sex ploidy")
    lgt_expr = hl.if_else(
        mt.locus.in_autosome(),
        mt.LGT,
        adjusted_sex_ploidy_expr(mt.locus, mt.LGT, mt.meta.sex_karyotype)
    )
    mt = mt.select_entries(
        "LA", "GQ", LGT=lgt_expr, adj=get_adj_expr(lgt_expr, mt.GQ, mt.DP, mt.LAD)
    )

    logger.info("Filtering to entries meeting GQ, DP and other 'adj' thresholds")
    mt = filter_to_adj(mt)

    logger.info("Performing sparse split multi")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Filter variants with at least one non-ref GT")
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    #logger.info(f"Saving checkpoint")
    #mt = mt.checkpoint(os.path.join(args.temp_bucket, "readviz_select_samples_checkpoint1.vds"),
    #                   overwrite=True, _read_if_exists=True)

    def sample_ordering_expr(mt):
        """For variants that are present in more than 10 samples (or whatever number args.num_samples is set to),
        this sample ordering determines which samples will be used for readviz. Sorting first by GQ ensures that
        samples with the highest genotype quality are shown. For common variants, many samples may have maximum GQ,
        so if we only sort by GQ, the samples for these variants will, in practice, be chosen based on the
        alphabetical order of their sample ids. This can be problematic because the same few samples would end up being
        selected for most common variants, and so a disproportionate number of common variants would need read data from
        these samples. One problem with would be an imbalance in computational load across samples during downstream
        steps of the readviz pipeline - which would cause some shards to take much longer than others.
        To avoid this, and distribute variants more evenly across samples, we add a random number as the second sort key
        so that the choice of samples will be randomized for each common variant.
        """

        return -mt.GQ, hl.rand_unif(0, 1, seed=1)

    logger.info(
        f"For each site, take up to {args.num_samples} samples from each genotype category: het, hom_var, or hemi"
    )
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
        help="Overwrite if output file already exists", action="store_true",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        help="Number of samples to take from each genotype category at each site",
        default=10,
    )
    parser.add_argument(
        "--temp-bucket",
        help="Bucket for intermediate files",
        default="gs://bw2-delete-after-15-days",
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
