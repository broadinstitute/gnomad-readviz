version 1.0

# GATK v35 command-line from gnomAD v3 GVCF:
# GATKCommandLine.HaplotypeCaller=<ID=HaplotypeCaller,Version=3.5-0-g36282e4,Date="Mon Nov 05 02:49:57 UTC 2018",Epoch=1541386197920,CommandLineOptions="analysis_type=HaplotypeCaller input_file=[/cromwell_root/broad-gotc-prod-cromwell-execution/CramToGvcfWorkflow/d80cadfe-69a4-4710-b163-94eb7e1f7342/call-CramToBam/79064-CVD.roundtrip.bam] showFullBamList=false read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[OverclippedRead] disable_read_filter=[] intervals=[/cromwell_root/broad-gotc-prod-cromwell-execution/CramToGvcfWorkflow/d80cadfe-69a4-4710-b163-94eb7e1f7342/call-ScatterIntervalList/glob-cb4648beeaff920acb03de7603c06f98/1scattered.interval_list] excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=/cromwell_root/broad-references/hg38/v0/Homo_sapiens_assembly38.fasta nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=500 baq=OFF baqGapOpenPenalty=40.0 refactor_NDN_cigar_string=false fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 static_quantized_quals=null round_down_quantized=false disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false no_cmdline_in_header=false sites_only=false never_trim_vcf_format_field=false bcf=false bam_compression=null simplifyBAM=false disable_bam_indexing=false generate_md5=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false variant_index_type=LINEAR variant_index_parameter=128000 reference_window_stop=0 logging_level=INFO log_to_file=null help=false version=false likelihoodCalculationEngine=PairHMM heterogeneousKmerSizeResolution=COMBO_MIN dbsnp=(RodBinding name= source=UNBOUND) dontTrimActiveRegions=false maxDiscARExtension=25 maxGGAARExtension=300 paddingAroundIndels=150 paddingAroundSNPs=20 comp=[] annotation=[StrandBiasBySample] excludeAnnotation=[ChromosomeCounts, FisherStrand, StrandOddsRatio, QualByDepth] group=[Standard, StandardHCAnnotation] debug=false useFilteredReadsForAnnotations=false emitRefConfidence=GVCF bamOutput=null bamWriterType=CALLED_HAPLOTYPES disableOptimizations=false annotateNDA=false heterozygosity=0.001 indel_heterozygosity=1.25E-4 standard_min_confidence_threshold_for_calling=-0.0 standard_min_confidence_threshold_for_emitting=-0.0 max_alternate_alleles=3 input_prior=[] sample_ploidy=2 genotyping_mode=DISCOVERY alleles=(RodBinding name= source=UNBOUND) contamination_fraction_to_filter=0.0 contamination_fraction_per_sample_file=null p_nonref_model=null exactcallslog=null output_mode=EMIT_VARIANTS_ONLY allSitePLs=true gcpHMM=10 pair_hmm_implementation=VECTOR_LOGLESS_CACHING pair_hmm_sub_implementation=ENABLE_ALL always_load_vector_logless_PairHMM_lib=false phredScaledGlobalReadMismappingRate=45 noFpga=false sample_name=null kmerSize=[10, 25] dontIncreaseKmerSizesForCycles=false allowNonUniqueKmersInRef=false numPruningSamples=1 recoverDanglingHeads=false doNotRecoverDanglingBranches=false minDanglingBranchLength=4 consensus=false maxNumHaplotypesInPopulation=128 errorCorrectKmers=false minPruning=2 debugGraphTransformations=false allowCyclesInKmerGraphToGeneratePaths=false graphOutput=null kmerLengthForReadErrorCorrection=25 minObservationsForKmerToBeSolid=20 GVCFGQBands=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99] indelSizeToEliminateInRefModel=10 min_base_quality_score=10 includeUmappedReads=false useAllelesTrigger=false doNotRunPhysicalPhasing=false keepRG=null justDetermineActiveRegions=false dontGenotype=false dontUseSoftClippedBases=false captureAssemblyFailureBAM=false errorCorrectReads=false pcr_indel_model=CONSERVATIVE maxReadsInRegionPerSample=10000 minReadsPerAlignmentStart=10 mergeVariantsViaLD=false activityProfileOut=null activeRegionOut=null activeRegionIn=null activeRegionExtension=null forceActive=false activeRegionMaxSize=null bandPassSigma=null maxProbPropagationDistance=50 activeProbabilityThreshold=0.002 filter_is_too_short_value=30 do_not_require_softclips_both_ends=false min_mapping_quality_score=20 filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false">

# HpalotypeCaller command from wdl: https://github.com/broadinstitute/dsde-pipelines/blob/ed4973300d0772012405996b121151c9d21fffc8/genomes_in_the_cloud/cram_to_gvcf/CramToGvcfWf.wdl#L112-L124
task RunHaplotypeCallerBamout {

	input {
		File variants_tsv_bgz
		File input_cram
		File input_crai

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String requester_pays_project
		String output_prefix

		Int padding_around_variant

		Int disk_size = ceil(size(ref_fasta, "GiB") + 2 * size(input_cram, "GiB") + 10)
	}

	#parameter_meta {
	#    input_cram: {
	#        localization_optional: true
	#    }
	#}

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		# 1) Convert variants_tsv_bgz to sorted interval list

		zcat "~{variants_tsv_bgz}" | awk '{ OFS="\t" } { print("chr"$1, $2 - ~{padding_around_variant} - 1, $2 + ~{padding_around_variant}) }' > variant_windows.bed

		# bedtools sort -i variant_windows.bed | bedtools merge -d 200 > variant_windows.merged.bed

		# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
		# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
		# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
		# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

		[[ "$(dirname '~{input_cram}')" != "$(dirname '~{input_crai}')" ]] && cp "~{input_crai}"  $(dirname '~{input_cram}')

		java -Xms2g -jar /gatk/gatk.jar PrintReadsHeader \
			--gcs-project-for-requester-pays ~{requester_pays_project} \
			-R ~{ref_fasta} \
			-I "~{input_cram}" \
			-O header.bam

		java -Xms2g -jar /gatk/gatk.jar BedToIntervalList \
			--SORT true \
			--SEQUENCE_DICTIONARY header.bam \
			--INPUT variant_windows.bed \
			--OUTPUT variant_windows.interval_list

		# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

		java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+DisableAttachMechanism -XX:MaxHeapSize=1000m -Xmx7500m \
			-jar /gatk/GATK35.jar \
			-T HaplotypeCaller \
			-R ~{ref_fasta} \
			-I "~{input_cram}" \
			-L variant_windows.interval_list \
			--disable_auto_index_creation_and_locking_when_reading_rods \
			-ERC GVCF \
			--max_alternate_alleles 3 \
			-variant_index_parameter 128000 \
			-variant_index_type LINEAR \
			--read_filter OverclippedRead \
			-bamout "~{output_prefix}.bamout.bam" \
			-o "~{output_prefix}.gvcf"  |& grep -v "^DEBUG"

		bgzip "~{output_prefix}.gvcf"
		tabix "~{output_prefix}.gvcf.gz"

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_bamout_bam = "${output_prefix}.bamout.bam"
		File output_bamout_bai = "${output_prefix}.bamout.bai"
		File output_gvcf_gz = "${output_prefix}.gvcf.gz"
		File output_gvcf_gz_tbi = "${output_prefix}.gvcf.gz.tbi"

		# save small intermediate files for debugging
		File input_cram_header = "header.bam"
		File variant_windows_bed = "variant_windows.bed"
		File variant_windows_interval_list = "variant_windows.interval_list"
	}

	runtime {
		docker: "gcr.io/broad-mpg-gnomad/gnomad-readviz@sha256:7013fc57e3471617a314b08e2bcefe4711d401f83500c5c57e9a3e79ee8efebd"
		cpu: 1
		preemptible: 1
		memory: "8 GiB"
		disks: "local-disk ${disk_size} HDD"
		zones: "us-east1-b us-east1-c us-east1-d"
	}
}



task ConvertBamToCram {

	input {
		File input_bam
		File input_bai

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String output_prefix

		Int disk_size = ceil(2*size(input_bam, "GiB") + 10)
	}

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		samtools view -T ~{ref_fasta} -C "~{input_bam}" > "~{output_prefix}.bamout.cram"
		samtools index "~{output_prefix}.bamout.cram" "~{output_prefix}.bamout.cram.crai"

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_bamout_cram = "${output_prefix}.bamout.cram"
		File output_bamout_cram_crai = "${output_prefix}.bamout.cram.crai"
	}

	runtime {
		docker: "gcr.io/broad-mpg-gnomad/gnomad-readviz@sha256:7013fc57e3471617a314b08e2bcefe4711d401f83500c5c57e9a3e79ee8efebd"
		cpu: 1
		preemptible: 1
		memory: "2 GiB"
		disks: "local-disk ${disk_size} HDD"
		zones: "us-east1-b us-east1-c us-east1-d"
	}
}


workflow GetReadVizDataWorkflow {
	input {
		File variants_tsv_bgz
		File input_cram
		File input_crai

		File ref_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
		File ref_fasta_fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
		File ref_fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

		String requester_pays_project = "broad-mpg-gnomad"
		String output_prefix = sub(basename(input_cram), "\\.cram$", "")

		Int PADDING_AROUND_VARIANT = 200
	}

	call RunHaplotypeCallerBamout {
		input:
			variants_tsv_bgz = variants_tsv_bgz,
			input_cram = input_cram,
			input_crai = input_crai,
			ref_fasta = ref_fasta,
			ref_fasta_fai = ref_fasta_fai,
			ref_fasta_dict = ref_fasta_dict,
			requester_pays_project = requester_pays_project,
			output_prefix = output_prefix,
			padding_around_variant = PADDING_AROUND_VARIANT,
	}

	call ConvertBamToCram {
		input:
			input_bam = RunHaplotypeCallerBamout.output_bamout_bam,
			input_bai = RunHaplotypeCallerBamout.output_bamout_bai,
			ref_fasta = ref_fasta,
			ref_fasta_fai = ref_fasta_fai,
			ref_fasta_dict = ref_fasta_dict,
			output_prefix = output_prefix,
	}

	output {
		File output_bamout_cram = ConvertBamToCram.output_bamout_cram
		File output_bamout_cram_crai = ConvertBamToCram.output_bamout_cram_crai
		File output_gvcf_gz = RunHaplotypeCallerBamout.output_gvcf_gz
		File output_gvcf_gz_tbi = RunHaplotypeCallerBamout.output_gvcf_gz_tbi

		# save small intermediate files for debugging
		File input_cram_header = RunHaplotypeCallerBamout.input_cram_header
		File variant_windows_bed = RunHaplotypeCallerBamout.variant_windows_bed
		File variant_windows_interval_list = RunHaplotypeCallerBamout.variant_windows_interval_list
	}
}
