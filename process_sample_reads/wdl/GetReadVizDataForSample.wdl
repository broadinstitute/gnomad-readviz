version 1.0

# For each sample
task GetReadVizDataForSample {

	input {
		File variants_tsv_bgz
        String input_cram
		String input_crai

		File ref_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
		File ref_fasta_fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
		File ref_fasta_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

		String requester_pays_project = "broad-mpg-gnomad"
		String output_prefix = sub(basename(input_cram), "\\.cram$", "")

        Int PADDING_AROUND_VARIANT = 200

		Int disk_size = ceil(size(ref_fasta, "GB") + 10)
	}

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

        # 1) Convert variants_tsv_bgz to .bed

        zcat ~{variants_tsv_bgz} | awk '{ OFS="\t" } { print($1, $2 - ~{PADDING_AROUND_VARIANT} - 1, $2 + ~{PADDING_AROUND_VARIANT}) }' > variant_windows.bed

		# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
		# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
		# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
		# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

		java -Xms2g -jar /gatk/gatk.jar PrintReadsHeader \
			--gcs-project-for-requester-pays ~{requester_pays_project} \
			-R ~{ref_fasta} \
			-I "~{input_cram}" \
			--read-index "~{input_crai}" \
			-O header.bam

		java -Xms2g -jar /gatk/gatk.jar BedToIntervalList \
			--SORT true \
			--SEQUENCE_DICTIONARY header.bam \
			--INPUT variant_windows.bed \
			--OUTPUT variant_windows.interval_list

		# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

		java -Xms2g -jar /gatk/gatk.jar PrintReads \
		    --gcs-project-for-requester-pays ~{requester_pays_project} \
			-R ~{ref_fasta} \
			-I "~{input_cram}" \
			--read-index "~{input_crai}" \
			-L variant_windows.interval_list \
			-O "~{output_prefix}.raw.bam"

        # 3) Run HaplotypeCaller

        java -XX:+DisableAttachMechanism -XX:MaxHeapSize=1000m -Xmx7500m -jar /gatk/GATK35.jar -T HaplotypeCaller \
            -R ~{ref_fasta} \
            -I "~{output_prefix}.raw.bam" \
            -L variant_windows.interval_list \
            --disable_bam_indexing \
            --disable_auto_index_creation_and_locking_when_reading_rods \
            -bamout "~{output_prefix}.bamout.bam" \
            -o "~{output_prefix}.bamout.gvcf"

        # 4) Convert bam to cram

        samtools view -T ~{ref_fasta} -C "~{output_prefix}.bamout.bam" > "~{output_prefix}.bamout.cram"
        samtools index "~{output_prefix}.bamout.cram" "~{output_prefix}.bamout.cram.crai"

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
    >>>

	output {
		File output_gvcf = "${output_prefix}.bamout.gvcf"
		File output_gvcf_idx = "${output_prefix}.bamout.gvcf.idx"
		File output_cram = "${output_prefix}.bamout.cram"
		File output_cram_crai = "${output_prefix}.bamout.cram.crai"

        # save small intermediate files for debugging
        #File header_bam = "header.bam"
        #File variant_windows_bed = "variant_windows.bed"
        #File variant_windows_interval_list = "variant_windows.interval_list"
	}

	runtime {
		docker: "weisburd/gnomad-readviz@sha256:37b7a91920513883e993d589faf378d7c057c8f1f777cfb907332c1443302480"
		cpu: 1
		preemptible: 2
		memory: "4 GB"
		disks: "local-disk ${disk_size} HDD"
	}
}


workflow PrintReadVizReadsWorkflow {
	input {
		File variants_tsv_bgz
        String input_cram
		String input_crai
	}

	call GetReadVizDataForSample {
		input:
			variants_tsv_bgz = variants_tsv_bgz,
			input_cram = input_cram,
			input_crai = input_crai,
	}

	output {
		File output_cram = GetReadVizDataForSample.output_cram
		File output_cram_crai = GetReadVizDataForSample.output_cram_crai
	}
}
