version 1.0

task PrintReadVizIntervals {

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

		Int disk_size = ceil(size(ref_fasta, "GiB") + 2*size(input_cram, "GiB") + 10)
	}

    parameter_meta {
        input_cram: {
            localization_optional: true
        }
    }

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		# 1) Convert variants_tsv_bgz to sorted interval list

		zcat ~{variants_tsv_bgz} | awk '{ OFS="\t" } { print("chr"$1, $2 - ~{padding_around_variant} - 1, $2 + ~{padding_around_variant}) }' > variant_windows.bed

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

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_raw_bam = "${output_prefix}.raw.bam"
		File output_raw_bai = "${output_prefix}.raw.bai"
		File variant_windows_interval_list = "variant_windows.interval_list"

		# save small intermediate files for debugging
		File input_cram_header = "header.bam"
		File variant_windows_bed = "variant_windows.bed"
	}

	runtime {
		docker: "gcr.io/broad-mpg-gnomad/gnomad-readviz@sha256:7013fc57e3471617a314b08e2bcefe4711d401f83500c5c57e9a3e79ee8efebd"
		cpu: 1
		preemptible: 1
		memory: "4 GiB"
		disks: "local-disk ${disk_size} HDD"
		zones: "us-east1-b us-east1-c us-east1-d"
	}
}

task RunHaplotypeCallerBamout {

	input {
		File input_bam
		File input_bai
		File variant_windows_interval_list

		File ref_fasta
		File ref_fasta_fai
		File ref_fasta_dict

		String output_prefix

		Int disk_size = ceil(size(ref_fasta, "GiB") + 2*size(input_bam, "GiB") + 10)
	}

	command <<<

		echo --------------; echo "Start - time: $(date)"; set -euxo pipefail; df -kh; echo --------------

		java -XX:+DisableAttachMechanism -XX:MaxHeapSize=1000m -Xmx7500m -jar /gatk/GATK35.jar -T HaplotypeCaller \
			-R ~{ref_fasta} \
			-I "~{input_bam}" \
			-L ~{variant_windows_interval_list} \
			--disable_auto_index_creation_and_locking_when_reading_rods \
			-bamout "~{output_prefix}.bamout.bam" \
			-o "~{output_prefix}.gvcf"

		ls -lh
		echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------
	>>>

	output {
		File output_bamout_bam = "${output_prefix}.bamout.bam"
		File output_bamout_bai = "${output_prefix}.bamout.bai"
		File output_gvcf = "${output_prefix}.gvcf"
		File output_gvcf_idx = "${output_prefix}.gvcf.idx"
	}

	runtime {
		docker: "gcr.io/broad-mpg-gnomad/gnomad-readviz@sha256:7013fc57e3471617a314b08e2bcefe4711d401f83500c5c57e9a3e79ee8efebd"
		cpu: 1
		preemptible: 2
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


workflow PrintReadVizReadsWorkflow {
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

	call PrintReadVizIntervals {
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

	call RunHaplotypeCallerBamout {
		input:
			input_bam = PrintReadVizIntervals.output_raw_bam,
			input_bai = PrintReadVizIntervals.output_raw_bai,
			variant_windows_interval_list = PrintReadVizIntervals.variant_windows_interval_list,
			ref_fasta = ref_fasta,
			ref_fasta_fai = ref_fasta_fai,
			ref_fasta_dict = ref_fasta_dict,
			output_prefix = output_prefix,
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
		File output_raw_bam = PrintReadVizIntervals.output_raw_bam
		File output_raw_bai = PrintReadVizIntervals.output_raw_bai

		File output_bamout_cram = ConvertBamToCram.output_bamout_cram
		File output_bamout_cram_crai = ConvertBamToCram.output_bamout_cram_crai
		File output_gvcf = RunHaplotypeCallerBamout.output_gvcf
		File output_gvcf_idx = RunHaplotypeCallerBamout.output_gvcf_idx

		# save small intermediate files for debugging
		File input_cram_header = PrintReadVizIntervals.input_cram_header
		File variant_windows_bed = PrintReadVizIntervals.variant_windows_bed
		File variant_windows_interval_list = PrintReadVizIntervals.variant_windows_interval_list
	}
}
