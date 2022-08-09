version 1.0

import "RDTest.wdl" as rdt
import "BAFTest.wdl" as baft
import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
import "Utils.wdl" as util
import "GenerateBatchMetricsMetrics.wdl" as metrics
import "GeneratePESRMetrics.wdl" as pesr_metrics

workflow GenerateBatchMetricsAlgorithm {
	input {
		String batch
		String algorithm

		File vcf

		File? pe_file
		File? sr_file
		File? baf_file
		File? coveragefile

		File medianfile
		File sample_list
		File female_list
		File male_list
		File ped_file

		File mean_coverage_file
		File ploidy_table

		Int records_per_shard_pesr
		Int records_per_shard_depth

		String? additional_gatk_args_pesr_metrics

		File reference_dict
		String chr_x
		String chr_y

		Float? java_mem_fraction_pesr_metrics

		String gatk_docker
		String sv_base_mini_docker
		String sv_pipeline_docker


		String sv_pipeline_docker
		String sv_pipeline_rdtest_docker
		String sv_base_mini_docker
		String linux_docker

		RuntimeAttr? runtime_attr_scatter_pesr_metrics
		RuntimeAttr? runtime_attr_agg_pesr
		RuntimeAttr? runtime_override_generate_metrics
		RuntimeAttr? runtime_attr_annotate_overlap
		RuntimeAttr? runtime_attr_aggregate_tests

		RuntimeAttr? runtime_attr_rdtest
		RuntimeAttr? runtime_attr_split_rd_vcf
		RuntimeAttr? runtime_attr_merge_allo
		RuntimeAttr? runtime_attr_merge_stats
		RuntimeAttr? runtime_attr_get_male_only

		Int common_cnv_size_cutoff

		File rmsk
		File segdups
		File ped_file
		File autosome_contigs
		File allosome_contigs
		File ref_dict

		# Module metrics parameters
		# Run module metrics workflow at the end - on by default
		Boolean? run_module_metrics

	}

	call pesr_metrics.GeneratePESRBAFMetrics {
		input:
			vcf=vcf,
			prefix="~{batch}.generate_pesr_metrics.manta",
			mean_coverage_file=mean_coverage_file,
			ploidy_table=ploidy_table,
			pe_file=pe_file,
			sr_file=sr_file,
			baf_file=baf_file,
			records_per_shard=records_per_shard_pesr,
			additional_gatk_args=additional_gatk_args_pesr_metrics,
			chr_x=chr_x,
			chr_y=chr_y,
			java_mem_fraction=java_mem_fraction_pesr_metrics,
			gatk_docker=gatk_docker,
			sv_base_mini_docker=sv_base_mini_docker,
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_override_concat=runtime_override_generate_metrics
	}

	if (defined(coveragefile)) {
		call GetMaleOnlyVariantIDs {
			input:
				vcf = vcf,
				female_samples = female_list,
				male_samples = male_list,
				contig = select_first([chr_x, "chrX"]),
				sv_pipeline_docker = sv_pipeline_docker,
				runtime_attr_override = runtime_attr_get_male_only
		}
		call rdt.RDTest {
			input:
				vcf = vcf,
				algorithm = algorithm,
				coveragefile = coveragefile,
				medianfile = medianfile,
				ped_file = ped_file,
				autosome_contigs = autosome_contigs,
				split_size = records_per_shard_depth,
				flags = "",
				allosome_contigs = allosome_contigs,
				ref_dict = ref_dict,
				batch = batch,
				samples = sample_list,
				male_samples = male_list,
				female_samples = female_list,
				male_only_variant_ids = GetMaleOnlyVariantIDs.male_only_variant_ids,
				sv_pipeline_docker = sv_pipeline_docker,
				sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
				linux_docker = linux_docker,
				runtime_attr_rdtest = runtime_attr_rdtest,
				runtime_attr_split_rd_vcf = runtime_attr_split_rd_vcf,
				runtime_attr_merge_allo = runtime_attr_merge_allo,
				runtime_attr_merge_stats = runtime_attr_merge_stats
		}
	}

	call SVRegionOverlap {
		input:
			vcf = GeneratePESRBAFMetrics.out,
			vcf_index = GeneratePESRBAFMetrics.out_index,
			reference_dict = ref_dict,
			output_prefix = "~{batch}.~{algorithm}.annotate_overlap",
			region_files = [segdups, rmsk],
			region_names = ["SEGDUP", "RMSK"],
			gatk_docker = gatk_docker,
			runtime_attr_override = runtime_attr_annotate_overlap
	}

	call AggregateTests {
		input:
			vcf = GeneratePESRBAFMetrics.out,
			prefix = "~{batch}.~{algorithm}",
			rdtest = RDTest.rdtest,
			sv_pipeline_docker = sv_pipeline_docker,
			runtime_attr_override = runtime_attr_aggregate_tests
	}

	output {
		File out = AggregateTests.out
		File out_index = AggregateTests.out_index
	}

}

task GetMaleOnlyVariantIDs {
	input {
		File vcf
		File female_samples
		File male_samples
		String contig
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 3.75,
															 disk_gb: 10,
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
		File male_only_variant_ids = "male_only_variant_ids.txt"
	}
	command <<<
		set -euxo pipefail
		bcftools view -t ~{contig} -S ~{male_samples} ~{vcf} | bcftools view --min-ac 1 | bcftools query -f '%ID\n' > variant_ids_in_males.txt
		bcftools view -t ~{contig} -S ~{female_samples} ~{vcf} | bcftools view --min-ac 1 | bcftools query -f '%ID\n' > variant_ids_in_females.txt
		awk 'NR==FNR{a[$0];next} !($0 in a)' variant_ids_in_females.txt variant_ids_in_males.txt > male_only_variant_ids.txt
	>>>
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: sv_pipeline_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task SVRegionOverlap {
	input {
		File vcf
		File vcf_index
		File reference_dict
		String output_prefix
		Array[File] region_files
		Array[String] region_names

		String? region_set_rule
		String? region_merging_rule
		Int? region_padding

		Boolean? suppress_overlap_fraction
		Boolean? suppress_endpoint_counts

		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 3.75,
															 disk_gb: 10,
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
	Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

	output {
		File out = "~{output_prefix}.vcf.gz"
		File out_index = "~{output_prefix}.vcf.gz.tbi"
	}
	command <<<

		set -euo pipefail
		gatk --java-options -Xmx~{java_mem_mb}M SVRegionOverlap \
			-V ~{vcf} \
			-O ~{output_prefix}.vcf.gz \
			--sequence-dictionary ~{reference_dict} \
			--region-file ~{sep=" --region-file " region_files} \
			--region-name ~{sep=" --region-name " region_names} \
			~{"--region-set-rule " + region_set_rule} \
			~{"--region-merging-rule " + region_merging_rule} \
			~{"--region-padding " + region_padding} \
			~{"--suppress-overlap-fraction " + suppress_overlap_fraction} \
			~{"--suppress-endpoint-counts " + suppress_endpoint_counts}

	>>>
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: gatk_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}

task AggregateTests {
	input {
		File vcf
		String prefix
		File? rdtest
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 7.5,
															 disk_gb: 10,
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
		File out = "~{prefix}.vcf.gz"
		File out_index = "~{prefix}.vcf.gz.tbi"
	}
	command <<<
		/opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py \
			-v ~{vcf} \
			~{"-r " + rdtest} \
			~{prefix}.vcf.gz
		tabix ~{prefix}.vcf.gz
	>>>
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: sv_pipeline_docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}