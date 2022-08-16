version 1.0

import "RDTest.wdl" as rdt
import "BAFTest.wdl" as baft
import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
import "TasksMakeCohortVcf.wdl" as taskscohort
import "Utils.wdl" as util

workflow GenerateBatchMetricsAlgorithm {
	input {
		String batch
		String algorithm

		File vcf
		File? pe_file
		File? sr_file
		File? baf_file
		File? rd_file  # Runs RdTest iff provided

		File mean_coverage_file
		File median_file
		File ped_file
		File ploidy_table

		File sample_list
		File female_list
		File male_list

		Int records_per_shard_pesr
		Int records_per_shard_depth

		String? additional_gatk_args

		File? svtk_to_gatk_script

		File reference_dict
		String chr_x
		String chr_y

		File rmsk
		File segdups
		File ped_file
		File autosome_contigs
		File allosome_contigs

		Float? java_mem_fraction

		String gatk_docker
		String sv_base_mini_docker
		String sv_pipeline_docker
		String sv_pipeline_docker
		String sv_pipeline_rdtest_docker
		String sv_base_mini_docker
		String linux_docker

		RuntimeAttr? runtime_attr_scatter_vcf
		RuntimeAttr? runtime_attr_format
		RuntimeAttr? runtime_attr_agg
		RuntimeAttr? runtime_attr_annotate_overlap
		RuntimeAttr? runtime_attr_aggregate_tests
		RuntimeAttr? runtime_attr_concat_vcfs

		RuntimeAttr? runtime_attr_rdtest
		RuntimeAttr? runtime_attr_split_rd_vcf
		RuntimeAttr? runtime_attr_merge_allo
		RuntimeAttr? runtime_attr_merge_stats
		RuntimeAttr? runtime_attr_get_male_only

		# Module metrics parameters
		# Run module metrics workflow at the end - on by default
		Boolean? run_module_metrics

	}

	String prefix = "~{batch}.generate_batch_metrics.~{algorithm}"

	call taskscohort.ScatterVcf {
		input:
			vcf=vcf,
			records_per_shard = records_per_shard_pesr,
			prefix = "~{prefix}.scatter_vcf",
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_attr_scatter_vcf
	}

	scatter ( i in range(length(ScatterVcf.shards)) ) {
		call FormatVcfForGatk {
			input:
				vcf=ScatterVcf.shards[i],
				ploidy_table=ploidy_table,
				output_prefix="~{prefix}.format.shard_~{i}",
				script=svtk_to_gatk_script,
				sv_pipeline_docker=sv_pipeline_docker,
				runtime_attr_override=runtime_attr_format
		}
		call SVRegionOverlap {
			input:
				vcf = FormatVcfForGatk.out,
				vcf_index = FormatVcfForGatk.out_index,
				reference_dict = reference_dict,
				output_prefix = "~{prefix}.region_overlap.shard_~{i}",
				region_files = [segdups, rmsk],
				region_file_indexes = [segdups + ".tbi", rmsk + ".tbi"],
				region_names = ["SEGDUP", "RMSK"],
				java_mem_fraction=java_mem_fraction,
				gatk_docker = gatk_docker,
				runtime_attr_override = runtime_attr_annotate_overlap
		}
		call Aggregate {
			input:
				vcf = SVRegionOverlap.out,
				vcf_index = SVRegionOverlap.out_index,
				output_prefix = "~{prefix}.aggregate.shard_~{i}",
				mean_coverage_file = mean_coverage_file,
				ploidy_table=ploidy_table,
				pe_file = pe_file,
				sr_file = sr_file,
				baf_file = baf_file,
				chr_x = chr_x,
				chr_y = chr_y,
				additional_args=additional_gatk_args,
				java_mem_fraction = java_mem_fraction,
				gatk_docker = gatk_docker,
				runtime_attr_override = runtime_attr_agg
		}
	}

	call taskscohort.ConcatVcfs {
		input:
			vcfs=Aggregate.out,
			vcfs_idx=Aggregate.out_index,
			naive=true,
			outfile_prefix="~{prefix}.concat",
			sv_base_mini_docker=sv_base_mini_docker,
			runtime_attr_override=runtime_attr_concat_vcfs
	}

	if (defined(rd_file)) {
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
				coveragefile = select_first([rd_file]),
				medianfile = median_file,
				ped_file = ped_file,
				autosome_contigs = autosome_contigs,
				split_size = records_per_shard_depth,
				flags = "",
				allosome_contigs = allosome_contigs,
				ref_dict = reference_dict,
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

	call AggregateTests {
		input:
			vcf = ConcatVcfs.concat_vcf,
			vcf_index = ConcatVcfs.concat_vcf_idx,
			prefix = "~{prefix}",
			rdtest = RDTest.rdtest,
			sv_pipeline_docker = sv_pipeline_docker,
			runtime_attr_override = runtime_attr_aggregate_tests
	}

	output {
		File out = AggregateTests.out
	}

}

task Aggregate {
	input {
		File vcf
		File vcf_index
		String output_prefix

		File mean_coverage_file
		File ploidy_table
		File? pe_file
		File? sr_file
		File? baf_file

		String chr_x
		String chr_y

		String? additional_args

		Float? java_mem_fraction
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	parameter_meta {
		pe_file: {
							 localization_optional: true
						 }
		sr_file: {
							 localization_optional: true
						 }
		baf_file: {
							 localization_optional: true
						 }
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 3.75,
															 disk_gb: ceil(10 + size(vcf, "GB") * 2.5),
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
		File out = "~{output_prefix}.vcf.gz"
		File out_index = "~{output_prefix}.vcf.gz.tbi"
	}
	command <<<
		set -euo pipefail

		function getJavaMem() {
			# get JVM memory in MiB by getting total memory from /proc/meminfo
			# and multiplying by java_mem_fraction
			cat /proc/meminfo \
				| awk -v MEM_FIELD="$1" '{
					f[substr($1, 1, length($1)-1)] = $2
				} END {
					printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
				}'
		}
		JVM_MAX_MEM=$(getJavaMem MemTotal)
		echo "JVM memory: $JVM_MAX_MEM"

		gatk --java-options "-Xmx${JVM_MAX_MEM}" AggregatePairedEndAndSplitReadEvidence \
			-V ~{vcf} \
			-O ~{output_prefix}.vcf.gz \
			--sample-coverage ~{mean_coverage_file} \
			--ploidy-table ~{ploidy_table} \
			--x-chromosome-name ~{chr_x} \
			--y-chromosome-name ~{chr_y} \
			~{"--discordant-pairs-file " + pe_file} \
			~{"--split-reads-file " + sr_file} \
			~{"--baf-file " + baf_file} \
			~{additional_args}
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
															 disk_gb: ceil(10 + size(vcf, "GB")),
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

task FormatVcfForGatk {
	input {
		File vcf
		File ploidy_table
		File? script
		String? remove_infos
		String? remove_formats
		String output_prefix
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 3.75,
															 disk_gb: ceil(10 + size(vcf, "GB") * 2.0),
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
		File out = "~{output_prefix}.vcf.gz"
		File out_index = "~{output_prefix}.vcf.gz.tbi"
	}
	command <<<
		set -euo pipefail
		python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
			--vcf ~{vcf} \
			--out ~{output_prefix}.vcf.gz \
			--ploidy-table ~{ploidy_table} \
			~{"--remove-infos " + remove_infos} \
			~{"--remove-formats " + remove_formats}
		tabix ~{output_prefix}.vcf.gz
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
		Array[File] region_file_indexes
		Array[String] region_names

		String? region_set_rule
		String? region_merging_rule
		Int? region_padding

		Boolean? suppress_overlap_fraction
		Boolean? suppress_endpoint_counts

		Float? java_mem_fraction

		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 3.75,
															 disk_gb: ceil(10 + size(vcf, "GB") * 2.0),
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
		File out = "~{output_prefix}.vcf.gz"
		File out_index = "~{output_prefix}.vcf.gz.tbi"
	}
	command <<<

		set -euo pipefail

		function getJavaMem() {
			# get JVM memory in MiB by getting total memory from /proc/meminfo
			# and multiplying by java_mem_fraction
			cat /proc/meminfo \
				| awk -v MEM_FIELD="$1" '{
					f[substr($1, 1, length($1)-1)] = $2
				} END {
					printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
				}'
			}
		JVM_MAX_MEM=$(getJavaMem MemTotal)
		echo "JVM memory: $JVM_MAX_MEM"

		gatk --java-options "-Xmx${JVM_MAX_MEM}" SVRegionOverlap \
			-V ~{vcf} \
			-O ~{output_prefix}.vcf.gz \
			--sequence-dictionary ~{reference_dict} \
			--region-file ~{sep=" --region-file " region_files} \
			--region-name ~{sep=" --region-name " region_names} \
			~{"--region-set-rule " + region_set_rule} \
			~{"--region-merging-rule " + region_merging_rule} \
			~{"--region-padding " + region_padding} \
			--suppress-overlap-fraction ~{default="false" suppress_overlap_fraction} \
			--suppress-endpoint-counts ~{default="false" suppress_endpoint_counts}

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
		File vcf_index
		String prefix
		File? rdtest
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
															 cpu_cores: 1,
															 mem_gb: 7.5,
															 disk_gb: ceil(50 + size(vcf, "GB")),
															 boot_disk_gb: 10,
															 preemptible_tries: 3,
															 max_retries: 1
														 }
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

	output {
		File out = "~{prefix}.metrics.tsv"
	}
	command <<<
		/opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py \
			-v ~{vcf} \
			~{"-r " + rdtest} \
			~{prefix}.metrics.tsv
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