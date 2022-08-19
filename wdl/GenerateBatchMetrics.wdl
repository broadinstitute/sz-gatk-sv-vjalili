version 1.0

import "RDTest.wdl" as rdt
import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
import "Utils.wdl" as util
import "GenerateBatchMetricsAlgorithm.wdl" as gbma
import "TestUtils.wdl" as tu

workflow GenerateBatchMetrics {
  input {
    String batch

    File depth_vcf
    File? melt_vcf
    File? scramble_vcf
    File? wham_vcf
    File? manta_vcf

    File pe_file
    File sr_file
    File baf_file
    File rd_file

    File median_file
    File mean_coverage_file
    File ploidy_table

    Int records_per_shard_agg
    Int records_per_shard_agg_sr
    Int records_per_shard_rdtest

    String? additional_gatk_args_agg

    File? svtk_to_gatk_script

    String chr_x
    String chr_y

    Float? java_mem_fraction

    File rmsk
    File segdups
    File ped_file
    File autosome_contigs
    File allosome_contigs
    File reference_dict

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    File? primary_contigs_list  # required if run_module_metrics = true

    String gatk_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_aggregate_callers
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_agg
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
    RuntimeAttr? runtime_attr_get_male_only
    RuntimeAttr? runtime_attr_metrics_file_metrics
    RuntimeAttr? runtime_attr_annotate_overlap
  }

  call util.GetSampleIdsFromVcf {
    input:
      vcf = depth_vcf,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  call util.SubsetPedFile {
    input:
      ped_file = ped_file,
      sample_list = GetSampleIdsFromVcf.out_file,
      subset_name = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_ped
  }

  call GetSampleLists {
    input:
      ped_file = SubsetPedFile.ped_subset_file,
      samples_list = GetSampleIdsFromVcf.out_file,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_sample_list
  }

  Array[File?] vcfs_ = [depth_vcf, manta_vcf, melt_vcf, scramble_vcf, wham_vcf]
  Array[String] algorithms_ = ["depth", "manta", "melt", "scramble", "wham"]
  scatter (i in range(length(algorithms_))) {
    if (defined(vcfs_[i])) {
      if (algorithms_[i] == "depth" || algorithms_[i] == "manta" || algorithms_[i] == "wham") {
        File baf_file_ = baf_file
        File rd_file_ = rd_file
      }
      if (algorithms_[i] == "manta" || algorithms_[i] == "wham") {
        File pe_file_ = pe_file
      }
      if (algorithms_[i] == "manta" || algorithms_[i] == "wham" || algorithms_[i] == "melt" || algorithms_[i] == "scramble") {
        File sr_file_ = sr_file
      }

      call gbma.GenerateBatchMetricsAlgorithm {
        input:
          vcf=select_first([vcfs_[i]]),
          batch=batch,
          algorithm=algorithms_[i],
          pe_file=pe_file_,
          sr_file=sr_file_,
          baf_file=baf_file_,
          rd_file=rd_file_,
          mean_coverage_file=mean_coverage_file,
          median_file=median_file,
          ped_file=ped_file,
          ploidy_table=ploidy_table,
          sample_list = GetSampleIdsFromVcf.out_file,
          male_list = GetSampleLists.male_samples,
          female_list = GetSampleLists.female_samples,
          records_per_shard_pesr=if defined(sr_file_) then records_per_shard_agg_sr else records_per_shard_agg,
          records_per_shard_depth=records_per_shard_rdtest,
          additional_gatk_args=additional_gatk_args_agg,
          chr_x=chr_x,
          chr_y=chr_y,
          rmsk=rmsk,
          segdups=segdups,
          autosome_contigs=autosome_contigs,
          allosome_contigs=allosome_contigs,
          reference_dict=reference_dict,
          java_mem_fraction=java_mem_fraction,
          svtk_to_gatk_script=svtk_to_gatk_script,
          gatk_docker=gatk_docker,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
          linux_docker=linux_docker,
          runtime_attr_scatter_vcf=runtime_attr_scatter_vcf,
          runtime_attr_agg=runtime_attr_agg,
          runtime_attr_format=runtime_attr_format,
          runtime_attr_rdtest=runtime_attr_rdtest,
          runtime_attr_concat_vcfs=runtime_attr_concat_vcfs,
          runtime_attr_annotate_overlap=runtime_attr_annotate_overlap
      }
    }
  }

  call AggregateCallers {
    input:
      batch = batch,
      input_metrics = select_all(GenerateBatchMetricsAlgorithm.out),
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_aggregate_callers
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call tu.MetricsFileMetrics {
      input:
        metrics_file = AggregateCallers.metrics,
        contig_list = select_first([primary_contigs_list]),
        common = false,
        prefix = "GenerateBatchMetrics.~{batch}",
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_metrics_file_metrics
    }
  }

  output {
    File metrics = AggregateCallers.metrics
    File? metrics_file_batchmetrics = MetricsFileMetrics.out
  }
}

task GetSampleLists {
  input {
    File ped_file
    File samples_list
    String sv_base_docker
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
    File male_samples = "male.list"
    File female_samples = "female.list"
  }
  command <<<

    set -eu
    awk -v sex=1 '($5==sex) {print $2}' ~{ped_file} > ped_males.list
    awk -v sex=2 '($5==sex) {print $2}' ~{ped_file} > ped_females.list

    python3 <<CODE
    with open("ped_males.list",'r') as ped_m, open("ped_females.list",'r') as ped_f:
      male_samples = set([x.strip() for x in ped_m.readlines() if x.strip()])
      female_samples = set([x.strip() for x in ped_f.readlines() if x.strip()])
      with open("male.list", 'w') as samples_m, open("female.list",'w') as samples_f, open("~{samples_list}",'r') as samples:
        for line in samples:
          if line.strip():
            if (line.strip() in male_samples):
              samples_m.write(line)
            if (line.strip() in female_samples):
              samples_f.write(line)
    CODE
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task AggregateCallers {
  input {
    String batch
    Array[File] input_metrics
    String sv_pipeline_base_docker
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

  String output_file = "${batch}.metrics"

  output {
    File metrics = "~{output_file}"
  }
  command <<<

    set -eu
    python3 <<CODE
    import pandas as pd
    metrics = ["~{sep='", "' input_metrics}"]
    dfs=[]
    for df in metrics:
      dfs.append(pd.read_table(df))
    df = pd.concat(dfs)
    df.to_csv("~{output_file}", index=False, sep='\t')
    CODE
        
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
