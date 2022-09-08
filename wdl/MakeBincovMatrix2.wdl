version 1.0

import "Structs.wdl"

task MakeBincovMatrix {
  input {
    Array[File] rd_files
    File? bincov_matrix
    String batch
    Int? binsize
    String sv_base_docker
    RuntimeAttr? runtime_attr_override
  }

  # Only compressed files are stored (if localization_optional, then only output file is stored),
  # so this is a reasonably conservative estimate for disk:
  Int disk_gb = 10 + ceil(2.0 * (size(rd_files, "GiB") + size(bincov_matrix, "GiB")))
  # Some memory is used up by the named pipes. Not a lot, but allocate in case the batch is huge:
  Float mem_gb = 3.0 + 0.003 * length(rd_files)
  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: mem_gb,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged_bincov = "~{batch}.rd.txt.gz"
    File merged_bincov_idx = "~{batch}.rd.txt.gz.tbi"
  }

  command <<<
    set -Eeuo pipefail

    # determine bin size, and drop all bins not exactly equal to this size
    if ~{defined(binsize)}; then
      # use the provided bin size
      binsize=~{binsize}
    else
      # use the most common bin size from the bins
      binsize=$(
        sed -n '2,1000p' ~{rd_files[0]} | awk '{ print $3-$2 }' \
        | sort | uniq -c | sort -nrk1,1 \
        | sed -n '1p' | awk '{ print $2 }'
      )
    fi

    # Use loops rather than WDL sep feature in case there are enough samples to exceed bash line limits
    # Use named pipes to stream unzipped column files in memory
    mkdir -p fifos
    mkfifo fifos/++++++ # peculiar filename to glob before the numbered files
    zcat ~{rd_files[0]} | awk -v bsz=$binsize '{FS=OFS="\t"}{if($3-$2==bsz)print $0}' | cut -f1-3 > fifos/++++++ &

    FILE_NUM=0
    while read -r RD_FILE; do
      FIFO=$(printf "fifos/%06d" $FILE_NUM) # filename is just a bare 6-digit number
      mkfifo $FIFO
      zcat $RD_FILE | awk -v bsz=$binsize '{FS=OFS="\t"}{if($3-$2==bsz)print $0}' | cut -f4-  > $FIFO &
      ((++FILE_NUM))
    done < ~{write_lines(rd_files)}

    if ~{defined(bincov_matrix)}' then
      mkfifo fifos/aaaaaa # peculiar filename to glob after the numbered files
      zcat ~{bincov_matrix} | awk -v bsz=$binsize '{FS=OFS="\t"}{if($3-$2==bsz)print $0}' | cut -f4- > fifos/aaaaaa &
    fi

    # paste unzipped files and compress
    cd fifos
    paste * | bgzip > "../~{batch}.rd.txt.gz"
    cd ..
    tabix -p bed "~{batch}.rd.txt.gz"
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
