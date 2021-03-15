cwlVersion: v1.0
class: Workflow
doc: "{'extract_fastq': 'http://218.77.58.141:9619','extract_fastq2': 'http://218.77.58.141:9619','extract_fastq3': 'http://121.46.19.86:9619','extract_fastq4': 'http://121.46.19.86:9619', 'fastx_quality_stats': 'http://218.77.58.141:9619','fastx_quality_stats2': 'http://218.77.58.141:9619','fastx_quality_stats3': 'http://121.46.19.86:9619','fastx_quality_stats4': 'http://121.46.19.86:9619','bowtie_aligner': 'http://218.77.58.141:9619','bowtie_aligner2': 'http://218.77.58.141:9619', 'bowtie_aligner3': 'http://121.46.19.86:9619', 'bowtie_aligner4': 'http://121.46.19.86:9619', 'samtools_sort_index': 'http://218.77.58.141:9619','samtools_sort_index2': 'http://218.77.58.141:9619','samtools_sort_index3': 'http://121.46.19.86:9619','samtools_sort_index4': 'http://121.46.19.86:9619', 'samtools_rmdup': 'http://218.77.58.141:9619', 'samtools_rmdup2': 'http://218.77.58.141:9619', 'samtools_rmdup3': 'http://121.46.19.86:9619', 'samtools_rmdup4': 'http://121.46.19.86:9619', 'samtools_sort_index_after_rmdup': 'http://218.77.58.141:9619', 'samtools_sort_index_after_rmdup2': 'http://218.77.58.141:9619', 'samtools_sort_index_after_rmdup3': 'http://121.46.19.86:9619', 'samtools_sort_index_after_rmdup4': 'http://121.46.19.86:9619', 'macs2_callpeak': 'http://218.77.58.141:9619', 'macs2_callpeak2': 'http://218.77.58.141:9619','macs2_callpeak3': 'http://121.46.19.86:9619','macs2_callpeak4': 'http://121.46.19.86:9619', 'get_stat': 'http://218.77.58.141:9619', 'get_stat2': 'http://218.77.58.141:9619', 'get_stat3': 'http://121.46.19.86:9619', 'get_stat4': 'http://121.46.19.86:9619', 'island_intersect': 'http://218.77.58.141:9619', 'island_intersect2': 'http://218.77.58.141:9619', 'island_intersect3': 'http://121.46.19.86:9619', 'island_intersect4': 'http://121.46.19.86:9619', 'average_tag_density': 'http://218.77.58.141:9619', 'average_tag_density2': 'http://218.77.58.141:9619', 'average_tag_density3': 'http://121.46.19.86:9619', 'average_tag_density4': 'http://121.46.19.86:9619' }"
requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


inputs:

  indices_folder:
    type: Directory
    label: "BOWTIE indices folder"
    doc: "Path to BOWTIE generated indices folder"

  annotation_file:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated input annotation file"

  genome_size:
    type: string
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

  chrom_length:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosome length file"

  control_file:
    type: File?
    default: null
    label: "Control BAM file"
    format: "http://edamontology.org/format_2572"
    doc: "Control BAM file file for MACS2 peak calling"

  broad_peak:
    type: boolean
    label: "Callpeak broad"
    doc: "Set to call broad peak for MACS2"

  fastq_file1:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"

  fastq_file2:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"

  fastq_file3:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"

  fastq_file4:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"
  
  exp_fragment_size:
    type: int?
    default: 150
    label: "Expected fragment size"
    doc: "Expected fragment size for MACS2"

  force_fragment_size:
    type: boolean?
    default: false
    label: "Force fragment size"
    doc: "Force MACS2 to use exp_fragment_size"

  clip_3p_end:
    type: int?
    default: 0
    label: "Clip from 3p end"
    doc: "Number of bases to clip from the 3p end"

  clip_5p_end:
    type: int?
    default: 0
    label: "Clip from 5p end"
    doc: "Number of bases to clip from the 5p end"

  remove_duplicates:
    type: boolean?
    default: false
    label: "Remove duplicates"
    doc: "Calls samtools rmdup to remove duplicates from sortesd BAM file"

  threads:
    type: int?
    default: 48
    doc: "Number of threads for those steps that support multithreading"
    label: "Number of threads"

outputs:

 # bigwig:
 #   type: File
 #   format: "http://edamontology.org/format_3006"
 #   label: "BigWig file"
 #   doc: "Generated BigWig file"
 #   outputSource: bam_to_bigwig/bigwig_file

  fastx_statistics:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats/statistics_file

  fastx_statistics2:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats2/statistics_file

  fastx_statistics3:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats3/statistics_file

  fastx_statistics4:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats4/statistics_file

  bowtie_log:
    type: File
    label: "BOWTIE alignment log"
    format: "http://edamontology.org/format_2330"
    doc: "BOWTIE generated alignment log"
    outputSource: bowtie_aligner/log_file

  bowtie_log2:
    type: File
    label: "BOWTIE alignment log"
    format: "http://edamontology.org/format_2330"
    doc: "BOWTIE generated alignment log"
    outputSource: bowtie_aligner2/log_file

  bowtie_log3:
    type: File
    label: "BOWTIE alignment log"
    format: "http://edamontology.org/format_2330"
    doc: "BOWTIE generated alignment log"
    outputSource: bowtie_aligner3/log_file

  bowtie_log4:
    type: File
    label: "BOWTIE alignment log"
    format: "http://edamontology.org/format_2330"
    doc: "BOWTIE generated alignment log"
    outputSource: bowtie_aligner4/log_file


  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_rmdup/bam_bai_pair

  bambai_pair2:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_rmdup2/bam_bai_pair

  bambai_pair3:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_rmdup3/bam_bai_pair

  bambai_pair4:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_rmdup4/bam_bai_pair


  get_stat_log:
    type: File?
    label: "Bowtie & Samtools Rmdup combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup log"
    outputSource: get_stat/output_file

  get_stat_log2:
    type: File?
    label: "Bowtie & Samtools Rmdup combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup log"
    outputSource: get_stat2/output_file

  get_stat_log3:
    type: File?
    label: "Bowtie & Samtools Rmdup combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup log"
    outputSource: get_stat3/output_file

  get_stat_log4:
    type: File?
    label: "Bowtie & Samtools Rmdup combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup log"
    outputSource: get_stat4/output_file

  macs2_called_peaks:
    type: File?
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "XLS file to include information about called peaks"
    outputSource: macs2_callpeak/peak_xls_file

  macs2_called_peaks2:
    type: File?
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "XLS file to include information about called peaks"
    outputSource: macs2_callpeak2/peak_xls_file

  macs2_called_peaks3:
    type: File?
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "XLS file to include information about called peaks"
    outputSource: macs2_callpeak3/peak_xls_file

  macs2_called_peaks4:
    type: File?
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "XLS file to include information about called peaks"
    outputSource: macs2_callpeak4/peak_xls_file

  macs2_narrow_peaks:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak/narrow_peak_file

  macs2_narrow_peaks2:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak2/narrow_peak_file

  macs2_narrow_peaks3:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak3/narrow_peak_file

  macs2_narrow_peaks4:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak4/narrow_peak_file

  macs2_broad_peaks:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak/broad_peak_file

  macs2_broad_peaks2:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak2/broad_peak_file

  macs2_broad_peaks3:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak3/broad_peak_file

  macs2_broad_peaks4:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak4/broad_peak_file

  macs2_peak_summits:
    type: File?
    label: "Peak summits"
    format: "http://edamontology.org/format_3003"
    doc: "Contains the peak summits locations for every peaks"
    outputSource: macs2_callpeak/peak_summits_file

  macs2_peak_summits2:
    type: File?
    label: "Peak summits"
    format: "http://edamontology.org/format_3003"
    doc: "Contains the peak summits locations for every peaks"
    outputSource: macs2_callpeak2/peak_summits_file

  macs2_peak_summits3:
    type: File?
    label: "Peak summits"
    format: "http://edamontology.org/format_3003"
    doc: "Contains the peak summits locations for every peaks"
    outputSource: macs2_callpeak3/peak_summits_file

  macs2_peak_summits4:
    type: File?
    label: "Peak summits"
    format: "http://edamontology.org/format_3003"
    doc: "Contains the peak summits locations for every peaks"
    outputSource: macs2_callpeak4/peak_summits_file

  macs2_moder_r:
    type: File?
    label: "MACS2 generated R script"
    format: "http://edamontology.org/format_2330"
    doc: "R script to produce a PDF image about the model based on your data"
    outputSource: macs2_callpeak/moder_r_file
 
  macs2_moder_r2:
    type: File?
    label: "MACS2 generated R script"
    format: "http://edamontology.org/format_2330"
    doc: "R script to produce a PDF image about the model based on your data"
    outputSource: macs2_callpeak2/moder_r_file

  macs2_moder_r3:
    type: File?
    label: "MACS2 generated R script"
    format: "http://edamontology.org/format_2330"
    doc: "R script to produce a PDF image about the model based on your data"
    outputSource: macs2_callpeak3/moder_r_file

  macs2_moder_r4:
    type: File?
    label: "MACS2 generated R script"
    format: "http://edamontology.org/format_2330"
    doc: "R script to produce a PDF image about the model based on your data"
    outputSource: macs2_callpeak4/moder_r_file

  macs2_gapped_peak:
    type: File?
    label: "Gapped peak"
    format: "http://edamontology.org/format_3586"
    doc: "Contains both the broad region and narrow peaks"
    outputSource: macs2_callpeak/gapped_peak_file

  macs2_gapped_peak2:
    type: File?
    label: "Gapped peak"
    format: "http://edamontology.org/format_3586"
    doc: "Contains both the broad region and narrow peaks"
    outputSource: macs2_callpeak2/gapped_peak_file

  macs2_gapped_peak3:
    type: File?
    label: "Gapped peak"
    format: "http://edamontology.org/format_3586"
    doc: "Contains both the broad region and narrow peaks"
    outputSource: macs2_callpeak3/gapped_peak_file

  macs2_gapped_peak4:
    type: File?
    label: "Gapped peak"
    format: "http://edamontology.org/format_3586"
    doc: "Contains both the broad region and narrow peaks"
    outputSource: macs2_callpeak4/gapped_peak_file

  macs2_log:
    type: File?
    label: "MACS2 log"
    format: "http://edamontology.org/format_2330"
    doc: "MACS2 output log"
    outputSource: macs2_callpeak/macs_log

  macs2_log2:
    type: File?
    label: "MACS2 log"
    format: "http://edamontology.org/format_2330"
    doc: "MACS2 output log"
    outputSource: macs2_callpeak2/macs_log

  macs2_log3:
    type: File?
    label: "MACS2 log"
    format: "http://edamontology.org/format_2330"
    doc: "MACS2 output log"
    outputSource: macs2_callpeak3/macs_log

  macs2_log4:
    type: File?
    label: "MACS2 log"
    format: "http://edamontology.org/format_2330"
    doc: "MACS2 output log"
    outputSource: macs2_callpeak4/macs_log

  iaintersect_log:
    type: File
    label: "Island intersect log"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated log"
    outputSource: island_intersect/log_file

  iaintersect_log2:
    type: File
    label: "Island intersect log"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated log"
    outputSource: island_intersect2/log_file

  iaintersect_log3:
    type: File
    label: "Island intersect log"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated log"
    outputSource: island_intersect3/log_file

  iaintersect_log4:
    type: File
    label: "Island intersect log"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated log"
    outputSource: island_intersect4/log_file

  iaintersect_result:
    type: File
    label: "Island intersect results"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated results"
    outputSource: island_intersect/result_file

  iaintersect_result2:
    type: File
    label: "Island intersect results"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated results"
    outputSource: island_intersect2/result_file

  iaintersect_result3:
    type: File
    label: "Island intersect results"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated results"
    outputSource: island_intersect3/result_file

  iaintersect_result4:
    type: File
    label: "Island intersect results"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated results"
    outputSource: island_intersect4/result_file

  atdp_log:
    type: File
    label: "ATDP log"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated log"
    outputSource: average_tag_density/log_file

  atdp_log2:
    type: File
    label: "ATDP log"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated log"
    outputSource: average_tag_density2/log_file

  atdp_log3:
    type: File
    label: "ATDP log"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated log"
    outputSource: average_tag_density3/log_file

  atdp_log4:
    type: File
    label: "ATDP log"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated log"
    outputSource: average_tag_density4/log_file

  atdp_result:
    type: File
    label: "ATDP results"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: average_tag_density/result_file

  atdp_result2:
    type: File
    label: "ATDP results"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: average_tag_density2/result_file

  atdp_result3:
    type: File
    label: "ATDP results"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: average_tag_density3/result_file

  atdp_result4:
    type: File
    label: "ATDP results"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: average_tag_density4/result_file

steps:

  extract_fastq:
    run: ./tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file1
    out: [fastq_file]
 
  extract_fastq2:
    run: ./tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file2
    out: [fastq_file]

  extract_fastq3:
    run: ./tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file3
    out: [fastq_file]

  extract_fastq4:
    run: ./tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file4
    out: [fastq_file]

  fastx_quality_stats:
    run: ./tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq/fastq_file
    out: [statistics_file]
  
  fastx_quality_stats2:
    run: ./tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq2/fastq_file
    out: [statistics_file]
 
  fastx_quality_stats3:
    run: ./tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq3/fastq_file
    out: [statistics_file]
   
  fastx_quality_stats4:
    run: ./tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq4/fastq_file
    out: [statistics_file]
  
  bowtie_aligner:
    run: ./tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq/fastq_file
      indices_folder: indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
      q:
        default: true
    out: [sam_file, log_file]

  bowtie_aligner2:
    run: ./tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq2/fastq_file
      indices_folder: indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
      q:
        default: true
    out: [sam_file, log_file]

  bowtie_aligner3:
    run: ./tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq3/fastq_file
      indices_folder: indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
      q:
        default: true
    out: [sam_file, log_file]

  bowtie_aligner4:
    run: ./tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq4/fastq_file
      indices_folder: indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
      q:
        default: true
    out: [sam_file, log_file]

  samtools_sort_index:
    run: ./tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner/sam_file
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index2:
    run: ./tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner2/sam_file
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index3:
    run: ./tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner3/sam_file
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index4:
    run: ./tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner4/sam_file
      threads: threads
    out: [bam_bai_pair]

  samtools_rmdup:
    run: ./tools/samtools-rmdup.cwl
    in:
      trigger: remove_duplicates
      bam_file: samtools_sort_index/bam_bai_pair
      single_end:
        default: true
    out: [rmdup_output, rmdup_log]

  samtools_rmdup2:
    run: ./tools/samtools-rmdup.cwl
    in:
      trigger: remove_duplicates
      bam_file: samtools_sort_index2/bam_bai_pair
      single_end:
        default: true
    out: [rmdup_output, rmdup_log]

  samtools_rmdup3:
    run: ./tools/samtools-rmdup.cwl
    in:
      trigger: remove_duplicates
      bam_file: samtools_sort_index3/bam_bai_pair
      single_end:
        default: true
    out: [rmdup_output, rmdup_log]

  samtools_rmdup4:
    run: ./tools/samtools-rmdup.cwl
    in:
      trigger: remove_duplicates
      bam_file: samtools_sort_index4/bam_bai_pair
      single_end:
        default: true
    out: [rmdup_output, rmdup_log]

  samtools_sort_index_after_rmdup:
    run: ./tools/samtools-sort-index.cwl
    in:
      trigger: remove_duplicates
      sort_input: samtools_rmdup/rmdup_output
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index_after_rmdup2:
    run: ./tools/samtools-sort-index.cwl
    in:
      trigger: remove_duplicates
      sort_input: samtools_rmdup2/rmdup_output
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index_after_rmdup3:
    run: ./tools/samtools-sort-index.cwl
    in:
      trigger: remove_duplicates
      sort_input: samtools_rmdup3/rmdup_output
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index_after_rmdup4:
    run: ./tools/samtools-sort-index.cwl
    in:
      trigger: remove_duplicates
      sort_input: samtools_rmdup4/rmdup_output
      threads: threads
    out: [bam_bai_pair]

  macs2_callpeak:
    run: ./tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: samtools_sort_index_after_rmdup/bam_bai_pair
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format_mode:
        default: BAM
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - peak_summits_file
      - broad_peak_file
      - moder_r_file
      - gapped_peak_file
      - treat_pileup_bdg_file
      - control_lambda_bdg_file
      - macs_log
      - macs2_stat_file
      - macs2_fragments_calculated

  macs2_callpeak2:
    run: ./tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: samtools_sort_index_after_rmdup2/bam_bai_pair
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format_mode:
        default: BAM
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - peak_summits_file
      - broad_peak_file
      - moder_r_file
      - gapped_peak_file
      - treat_pileup_bdg_file
      - control_lambda_bdg_file
      - macs_log
      - macs2_stat_file
      - macs2_fragments_calculated

  macs2_callpeak3:
    run: ./tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: samtools_sort_index_after_rmdup3/bam_bai_pair
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format_mode:
        default: BAM
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - peak_summits_file
      - broad_peak_file
      - moder_r_file
      - gapped_peak_file
      - treat_pileup_bdg_file
      - control_lambda_bdg_file
      - macs_log
      - macs2_stat_file
      - macs2_fragments_calculated

  macs2_callpeak4:
    run: ./tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: samtools_sort_index_after_rmdup4/bam_bai_pair
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format_mode:
        default: BAM
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - peak_summits_file
      - broad_peak_file
      - moder_r_file
      - gapped_peak_file
      - treat_pileup_bdg_file
      - control_lambda_bdg_file
      - macs_log
      - macs2_stat_file
      - macs2_fragments_calculated

  get_stat:
      run: ./tools/python-get-stat-chipseq.cwl
      in:
        bowtie_log: bowtie_aligner/log_file
        rmdup_log: samtools_rmdup/rmdup_log
      out:
        - output_file
        - mapped_reads

  get_stat2:
      run: ./tools/python-get-stat-chipseq.cwl
      in:
        bowtie_log: bowtie_aligner2/log_file
        rmdup_log: samtools_rmdup2/rmdup_log
      out:
        - output_file
        - mapped_reads

  get_stat3:
      run: ./tools/python-get-stat-chipseq.cwl
      in:
        bowtie_log: bowtie_aligner3/log_file
        rmdup_log: samtools_rmdup3/rmdup_log
      out:
        - output_file
        - mapped_reads

  get_stat4:
      run: ./tools/python-get-stat-chipseq.cwl
      in:
        bowtie_log: bowtie_aligner4/log_file
        rmdup_log: samtools_rmdup4/rmdup_log
      out:
        - output_file
        - mapped_reads


  island_intersect:
      run: ./tools/iaintersect.cwl
      in:
        input_filename: macs2_callpeak/peak_xls_file
        annotation_filename: annotation_file
        promoter_bp:
          default: 1000
      out: [result_file, log_file]

  island_intersect2:
      run: ./tools/iaintersect.cwl
      in:
        input_filename: macs2_callpeak2/peak_xls_file
        annotation_filename: annotation_file
        promoter_bp:
          default: 1000
      out: [result_file, log_file]

  island_intersect3:
      run: ./tools/iaintersect.cwl
      in:
        input_filename: macs2_callpeak3/peak_xls_file
        annotation_filename: annotation_file
        promoter_bp:
          default: 1000
      out: [result_file, log_file]

  island_intersect4:
      run: ./tools/iaintersect.cwl
      in:
        input_filename: macs2_callpeak4/peak_xls_file
        annotation_filename: annotation_file
        promoter_bp:
          default: 1000
      out: [result_file, log_file]

  average_tag_density:
      run: ./tools/atdp.cwl
      in:
        input_file: samtools_sort_index_after_rmdup/bam_bai_pair
        annotation_filename: annotation_file
        fragmentsize_bp: macs2_callpeak/macs2_fragments_calculated
        avd_window_bp:
          default: 5000
        avd_smooth_bp:
          default: 50
        ignore_chr:
          default: chrM
        double_chr:
          default: "chrX chrY"
        avd_heat_window_bp:
          default: 200
        mapped_reads: get_stat/mapped_reads
      out: [result_file, log_file]

  average_tag_density2:
      run: ./tools/atdp.cwl
      in:
        input_file: samtools_sort_index_after_rmdup2/bam_bai_pair
        annotation_filename: annotation_file
        fragmentsize_bp: macs2_callpeak2/macs2_fragments_calculated
        avd_window_bp:
          default: 5000
        avd_smooth_bp:
          default: 50
        ignore_chr:
          default: chrM
        double_chr:
          default: "chrX chrY"
        avd_heat_window_bp:
          default: 200
        mapped_reads: get_stat2/mapped_reads
      out: [result_file, log_file]

  average_tag_density3:
      run: ./tools/atdp.cwl
      in:
        input_file: samtools_sort_index_after_rmdup3/bam_bai_pair
        annotation_filename: annotation_file
        fragmentsize_bp: macs2_callpeak3/macs2_fragments_calculated
        avd_window_bp:
          default: 5000
        avd_smooth_bp:
          default: 50
        ignore_chr:
          default: chrM
        double_chr:
          default: "chrX chrY"
        avd_heat_window_bp:
          default: 200
        mapped_reads: get_stat3/mapped_reads
      out: [result_file, log_file]

  average_tag_density4:
      run: ./tools/atdp.cwl
      in:
        input_file: samtools_sort_index_after_rmdup4/bam_bai_pair
        annotation_filename: annotation_file
        fragmentsize_bp: macs2_callpeak4/macs2_fragments_calculated
        avd_window_bp:
          default: 5000
        avd_smooth_bp:
          default: 50
        ignore_chr:
          default: chrM
        double_chr:
          default: "chrX chrY"
        avd_heat_window_bp:
          default: 200
        mapped_reads: get_stat4/mapped_reads
      out: [result_file, log_file]
