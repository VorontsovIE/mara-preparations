source ./pipeline_core.sh

stage_00_1() {
  pushd source_data/genome
      wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      gzip -d hg38.fa.gz
  popd
}

stage_00_2() {
  pushd source_data/motifs
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_annotation.jsonl
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pwm.tar.gz
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pcm.tar.gz
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pfm.tar.gz
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_thresholds.tar.gz
      for FN in $(ls *.tar.gz); do tar -zxf "${FN}"; done
  popd
}

stage_00_3() {
  pushd app
      wget https://raw.githubusercontent.com/autosome-ru/sarus/master/releases/sarus-2.1.0.jar
      ln -s sarus-2.1.0.jar sarus.jar
  popd
}

stage_00() {
  stage_00_1
  stage_00_2
  stage_00_3
}

stage_01() {
  mkdir -p ./stages/stage_01.preprocess/
  mkdir -p ./stages/stage_01/

  ln -s ../../source_data/Credible_TSSclusters_regions.bed stages/stage_01/
  ln -s ../../source_data/Credible_TSSclusters_regions_summit.bed stages/stage_01/
  ln -s ../../source_data/Credible_TSSclusters_TMMnormalized_logCPM.tsv stages/stage_01/

  ln -s ../../source_data/hg38_promoters_v1.bed stages/stage_01.preprocess/

  # chr1    65418   65419   promoter        65418   +   ——→   hg38_v1_chr1_+_65419_65419
  cat stages/stage_01.preprocess/hg38_promoters_v1.bed \
    | awk -F $'\t' -e '{NAME="hg38_v1_"$1"_"$6"_"($2+1)"_"$3; print $1 "\t" $2 "\t" $3 "\t" NAME "\t" $5 "\t" $6}' \
    > stages/stage_01/hg38_promoters_v1.bed
}

stage_02() {
  mkdir -p ./stages/stage_02/
  make_flanks 250u 50d  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 400u 100d ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 300u 100d ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 200u 50d  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters

  make_flanks 300u 50d  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 400u 50d  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 250u 25d  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters

  make_flanks 250u 50d  ./stages/stage_01/hg38_promoters_v1.bed ./stages/stage_02/hg38_promoters_v1
  
  # upstream
  make_flanks 250u   0  ./stages/stage_01/Credible_TSSclusters_regions_summit.bed ./stages/stage_02/TSS_cluster_summit
  make_flanks 250u 10d  ./stages/stage_01/Credible_TSSclusters_regions_summit.bed ./stages/stage_02/TSS_cluster_summit
  make_flanks 250u 20d  ./stages/stage_01/Credible_TSSclusters_regions_summit.bed ./stages/stage_02/TSS_cluster_summit

  # downstream
  make_flanks 10u  50d  ./stages/stage_01/Credible_TSSclusters_regions_summit.bed ./stages/stage_02/TSS_cluster_summit
  make_flanks 0    50d  ./stages/stage_01/Credible_TSSclusters_regions_summit.bed ./stages/stage_02/TSS_cluster_summit
  make_flanks 10d  50d  ./stages/stage_01/Credible_TSSclusters_regions_summit.bed ./stages/stage_02/TSS_cluster_summit
}

stage_03() {
  mkdir -p ./stages/stage_03/
  flanks_fasta 250u 50d  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 400u 100d ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 300u 100d ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 200u 50d  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters

  flanks_fasta 300u 50d  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 400u 50d  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 250u 25d  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters

  flanks_fasta 250u 50d  ./stages/stage_02/hg38_promoters_v1  ./stages/stage_03/hg38_promoters_v1

  # upstream
  flanks_fasta 250u 0  ./stages/stage_02/TSS_cluster_summit  ./stages/stage_03/TSS_cluster_summit
  flanks_fasta 250u 10d  ./stages/stage_02/TSS_cluster_summit  ./stages/stage_03/TSS_cluster_summit
  flanks_fasta 250u 20d  ./stages/stage_02/TSS_cluster_summit  ./stages/stage_03/TSS_cluster_summit

  # downstream
  flanks_fasta 10u 50d  ./stages/stage_02/TSS_cluster_summit  ./stages/stage_03/TSS_cluster_summit
  flanks_fasta 0   50d  ./stages/stage_02/TSS_cluster_summit  ./stages/stage_03/TSS_cluster_summit
  flanks_fasta 10d 50d  ./stages/stage_02/TSS_cluster_summit  ./stages/stage_03/TSS_cluster_summit
}

stage_04() {
  NUM_THREADS=${1:-10} # 10 threads by default

  mkdir -p ./stages/stage_04/
  (
    motif_occupancies_flanks_cmd 250u 50d  ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
    motif_occupancies_flanks_cmd 400u 100d ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
    motif_occupancies_flanks_cmd 300u 100d ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
    motif_occupancies_flanks_cmd 200u 50d  ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04

    motif_occupancies_flanks_cmd 250u 50d  ./stages/stage_03/hg38_promoters_v1  hg38_promoters_v1  ./stages/stage_04
  ) | parallel -j ${NUM_THREADS}

  (
    motif_occupancies_flanks_cmd 300u 50d  ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
    motif_occupancies_flanks_cmd 400u 50d  ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
    motif_occupancies_flanks_cmd 250u 25d  ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
  ) | parallel -j ${NUM_THREADS}

  (
    motif_besthits_flanks_cmd 250u 50d  ./stages/stage_03/TSS_clusters  TSS  ./stages/stage_04
  ) | parallel -j ${NUM_THREADS}

  (
   motif_occupancies_flanks_cmd 250u 50d  ./stages/stage_03/hg38_promoters_v1  hg38_promoters_v1  ./stages/stage_04
   motif_besthits_flanks_cmd 250u 50d  ./stages/stage_03/hg38_promoters_v1  hg38_promoters_v1  ./stages/stage_04
  ) | parallel -j ${NUM_THREADS}


  (
    # upstream
    motif_occupancies_flanks_cmd 250u 0    ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_occupancies_flanks_cmd 250u 10d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_occupancies_flanks_cmd 250u 20d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04

    # downstream
    motif_occupancies_flanks_cmd 10u 50d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_occupancies_flanks_cmd 0   50d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_occupancies_flanks_cmd 10d 50d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
  ) | parallel -j ${NUM_THREADS}

  (
    # upstream
    motif_besthits_flanks_cmd 250u 0    ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_besthits_flanks_cmd 250u 10d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_besthits_flanks_cmd 250u 20d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04

    # downstream
    motif_besthits_flanks_cmd 10u 50d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_besthits_flanks_cmd 0   50d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
    motif_besthits_flanks_cmd 10d 50d  ./stages/stage_03/TSS_cluster_summit  TSS_summit  ./stages/stage_04
  ) | parallel -j ${NUM_THREADS}
}

stage_00 # downloading data

stage_01
stage_02
stage_03
stage_04 175 # adjust for number of threads available

