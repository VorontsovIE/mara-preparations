set -euo pipefail

source ./pipeline_core.sh

stage_00_1() {
  pushd source_data/genome
      wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      gzip -d hg38.fa.gz
  popd
}

stage_00_2() {
  rm -rf source_data/motifs/
  mkdir -p source_data/motifs/
  pushd source_data/motifs
      wget https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13CORE-CLUSTERED/H13CORE-CLUSTERED_annotation.jsonl
      wget https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13CORE-CLUSTERED/H13CORE-CLUSTERED_pwm.tar.gz
      wget https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13CORE-CLUSTERED/H13CORE-CLUSTERED_pcm.tar.gz
      wget https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13CORE-CLUSTERED/H13CORE-CLUSTERED_pfm.tar.gz
      wget https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13CORE-CLUSTERED/H13CORE-CLUSTERED_thresholds.tar.gz
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

  ln -s ../../source_data/data-v3/hg38_fair+new_CAGE_peaks_phase1and2.bed stages/stage_01/
}

stage_02() {
  mkdir -p ./stages/stage_02/
  # 7-th column used for peak summit
  make_flanks 250u 10d  ./stages/stage_01/hg38_fair+new_CAGE_peaks_phase1and2.bed  ./stages/stage_02/hg38_fair+new_CAGE_peaks_phase1and2  7
}

stage_03() {
  mkdir -p ./stages/stage_03/
  flanks_fasta 250u 10d ./stages/stage_02/hg38_fair+new_CAGE_peaks_phase1and2  ./stages/stage_03/hg38_fair+new_CAGE_peaks_phase1and2.fa
}

stage_04() {
  NUM_THREADS=${1:-10} # 10 threads by default

  mkdir -p ./stages/stage_04/
  (
    motif_occupancies_flanks_cmd 250u 10d  ./stages/stage_03/hg38_fair+new_CAGE_peaks_phase1and2.fa  hg38_fair+new_CAGE_peaks_phase1and2.bed
    # motif_besthits_flanks_cmd 250u 10d  ./stages/stage_03/hg38_fair+new_CAGE_peaks_phase1and2.fa  hg38_fair+new_CAGE_peaks_phase1and2.bed
  ) | parallel -j ${NUM_THREADS}
}

# stage_00 # downloading data

stage_01
stage_02
stage_03
stage_04 100 # adjust for number of threads available
