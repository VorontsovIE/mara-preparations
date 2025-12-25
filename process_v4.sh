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
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_annotation.jsonl
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pwm.tar.gz
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pcm.tar.gz
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_pfm.tar.gz
      wget https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE-CLUSTERED/H14CORE-CLUSTERED_thresholds.tar.gz
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
  mkdir -p ./stages_fantomus/stage_01.preprocess/
  mkdir -p ./stages_fantomus/stage_01/

  ln -s ../../source_data/data-v4/CAGE_clusters_FANTOMUS_annotated.bed stages_fantomus/stage_01/
}

stage_02() {
  mkdir -p ./stages_fantomus/stage_02/
  # 5-th column used for peak summit
  make_flanks 250u 10d  <( cat ./stages_fantomus/stage_01/CAGE_clusters_FANTOMUS_annotated.bed | awk -F $'\t' -e '((NR != 1) && ($7 == "TSS")) { print $1 "\t" $2 "\t" $3 "\t" $4 "@" $9 "\t" $5 "\t" $6 }' ) ./stages_fantomus/stage_02/CAGE_clusters_FANTOMUS_annotated  5
}

stage_03() {
  mkdir -p ./stages_fantomus/stage_03/
  flanks_fasta 250u 10d ./stages_fantomus/stage_02/CAGE_clusters_FANTOMUS_annotated  ./stages_fantomus/stage_03/CAGE_clusters_FANTOMUS_annotated.fa
}

stage_04() {
  NUM_THREADS=${1:-10} # 10 threads by default

  mkdir -p ./stages_fantomus/stage_04/
  (
    motif_occupancies_flanks_cmd 250u 10d  ./stages_fantomus/stage_03/CAGE_clusters_FANTOMUS_annotated.fa  CAGE_clusters_FANTOMUS_annotated.bed ./stages_fantomus/stage_04
    # motif_besthits_flanks_cmd 250u 10d  ./stages_fantomus/stage_03/CAGE_clusters_FANTOMUS_annotated.fa  CAGE_clusters_FANTOMUS_annotated.bed ./stages_fantomus/stage_04
  ) | parallel -j ${NUM_THREADS}
}

# stage_00 # downloading data

stage_01
stage_02
stage_03
stage_04 100 # adjust for number of threads available
