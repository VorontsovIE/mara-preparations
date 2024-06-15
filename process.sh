make_flanks() {
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN="$3" # ./stages/stage_01/Credible_TSSclusters_regions.bed
  OUTPUT_FN_PREFIX="$4" # ./stages/stage_02/TSS_clusters
  cat ${INPUT_FN} | awk -F $'\t' -e "
      (\$6 == \"+\") { print \$1 \"\t\" (\$5 - ${FLANK_5}) \"\t\" (\$5 + ${FLANK_3} + 1)  \"\t\" \$4 \"\t\" \".\" \"\t\" \"+\" };
      (\$6 == \"-\") { print \$1 \"\t\" (\$5 - ${FLANK_3})  \"\t\" (\$5 + ${FLANK_5} + 1) \"\t\" \$4 \"\t\" \".\" \"\t\" \"-\" };
    " > ${OUTPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.bed
}

flanks_fasta() {
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_02/TSS_clusters
  OUTPUT_FN_PREFIX="$4" # ./stages/stage_03/TSS_clusters
  bedtools getfasta -bed ${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.bed \
                    -fi source_data/genome/hg38.fa \
                    -name+ \
                    -s \
      > ${OUTPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa
}

motif_occupancies_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_03/TSS_clusters
  OUTPUT_BN_PREFIX="$4" # TSS
  for MOTIF_FN in $( find source_data/motifs/pfm/ -xtype f -iname '*.pfm' ); do
      MOTIF_BN=$(basename -s .pfm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa"  \
              "${MOTIF_FN}" \
              pfm-sum-occupancy \
              --pfm-pseudocount 0.0001 \
              --naive \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb " \
          " > ./stages/stage_04/occupancy@${MOTIF_BN}@${OUTPUT_BN_PREFIX}_${FLANK_5}_${FLANK_3}.bed "
  done
}

motif_besthits_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_03/TSS_clusters
  OUTPUT_BN_PREFIX="$4" # TSS
  for MOTIF_FN in $( find source_data/motifs/pwm/ -xtype f -iname '*.pwm' ); do
      MOTIF_BN=$(basename -s .pwm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa"  \
              "${MOTIF_FN}" \
              besthit \
              --pvalues-file "source_data/motifs/thresholds/${MOTIF_BN}.thr" \
              --output-scoring-mode logpvalue \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb " \
          " > ./stages/stage_04/besthit@${MOTIF_BN}@${OUTPUT_BN_PREFIX}_${FLANK_5}_${FLANK_3}.bed "
  done
}

stage_00() {
  pushd source_data/genome
      wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      gzip -d hg38.fa.gz
  popd

  pushd source_data/motifs
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_annotation.jsonl
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pwm.tar.gz
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pcm.tar.gz
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_pfm.tar.gz
      wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12INVIVO/H12INVIVO_thresholds.tar.gz
      for FN in $(ls *.tar.gz); do tar -zxf "${FN}"; done
  popd

  pushd app
      wget https://raw.githubusercontent.com/autosome-ru/sarus/master/releases/sarus-2.1.0.jar
      ln -s sarus-2.1.0.jar sarus.jar
  popd
}

stage_01() {
  mkdir -p ./stages/stage_01.preprocess/
  mkdir -p ./stages/stage_01/

  ln -s ../../source_data/Credible_TSSclusters_regions.bed stages/stage_01/
  ln -s ../../source_data/Credible_TSSclusters_TMMnormalized_logCPM.tsv stages/stage_01/

  ln -s ../../source_data/hg38_promoters_v1.bed stages/stage_01.preprocess/

  # chr1    65418   65419   promoter        65418   +   ——→   hg38_v1_chr1_+_65419_65419
  cat stages/stage_01.preprocess/hg38_promoters_v1.bed \
    | awk -F $'\t' -e '{NAME="hg38_v1_"$1"_"$6"_"($2+1)"_"$3; print $1 "\t" $2 "\t" $3 "\t" NAME "\t" $5 "\t" $6}' \
    > stages/stage_01/hg38_promoters_v1.bed
}

stage_02() {
  mkdir -p ./stages/stage_02/
  make_flanks 250 50  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 400 100 ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 300 100 ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 200 50  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters

  make_flanks 300 50  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 400 50  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters
  make_flanks 250 25  ./stages/stage_01/Credible_TSSclusters_regions.bed  ./stages/stage_02/TSS_clusters

  make_flanks 250 50  ./stages/stage_01/hg38_promoters_v1.bed ./stages/stage_02/hg38_promoters_v1
}

stage_03() {
  mkdir -p ./stages/stage_03/
  flanks_fasta 250 50  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 400 100 ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 300 100 ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 200 50  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters

  flanks_fasta 300 50  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 400 50  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters
  flanks_fasta 250 25  ./stages/stage_02/TSS_clusters  ./stages/stage_03/TSS_clusters

  flanks_fasta 250 50  ./stages/stage_02/hg38_promoters_v1  ./stages/stage_03/hg38_promoters_v1
}

stage_04() {
  NUM_THREADS=${1:-10} # 10 threads by default

  mkdir -p ./stages/stage_04/
  (
    motif_occupancies_cmd 250 50  ./stages/stage_03/TSS_clusters  TSS
    motif_occupancies_cmd 400 100 ./stages/stage_03/TSS_clusters  TSS
    motif_occupancies_cmd 300 100 ./stages/stage_03/TSS_clusters  TSS
    motif_occupancies_cmd 200 50  ./stages/stage_03/TSS_clusters  TSS

    motif_occupancies_cmd 250 50  ./stages/stage_03/hg38_promoters_v1  hg38_promoters_v1
  ) | parallel -j ${NUM_THREADS}

  (
    motif_occupancies_cmd 300 50  ./stages/stage_03/TSS_clusters  TSS
    motif_occupancies_cmd 400 50  ./stages/stage_03/TSS_clusters  TSS
    motif_occupancies_cmd 250 25  ./stages/stage_03/TSS_clusters  TSS
  ) | parallel -j ${NUM_THREADS}

  motif_besthits_cmd 250 50  ./stages/stage_03/TSS_clusters  TSS | parallel -j ${NUM_THREADS}

  (
   motif_occupancies_cmd 250 50  ./stages/stage_03/hg38_promoters_v1  hg38_promoters_v1
   motif_besthits_cmd 250 50  ./stages/stage_03/hg38_promoters_v1  hg38_promoters_v1
  ) | parallel -j ${NUM_THREADS}
}

stage_00 # downloading data

stage_01
stage_02
stage_03
stage_04 175 # adjust for number of threads available

