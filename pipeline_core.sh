make_flanks() {
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN="$3" # ./stages/stage_01/Credible_TSSclusters_regions.bed
  OUTPUT_FN_PREFIX="$4" # ./stages/stage_02/TSS_clusters
  ruby make_flanks.rb ${FLANK_5} ${FLANK_3} ${INPUT_FN} > ${OUTPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.bed
}

flanks_fasta() {
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_02/TSS_clusters
  OUTPUT_FN_PREFIX="$4" # ./stages/stage_03/TSS_clusters
  make_fasta "${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.bed" \
             "${OUTPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa"
}

make_fasta() {
  INPUT_FN="$1" # ./stages/stage_02/TSS_clusters.bed
  OUTPUT_FN="$2" # ./stages/stage_03/TSS_clusters.fa
  bedtools getfasta -bed "${INPUT_FN}" \
                    -fi source_data/genome/hg38.fa \
                    -name+ \
                    -s \
      > "${OUTPUT_FN}"
}

motif_occupancies_flanks_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_03/TSS_clusters
  OUTPUT_BN_PREFIX="$4" # TSS
  motif_occupancies_cmd "${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa" \
                        "${OUTPUT_BN_PREFIX}_${FLANK_5}_${FLANK_3}.bed"
}

motif_occupancies_cmd(){
  INPUT_FN="$1" # ./stages/stage_03/TSS_clusters.fa
  OUTPUT_BN="$2" # TSS.bed
  for MOTIF_FN in $( find source_data/motifs/pfm/ -xtype f -iname '*.pfm' ); do
      MOTIF_BN=$(basename -s .pfm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ${INPUT_FN}"  \
              "${MOTIF_FN}" \
              pfm-sum-occupancy \
              --pfm-pseudocount 0.0001 \
              --naive \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb" \
          " > ./stages/stage_04/occupancy@${MOTIF_BN}@${OUTPUT_BN}"
  done
}


motif_besthits_logpval_flanks_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_03/TSS_clusters
  OUTPUT_BN_PREFIX="$4" # TSS
  motif_besthits_logpval_cmd "${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa" \
                             "${OUTPUT_BN_PREFIX}_${FLANK_5}_${FLANK_3}.bed"
}

motif_besthits_logpval_cmd(){
  INPUT_FN="$1" # ./stages/stage_03/TSS_clusters.fa
  OUTPUT_BN="$2" # TSS.bed
  for MOTIF_FN in $( find source_data/motifs/pwm/ -xtype f -iname '*.pwm' ); do
      MOTIF_BN=$(basename -s .pwm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ${INPUT_FN}"  \
              "${MOTIF_FN}" \
              besthit \
              --pvalues-file "source_data/motifs/thresholds/${MOTIF_BN}.thr" \
              --output-scoring-mode logpvalue \
              --add-flanks \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb " \
          " > ./stages/stage_04/besthit-logpval@${MOTIF_BN}@${OUTPUT_BN}"
  done
}

motif_besthits_score_flanks_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_03/TSS_clusters
  OUTPUT_BN_PREFIX="$4" # TSS
  motif_besthits_score_cmd "${INPUT_FN_PREFIX}_${FLANK_5}_${FLANK_3}_around_center.fa" \
                           "${OUTPUT_BN_PREFIX}_${FLANK_5}_${FLANK_3}.bed"
}

motif_besthits_score_cmd(){
  INPUT_FN="$1" # ./stages/stage_03/TSS_clusters.fa
  OUTPUT_BN="$2" # TSS.bed
  for MOTIF_FN in $( find source_data/motifs/pwm/ -xtype f -iname '*.pwm' ); do
      MOTIF_BN=$(basename -s .pwm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ${INPUT_FN}"  \
              "${MOTIF_FN}" \
              besthit \
              --output-scoring-mode score \
              --add-flanks \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb" \
          " > ./stages/stage_04/besthit-score@${MOTIF_BN}@${OUTPUT_BN}"
  done
}

motif_besthits_flanks_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  INPUT_FN_PREFIX="$3" # ./stages/stage_03/TSS_clusters
  OUTPUT_BN_PREFIX="$4" # TSS
  motif_besthits_score_flanks_cmd ${FLANK_5} ${FLANK_3} "${INPUT_FN_PREFIX}" "${OUTPUT_BN_PREFIX}"
  motif_besthits_logpval_flanks_cmd ${FLANK_5} ${FLANK_3} "${INPUT_FN_PREFIX}" "${OUTPUT_BN_PREFIX}"
}

motif_besthits_cmd(){
  INPUT_FN="$1" # ./stages/stage_03/TSS_clusters.fa
  OUTPUT_BN="$2" # TSS.bed
  motif_besthits_score_cmd "${INPUT_FN}" "${OUTPUT_BN}"
  motif_besthits_logpval_cmd "${INPUT_FN}" "${OUTPUT_BN}"
}
