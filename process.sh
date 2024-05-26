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


mkdir -p ./stages/stage_01/
ln -s ../../source_data/Credible_TSSclusters_regions.bed stages/stage_01/
ln -s ../../source_data/Credible_TSSclusters_TMMnormalized_logCPM.tsv stages/stage_01/

mkdir -p ./stages/stage_02/
make_flanks() {
  FLANK_5=$1
  FLANK_3=$2
  cat ./stages/stage_01/Credible_TSSclusters_regions.bed | awk -F $'\t' -e "
      (\$6 == \"+\") { print \$1 \"\t\" (\$5 - ${FLANK_5}) \"\t\" (\$5 + ${FLANK_3} + 1)  \"\t\" \$4 \"\t\" \".\" \"\t\" \"+\" };
      (\$6 == \"-\") { print \$1 \"\t\" (\$5 - ${FLANK_3})  \"\t\" (\$5 + ${FLANK_5} + 1) \"\t\" \$4 \"\t\" \".\" \"\t\" \"-\" };
    " > ./stages/stage_02/TSS_clusters_${FLANK_5}_${FLANK_3}_around_center.bed
}
make_flanks 250 50
make_flanks 400 100

make_flanks 300 100
make_flanks 200 50

mkdir -p ./stages/stage_03/
flanks_fasta() {
  FLANK_5=$1
  FLANK_3=$2
  bedtools getfasta -bed ./stages/stage_02/TSS_clusters_${FLANK_5}_${FLANK_3}_around_center.bed \
                    -fi source_data/genome/hg38.fa \
                    -name+ \
                    -s \
      > ./stages/stage_03/TSS_clusters_${FLANK_5}_${FLANK_3}_around_center.fa
}
flanks_fasta 250 50
flanks_fasta 400 100

flanks_fasta 300 100
flanks_fasta 200 50

mkdir -p ./stages/stage_04/
motif_occupancies_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  for MOTIF_FN in $( find source_data/motifs/pfm/ -xtype f -iname '*.pfm' ); do
      MOTIF_BN=$(basename -s .pfm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ./stages/stage_03/TSS_clusters_${FLANK_5}_${FLANK_3}_around_center.fa"  \
              "${MOTIF_FN}" \
              pfm-sum-occupancy \
              --pfm-pseudocount 0.0001 \
              --naive \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb " \
          " > ./stages/stage_04/occupancy@${MOTIF_BN}@TSS_${FLANK_5}_${FLANK_3}.bed "
  done
}

motif_besthits_cmd(){
  FLANK_5=$1
  FLANK_3=$2
  for MOTIF_FN in $( find source_data/motifs/pwm/ -xtype f -iname '*.pwm' ); do
      MOTIF_BN=$(basename -s .pwm "${MOTIF_FN}" )
      echo "java -cp app/sarus.jar ru.autosome.SARUS ./stages/stage_03/TSS_clusters_${FLANK_5}_${FLANK_3}_around_center.fa"  \
              "${MOTIF_FN}" \
              besthit \
              --pvalues-file "source_data/motifs/thresholds/${MOTIF_BN}.thr" \
              --output-scoring-mode logpvalue \
          " | ruby app/convert_sarus_scoresFasta_to_tsv.rb " \
          " > ./stages/stage_04/besthit@${MOTIF_BN}@TSS_${FLANK_5}_${FLANK_3}.bed "
  done
}


(
  motif_occupancies_cmd 250 50
  motif_occupancies_cmd 400 100

  motif_occupancies_cmd 300 100
  motif_occupancies_cmd 200 50
) | parallel -j 30


make_flanks 300 50
make_flanks 400 50
make_flanks 250 25
flanks_fasta 300 50
flanks_fasta 400 50
flanks_fasta 250 25

(
  motif_occupancies_cmd 300 50
  motif_occupancies_cmd 400 50
  motif_occupancies_cmd 250 25
) | parallel -j 175

motif_besthits_cmd 250 50 | parallel -j 175
