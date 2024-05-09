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
cat ./stages/stage_01/Credible_TSSclusters_regions.bed | awktab -e '
    ($6 == "+") { print $1 "\t" ($5 - 250) "\t" ($5 + 51)  "\t" $4 "\t" "." "\t" "+"};
    ($6 == "-") { print $1 "\t" ($5 - 50)  "\t" ($5 + 251) "\t" $4 "\t" "." "\t" "-" };
  ' > ./stages/stage_02/TSS_clusters_250_50_around_center.bed

mkdir -p ./stages/stage_03/
bedtools getfasta -bed ./stages/stage_02/TSS_clusters_250_50_around_center.bed \
                  -fi source_data/genome/hg38.fa \
                  -name+ \
                  -s \
    > ./stages/stage_03/TSS_clusters_250_50_around_center.fa

mkdir -p ./stages/stage_04/
for MOTIF_FN in $( find source_data/motifs/pfm/ -xtype f ); do
    MOTIF_BN=$(basename -s .pfm "${MOTIF_FN}" )
    echo "java -cp app/sarus.jar ru.autosome.SARUS ./stages/stage_03/TSS_clusters_250_50_around_center.fa"  \
            "${MOTIF_FN}" \
            pfm-sum-occupancy \
            --pfm-pseudocount 0.0001 \
            --naive \
        " | ruby app/convert_sarus_scoresFasta_to_tsv.rb " \
        " > ./stages/stage_04/occupancy@${MOTIF_BN}.tsv "
done | parallel
