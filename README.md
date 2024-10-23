# Preparation process
## Source data
Put two files into source_data folder:
* `Credible_TSSclusters_regions.bed`
* `Credible_TSSclusters_TMMnormalized_logCPM.tsv`

## Processing
See `process.sh`

* Take CAGE peaks [-250; +50] and [-400; +100] around peak center (as specified in bed file; it can be replaced with maxima).
  Coordinates are in the hg38 assembly.
* Hocomoco 12, in-vivo, full collection is used.
* Sum-occupancy (PFM pseudocount = 0.0001) is calculated for each motif, each interval.
