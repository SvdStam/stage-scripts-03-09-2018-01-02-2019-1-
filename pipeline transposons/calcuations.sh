#!bin/bash/
#author: Sanne van der Stam

##filter de file
#verwijder alle info behalve waarvan de database repeatmasker is
awk '{if ($2 == "repeatmasker") {print $0}}' Aedes-aegypti-LVP_AGWG_REPEATFEATURES_AaegL5.gff3 | 

#split de colom met Class en Name in 2 kolommen
awk -F";" '$1=$1' OFS="\t" | 

#verwijder alle classes die geen transposons zijn
awk '$10 !~ /Unknown|Satellite|Simple_repeat|Low_complecity|tRNA|rRNA|ARTEFECT/ {print $0}' |

#verwijder onnodige info en zet de file op de goede volgorde
awk '{print $1"\t"$4"\t"$5"\t"$9"\t"$10"\t"$7}' > transposonLocationsOrderStrand

#merge alle overlap samen 
mergeBed -i transposonLocationsOrderStrand > transposonLocationsOrderStrandMerged

########################################################################################################################

##percentage genomic space uitrekenen
##tegen hoeveel procent van het genoom worden transposons gemapt

#bereken de lengte van het genoom
lengthGenome=$(awk '{sum += $2} END {print sum}' genomeFile) 

#bereken de lengte van de transposons
lengthTransposons=$(awk '{print $1"\t"$2"\t"$3"\t"($4=$3-$2)}' transposonLocationsOrderStrandMerged |

#bereken de lengte van alle transposons samen
awk '{sum += $4} END {print sum}') 

#bereken het percentage genomic space
percentageGenomicSpace=$(awk "BEGIN {print($lengthTransposons/$lengthGenome)*100}")
echo "Genomic Space:" $percentageGenomicSpace

#######################################################################################################################

##bereken het percentage clusters die bestaan uit transposons

#bereken de lengte van de clusters
lengthCluster=$(awk '{sum += $4} END {print sum}' clusterFile.50.5000.3)

#verwijder alle spaties in de clusterFile
sed 's/ //g' clusterFile.50.5000.3 > clusterFile.50.5000.3NoSpaces

#map de transposons tegen de clusters
countBPTransposons=$(bedtools coverage -a clusterFile.50.5000.3NoSpaces -b transposonLocationsOrderStrandMerged |

#tel kolom 9 (het aantal base paren transposons) op
awk '{sum += $9} END {print sum}')

#bereken het percentage clusters die bestaan uit transposons
percentageClustersUitTransposons=$(awk "BEGIN {print($countBPTransposons/$lengthCluster)*100}")
echo "Cluster uit TE:" $percentageClustersUitTransposons

#######################################################################################################################

##bereken het percentage piRNAs die gemapt worden tegen transposons

#bereken de lengte van alle piRNAs samen
lengthPiRNA=$(wc -l piRNAFile | awk '{print $1}')

#map de piRNAs tegen de transposons
piRNAsMapTransposons=$(bedtools coverage -a transposonLocationsOrderStrandMerged -b piRNAFile |

#tel kolom 4 (het aantal piRNAs) op
awk '{sum += $4} END {print sum}')

#bereken het percentage piRNAs die tegen transposons gemapt worden
percentagePiRNAsTegenTransposons=$(awk "BEGIN {print($piRNAsMapTransposons/$lengthPiRNA)*100}")
echo "piRNA tegen TE:" $percentagePiRNAsTegenTransposons

#######################################################################################################################
##37% van de piRNAs worden gemapt tegen transposons
##hoeveel procent hiervan valt in clusters?

#vind overlap tussen piRNAs en transposons
bedtools intersect -a transposonLocationsOrderStrandMerged -b piRNAFile > piRNAsOverlapMetTransposons

#map de overlap van de piRNAs & transposons tegen de cluters
aantalPiRNATegenTransposonsInCluster=$(bedtools coverage -a clusterFile.50.5000.3NoSpaces -b piRNAsOverlapMetTransposons |

#tel $8 (aantal transposons) op
awk '{sum += $8} END {print sum}')

#bereken het percentage piRNAs die tegen transposons mappen en in clusters vallen
percentagePiRNAsTegenTransposonsInClusters=$(awk "BEGIN {print($aantalPiRNATegenTransposonsInCluster/$piRNAsMapTransposons)*100}")
echo "piRNA tegen TE in cluster:" $percentagePiRNAsTegenTransposonsInClusters

rm transposonLocationsOrderStrand
rm piRNAsOverlapMetTransposons
rm clusterFile.50.5000.3NoSpaces
rm transposonLocationsOrderStrandMerged

echo "Done" 
