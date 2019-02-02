#!bin/bash/
#author: Sanne van der Stam

#maak variabelen aan die je mee geeft tijdens het aanroepen van een script
while getopts "a:e:v:" options; do
  case $options in
    a)
      cutoffLength=$OPTARG;;
    e)
      EVEFile=$OPTARG;;
    v)
      viralTax=$OPTARG;;
    esac
done

#verwijder de header
tail -n +3 $EVEFile | 

#selecteer de kolommen die je wilt gebruiken
awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }' | 

#splits op [ om het virus te krijgen die je wilt vergelijken met de andere file
awk -F"[" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }' | 

#verwijder ]
awk '{gsub(/\]/,"")}1' |

#vervang de spatie met een @ zodat er geen spaties meer in de file staat en hij alleen gescheiden word door tabs
tr " " @ |

#wissel de laatste 2 kolommen om de spatie(@) weg te halen bij kolom 7
awk '{print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' | sed 's/.$//'| 

#wissel de laatste 2 kolommen weer terug
awk '{print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}' |
sort -k1,1 | uniq > NIRVS_AaegL5VirusSeperateColumn

#######################################################################################################################

#verwijder de header
tail -n +3 $viralTax | awk 'BEGIN {FS="\t"}{print $4"\t"$5}' | 

#splits op de komma
awk -F"," '{print $1"\t"$2"\t"$3"\t"$4}' |

#verander de spatie in een @ zodat er geen spaties meer in de file staat en hij alleen gescheiden word door tabs
tr " " @ > taxid10239Column4+5

#wanneer kolom 3 en 4 geen waarde bevatten verplaats de waarde van kolom 1 naar kolom 3 en kolom 2 naae kolom 4 en vul voor de andere waardes een . in -> deze waardes zijn niet bekent
#wanneer alleen kolom 4 geen waarde bevat, verplaatst alle waardes dan een plek naar rechts en vul voor kolom 1 een . in -> deze waarde is onbekent
awk 'BEGIN{FS="\t"}{ if ($3 == "" && $4 == "") print ($3="."),"\t"($4="."),"\t"$1"\t"$2; else if ($4 == "") print ($4="."),"\t"$1"\t"$2"\t"$3; else if ($2 == "") print $1"\t"($2="."),"\t"$3"\t"$4; else print $0}' taxid10239Column4+5 |

#pak alleen alle unieke waardes uit de taxid file zodat er geen dubbele waardes in komen te staan. 
#wanneer dit niet gebeurt kan 1 waarde in de NIRVS_AaegL5 file vaker voorkomen omdat dit virus vaker voorkomt in de taxid file
uniq > taxid10239Column4+5Uniq

#######################################################################################################################

#sorteer de taxid file op kolom 4 omdat je de files later joint op deze kolom
sort -k 4,4 taxid10239Column4+5Uniq >  taxid10239Column4+5UniqSorted

#sorteer de NIRVS_AaegL5 file op kolom 8 omdat je de files later joint op deze kolom
sort -k 8,8 NIRVS_AaegL5VirusSeperateColumn > NIRVS_AaegL5VirusSeperateColumnSorted

#join de 2 files. 
#voor de taxid file doe dit op basis van kolom 4 en voor de NIRVS_AaegL5 file op kolom 8.
#print de kolommen in de volgorde: de file waarop je joint(species), family, genus, type species, Name, chr, start, stop, length, strand, blast_hit.
#voor de virussen die wel in NIRVS_AaegL5 staan maar niet in taxid word UNKNOWN ingevuld. 
join -2 8 -1 4 -o 0,1.1,1.2,1.3,2.1,2.2,2.3,2.4,2.5,2.6,2.7 -e UNKNOWN taxid10239Column4+5UniqSorted NIRVS_AaegL5VirusSeperateColumnSorted -a2|

#verander de space delimited in tab delimited
tr ' ' "\t"  > NIRVS_AaegL5Taxonomie

#zet de file in goede volgorde
awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$5"\t"$10"\t"$11"\t"$1"\t"$2"\t"$3"\t"$4}' NIRVS_AaegL5Taxonomie > NIRVS_AaegL5TaxonomieOrder

#voeg Unidentified toe wanneer de Name van de Eve Un bevat. 
awk '{if ($5 ~ /Un/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9="Unidentified"),"\t"($10="Unidentified"),"\t"$11 } else{print $0}}' NIRVS_AaegL5TaxonomieOrder |

#wanneer kolom 10 een . bevat verander kolom 10 dan in TaxonomieUnknown. dit betekent dat de 2de rang van taxonomie niet
#bekent is maar de Eve wel gevonden is in de taxid file.
awk '{if ($10=="."){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"($10="TaxonomieUnknown"),"\t"$11}else{print $0}}' |

#wanneer kolom 9 gelijk is aan UNKNOWN en kolom 5 Rha bevat maak dan van kolom 9 Rhabdoviridae
awk '{if ($9=="UNKNOWN" && $5~/Rha/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9="Unclassified@Rhabdoviridae"),"\t"$10"\t"$11}else{print $0}}' | 

#wanneer kolom 9 gelijk is aan UNKNOWN en kolom 5 Mon bevat maak dan van kolom 9 Mononegavirales
awk '{if ($9=="UNKNOWN" && $5~/Mon/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9="Unclassified@Mononegavirales"),"\t"$10"\t"$11}else{print $0}}' |  

#wanneer kolom 9 gelijk is aan TaxonomieUnknown en kolom 5 Rha bevat maak dan van kolom 9 Rhabdoviridae
awk '{if ($9=="." && $5~/Rha/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9="Unclassified@Rhabdoviridae"),"\t"$10"\t"$11}else{print $0}}'| 

#wanneer kolom 9 gelijk is aan TaxonomieUnknown en kolom 5 Rha bevat maak dan van kolom 9 Rhabdoviridae
awk '{if ($9=="." && $5~/Mon/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9="Unclassified@Mononegavirales"),"\t"$10"\t"$11}else{print $0}}'|

#wanneer kolom 9 gelijk is aan UNKNOWN en kolom 5 Bun bevat maak dan van kolom 9 Bunyavirales
awk '{if ($9=="UNKNOWN" && $5~/Bun/){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9="Unclassified@Bunyavirales"),"\t"$10"\t"$11}else{print $0}}' |

#wanneer kolom 10 een genus bevat zet deze genus dan in kolom 12
awk '{if($10 !~ /UNKNOWN|Unidentified|TaxonomieUnknown/){print $0"\t"($12="Genus:" $10)}else{print $0}}' | 

#wanneer kolom 9 een Order bevat zet deze Order dan in kolom 12
awk '{if($9 ~ /Mononegavirales|Bunyavirales/ && $12==""){print $0"\t"($12="Order:" $9)}else{print $0}}' | 

#wanneer kolom 9 een Familie bevat zet deze Familie dan in kolom 12
awk '{if($9 ~ /Rhabdoviridae|Flaviviridae/ && $12==""){print $0"\t"($12="Family:" $9)}else{print $0}}' |

#wanneer kolom 12 geen waarde bevat zet dan Unidentified in kolom 12
awk '{if($12==""){print $0"\t"($12="Unidentified")}else{print $0}}' > NIRVS_AaegL5TaxonomieOrderFamily

########################################################################################################################

##voeg een cluster nummer toe aan de clusters
sort -k 4,4nr clusterFile.50.5000.3 |

#voeg een cluster nummer toe
awk '{print NR,$0}' |

#sorteer de file en voeg een # toe voor het cluster nummer
awk '{print $2"\t"$3"\t"$4"\t"$5"\t"($1="#"$1),"\t"$6"\t"$7"\t"$8}' > clusterFile.50.5000.3ClusterNumber

########################################################################################################################

#kijk welke Eves in de clusters liggen
bedtools coverage -a clusterFile.50.5000.3ClusterNumber -b NIRVS_AaegL5TaxonomieOrder |

#sorteer op het aantal Eves in een cluster
sort -k 9,9nr > EvesInClusters

#map de virussen tegen de clusters om te kijken welke Eves voorkomen in de clusters
bedtools intersect -a EvesInClusters -b NIRVS_AaegL5TaxonomieOrderFamily -wa -wb | bedtools groupby -g 1,2,3,5 -c 24 -o freqasc |

#splits de kolommen op komma om van de laatste kolom allemaal losse kolommen te maken
awk -F"," 'BEGIN{OFS="\t"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }'  > EVENamesInCluster

#######################################################################################################################

##maak een file aan met het cluster (chr, start, stop, number), % piRNAs op de + strand, % piRNAs op de - strand,
##% EVEs op de + strand, % EVEs op de - strand, % cluster dat uit EVEs bestaat

#kijk naar de overlap tussen de clusters en de piRNAs om te bepalen hoeveel piRNAs op welke strand liggen
bedtools intersect -a piRNAFile -b clusterFile.50.5000.3ClusterNumber -wb|

#geef aan welke colommen overgehouden moeten worden en welke kolom opgeteld moet worden
bedtools groupby -g 5,6,7,8,9,11 -c 4  -o freqasc | 

#splits de colom op , zodat de kolom met - en + 2 apparte kolommen worden
awk -F"," 'BEGIN{OFS="\t"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' |

#zet de kolom met + in kolom 5 en de kolom met - in 6 voor een beter overzicht
awk '{if($7 ~ /+/){print $0}else{print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7}}' |

#splits de plus en min kolom op :
awk 'BEGIN{FS=OFS=":"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'|

#verwijder de "+" en de "-"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$10}' |

#voeg een 0 toe als de column leeg is
awk 'BEGIN{FS="\t"}{if($8==""){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"($8=0)}else{print $0}}' |

#tel de plus en min column op
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($9=($7+$8))}' | 

#bereken het percentage op de plus en min strand
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"($10=($7/$9)*100),"\t"($11=($8/$9)*100)}' > piRNAInClusterOpStrand

#######################################################################################################################

#verander de volgorde van de file
awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8}' NIRVS_AaegL5VirusSeperateColumn | 

#sorteer op chromosoom
sort -k 1,1 > NIRVS_AaegL5VirusSeperateColumnSorted

#map de EVEs tegen de clusters
bedtools intersect -a NIRVS_AaegL5VirusSeperateColumnSorted -b clusterFile.50.5000.3ClusterNumber -wb | 

#maak een file om te kijken op welke strands de Eves liggen
bedtools groupby -g 9,10,11,12,13 -c 6  -o freqasc |

#sorteer op clusternumber
sort -k 5,5 > EVEsOpClusterSorted

#sorteer de file op clusternumber
sort -k 5,5 piRNAInClusterOpStrand > piRNAInClusterOpStrandSorted

#join de 2 files en hou alles over van de piRNA file en alleen de strand van de EVE file 
join -2 5 -1 5 -o 0,2.1,2.2,2.3,2.4,2.6,2.7,2.8,2.9,2.10,2.11,1.6 -e NoEVE EVEsOpClusterSorted piRNAInClusterOpStrandSorted -a2 | 

#zet een tab tussen alle kolommen
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' | 

sort -k 5,5nr | 

awk -v a=$cutoffLength '{if($5 >= a)print $0}' > strandFile

#voeg een header toe aam de file
sed  -i '1i ClusterNumber\tChr\tStartCluster\tStopCluster\tLengthCluster\tAantalPiRNA\tAantalPiRNAopPlus\tAantalPiRNAopMin\tTotaalAantalOpStrand\tPercentagePiRNAOpPlus\tPercentagePiRNAOpMin\tEVEsOpStrand' strandFile

#######################################################################################################################

#map de EVE names tegen de clusters om te kijken welke names in een cluster liggen
bedtools intersect -a EvesInClusters -b NIRVS_AaegL5TaxonomieOrderFamily -wa -wb | 

#print alleen de laatste kolom
awk -F '\t' '{print $24}' | 

#sorteer
sort | 

#tel hoevaak alle names voorkomen
uniq -c |
 
#sorteer op aantal
sort -nr |

#verander de spatie tussen de 2 kolommen in een tab
awk -F " " '{print $1"\t"$2}' |

#sorteer op de Names
sort -k 2,2 > numberNamesInCluster

#print alleen de laatste kolom
awk -F '\t' '{print $12}' NIRVS_AaegL5TaxonomieOrderFamily | 

#sorteer
sort | 

#tel hoevaak alle names voorkomen
uniq -c | 

#sorteer op aantal
sort -nr | 

#verander de spatie tussen de 2 kolommen in een tab
awk -F " " '{print $1"\t"$2}' |

#sorteer op de Names
sort -k 2,2 > numberNamesTotal

#join de files op Name en print daarbij de Name het aantal totaal en het aantal in cluster
join -2 2 -1 2 -o 0,1.1,2.1 -e 0 numberNamesTotal numberNamesInCluster -a1 |

#trek kolom 2 van kolom 3 af om het aantal buiten een cluster te weten
awk '{print $1"\t"$2"\t"$3"\t"($4=$2-$3), "\t"($5=130)}' |

#bereken het percentage binnen en buiten en voeg "inside" toe om voor het plotten te weten dat dit over in cluster gaat
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6=($3/$5*100)),"\t"($7=($4/$5*100)), "\t"($8="Inside")}' > numberNamesTotalInOut1

#join de files op Name en print daarbij de Name het aantal totaal en het aantal in cluster
join -2 2 -1 2 -o 0,1.1,2.1 -e 0 numberNamesTotal numberNamesInCluster -a1 |

#trek kolom 2 van kolom 3 af om het aantal buiten een cluster te weten
awk '{print $1"\t"$2"\t"$3"\t"($4=$2-$3), "\t"($5=83)}' |

#bereken het percentage binnen en buiten de clusters
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($6=($3/$5*100)),"\t"($7=($4/$5*100))}' | 

#verander de volgorde zodat binnen en buiten onder elkaar staan en voeg outside toe
awk '{print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$7"\t"$6, "\t"($8="Outside")}' > numberNamesTotalInOut

#voeg inside en outside samen
cat numberNamesTotalInOut1 >> numberNamesTotalInOut

#print alleen de kolommen die nodig zijn voor de plot
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$8}' numberNamesTotalInOut > numberNamesTotalInOutR

#voeg headers toe
sed  -i '1i Name\tTotal\tInOutCluster\tTotalEves\tPercentage\tInOrOut' numberNamesTotalInOutR

rm NIRVS_AaegL5VirusSeperateColumn
rm taxid10239Column4+5
rm clusterFile.50.5000.3ClusterNumber
rm EvesInClusters
rm NIRVS_AaegL5Taxonomie
rm NIRVS_AaegL5TaxonomieOrder
rm NIRVS_AaegL5TaxonomieOrderFamily
rm taxid10239Column4+5Uniq
rm taxid10239Column4+5UniqSorted
rm EVEsOpClusterSorted
rm NIRVS_AaegL5VirusSeperateColumnSorted
rm numberNamesInCluster
rm numberNamesTotal
rm numberNamesTotalInOut
rm numberNamesTotalInOut1
rm piRNAInClusterOpStrand
rm piRNAInClusterOpStrandSorted

"Done"
