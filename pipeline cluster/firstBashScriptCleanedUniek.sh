#!bin/bash/
#author: Sanne van der Stam

rm loopPossibilities

#maak variabelen aan die je mee geeft tijdens het aanroepen van een script
while getopts "a:d:c:i:u:g:" options; do
  case $options in
    a)
      aantalPiRNA=$OPTARG;;
    d)
      distance=$OPTARG;;
    c)
      cutoff=$OPTARG;;
    i)
      inputMulti=$OPTARG;;
    u)
      inputUnique=$OPTARG;;
    g)
      genomeFile=$OPTARG;;
    esac
done

##pipeline1: hou de reads die tussen de 23 en 32 nt lang zijn over zodat je weet dat dit piRNAs zijn en verander de
##lengte van de reads naar 1 nt
#kijk hoeveel nt lang deze read is.
awk 'BEGIN {OFS = "\t"}  { $7 = $3 - $2 }1' $inputMulti|

#verwijder de reads buiten 23 en 32 nt
awk '(NR>0) && (($7>=23) && ($7 <=32))' |

#kijk of de strand + of - is
awk '{ if($6=="+"){ $3= $2 + 1; print $1"\t"$2"\t"$3"\t"$6;} else{$2=$3-1; print $1"\t"$2"\t"$3"\t"$6}}' |

#Sorteer en verwijder de laatste regel aangezien dit een witregel is
sort -k1,1 -k2,2n | head -n -1 > piRNAFile 

##pipeline1Uniek: hou de reads die tussen de 23 en 32 nt lang zijn over zodat je weet dat dit piRNAs zijn en 
##verander de lengte van de reads naar 1 nt. doe dit alleen voor de unieke piRNAs
#kijk hoeveel nt lang deze read is.  
awk 'BEGIN {OFS = "\t"}  { $7 = $3 - $2 }1' $inputUnique| 

#verwijder de reads buiten 23 en 32 nt
awk '(NR>0) && (($7>=23) && ($7 <=32))' |

#kijk of de strand + of - is
awk '{ if($6=="+"){ $3= $2 + 1; print $1"\t"$2"\t"$3"\t"$6;} else{$2=$3-1; print $1"\t"$2"\t"$3"\t"$6}}' | 

#Sorteer en verwijder de laatste regel aangezien dit een witregel is
sort -k1,1 -k2,2n | head -n -1 > piRNAFileUnique

##pipeline2: maak windows van de genomeFile
#maak windows
bedtools makewindows -g $genomeFile -w 5000 | 

#Sorteer
sort -k 1,1 -k 2,2n > pipeline2

##pipeline3: kijk hoeveel piRNAs in een window liggen, assign vervolgens een cutoff voor het aantal piRNAs dat in een
##window moet liggen om te bepalen welke windows meegenomen worden. merge vervolgens de windows die in een binnen een
##bepaalde afstand van elkaar liggen.
#kijk hoeveel piRNAs overlappen in een window, minstens 51% moet in een window liggen om het mee te laten tellen
bedtools coverage -F 0.51 -a pipeline2 -b piRNAFile -counts | 

#filter de windows en hou de windows over die een minimaal aantal piRNAs hebben
#maak met -v een nieuwe variabele omdat awk stom doet wanneer je daar $aantalPiRNA op schrijft
awk -v a=$aantalPiRNA '(NR>0) && ($4>=a)' |

#merge de windows die overlappen of maar maximaal 20000nt tussen zit
mergeBed -d $distance -c 4 -o count,sum > pipeline3.$aantalPiRNA.$distance

#kijk hoeveel piRNAs overlappen in een window, minstens 51% moet in een window liggen om het mee te laten
#tellen. doe dit zodat je ook de piRNAs mee neemt die wel in het window liggen maar de cutoff niet 
#gehaald hebben
bedtools coverage -F 0.51 -a pipeline3.$aantalPiRNA.$distance -b piRNAFile -counts |

#verwijder de colomn die alleen het aantal piRNAs aangeeft die de cutoff gehaald hebben
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > pipeline3.1.$aantalPiRNA.$distance

##assign een cutoff voor het aantal unieke piRNAs dat minimaan in een cluster moet zitten om de clusters mee te nemen
#kijk hoeveel unieke piRNAs overlappen in een window, minstens 51% moet in een window liggen om het mee 
#te laten tellen. doe dit om de unieke piRNAs te vinden na de cutoff en het mergen van de windows
bedtools coverage -F 0.51 -a pipeline3.1.$aantalPiRNA.$distance -b piRNAFileUnique -counts |

#geef een cutoff aan voor het aantal unieke piRNAs die er minimaal in moeten zitten
awk -v c=$cutoff '{if($6>=c){print $0}}' > pipeline3.Uniek.$aantalPiRNA.$distance.$cutoff

##pipeline4: zoek de startpositie van het cluster in het window en bepaal op bazis van de + en - strand welke kolom de
##startpositie is
#zoek de startpositie in de window
closestBed -t first -a pipeline3.Uniek.$aantalPiRNA.$distance.$cutoff -b piRNAFile  |

#verwijder de onnodige info 
awk '{print ($1"\t"$4"\t"$5"\t"$6"\t"$8"\t"$9"\t"$10)}' | 

#wanneer de strand = "+", $start = $4
#wanneer de strand = "-", $start = $5
awk '{if($7 == "+"){print $1"\t"$2"\t"$3"\t"$4"\t"$5} else{print $1"\t"$2"\t"$3"\t"$4"\t"$6}}' > pipeline4.$aantalPiRNA.$distance.$cutoff

##pipeline5: zoek de stoppositie van het cluster in het window en bepaal op bazis van de + en - strand welke kolom de
##stoppositie is
#zoek de stoppositie in de window
closestBed -t last -a pipeline3.Uniek.$aantalPiRNA.$distance.$cutoff -b piRNAFile | 

#verwijder de onnodige info
awk '{print ($1"\t"$4"\t"$5"\t"$6"\t"$8"\t"$9"\t"$10)}' | 

#wanneer de strand = "+", $stop = $5
#wanneer de strand = "-", $stop = $6
awk '{if($7 == "+"){print $1"\t"$2"\t"$3"\t"$4"\t"$6} else{print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' > pipeline5.$aantalPiRNA.$distance.$cutoff

##pipeline6: zet de file in goede volgorde met chr, start, stop, rest en bepaal de grootte van het cluster
#voeg start en stop positie File samen
paste pipeline4.$aantalPiRNA.$distance.$cutoff pipeline5.$aantalPiRNA.$distance.$cutoff | 

#filter de onnodige informatie en zet de file in goede volgorde
awk '{ print $1"\t"$5"\t"$10"\t"$2"\t"$3"\t"$4}' |

#bepaal de grootte van het cluster
awk '{print $1"\t"$2"\t"$3"\t"($3-$2),"\t"$4"\t"$5"\t"$6}' |

#maak alle clustergroottes positie, verander clustergroottes die 1nt lang zijn van -1 naar 0 
awk '{if($4=="-1"){print $1"\t"$2"\t"($3+1),"\t"($4=0),"\t"$5"\t"$6"\t"$7} else{print $0}}' |

#sorteer op aantal piRNA 
sort -k6 -nr > clusterFile.$aantalPiRNA.$distance.$cutoff 

##bereken de info die nosig is voor de loop file: totaalAantalPiRNA, totaalAantalCluster, AantalPiRNAInCluster,
##%piRNAInCluster en %GenomicSpace
#maak een file met alleen de chromosoom naam en de lengte van het cluster om daarmee de genomicSpace te berekenen
awk '{print $1 "\t" $4}' clusterFile.$aantalPiRNA.$distance.$cutoff > genomicSpace.$aantalPiRNA.$distance.$cutoff

#tel kolom 2 op om totale lengte van alle chromosomen te berekenen
awk '{sum+=$2} END {print sum}' genomeFile > sumGenomeFile

#tel kolom 2 op om de totale lengte van alle clusters op de chromosomen te berekenen
awk '{sum+=$2} END {print sum}' genomicSpace.$aantalPiRNA.$distance.$cutoff > sumGenomicSpace.$aantalPiRNA.$distance.$cutoff

#tel het aantal rijen in pipeline1 om zo te bereken hoeveel piRNAs zich op de chromosomen bevinden in de ruwe data
totaalAantalPiRNA=$(wc -l < piRNAFile)

#tel het aantal rijen -> dit is gelijk aan het aantal clusters
totaalAantalClusters=$(wc -l < clusterFile.$aantalPiRNA.$distance.$cutoff)

#bereken het aantal piRNAs in een cluster door colom 6 van pipeline6 op te tellen
piRNAInCluster=$(awk '{sum+=$6} END {print sum}' clusterFile.$aantalPiRNA.$distance.$cutoff)

#bereken het percentage piRNAs dat zich in clusters bevind
percentagePiRNA=$(awk "BEGIN {print ($piRNAInCluster / $totaalAantalPiRNA)*100}")

#bereken het percentage genomicSpace -> het percentage die de cluster innemen van een chromosoom
paste sumGenomicSpace.$aantalPiRNA.$distance.$cutoff sumGenomeFile > sumGenomeFileGenomicSpace

percentageGenomicSpace=$(awk '{print ($1/$2)*100}' sumGenomeFileGenomicSpace)

#maak een bestand waarin het aantalPiRNA & de gemergde afstand, het totaal aantal piRNAs het totaal aantal clusters
#het aantal piRNAs in een cluster, het percentage piRNA in een cluster en het percentage genomicSpace in opgeslagen
#wordt
echo "$aantalPiRNA,$distance,$cutoff	$totaalAantalPiRNA	$totaalAantalClusters	$piRNAInCluster	$percentagePiRNA	$percentageGenomicSpace" >> loopPossibilities

#verwijder de files die niet nodig zijn als output maar alleen als tussenstap
rm pipeline2
rm pipeline3.50.50000
rm pipeline3.1.50.50000
rm pipeline3.Uniek.50.50000.3
rm pipeline4.50.50000.3
rm genomicSpace.50.50000.3
rm pipeline5.50.50000.3
rm sumGenomeFile
rm sumGenomeFileGenomicSpace
rm sumGenomicSpace.50.50000.3 

#print klaar om bij te houden hoevaak het script al gerunt heeft 
echo "klaar"
