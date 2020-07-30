#!/bin/bash

#Set up to run in parallel across each sample directory.

for directoryname in /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1031 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-009 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-390 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-263 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-233 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-910 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-519 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1070 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-411 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-040 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-881 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-187 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1122 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-362 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1441 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-526 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1497 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-686 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-355 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1174 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1369 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-496 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-1027 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-886 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-724 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-115 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-602 /gscmnt/gc2547/griffithlab/matthewmosior/MCL_mosdepth/samples/TWGE-08-0075-930 
do
	#Change directory.
	cd $directoryname
	
	#Grab the normal file.
        normalthresholdsfile=$(find . -name "TWGE*normal*thresholds.bed.gz")
	normalcorrectedthresholdsfile=$(echo $normalthresholdsfile | sed 's/\.\///')

	#Grab the tumor file.
        tumorthresholdsfile=$(find . -name "TWGE*tumor*thresholds.bed.gz")
        tumorcorrectedthresholdsfile=$(echo $tumorthresholdsfile | sed 's/\.\///')

	#Remove .gz from normalthresholdsfile.
	newnormalthresholdsfile=$(echo $normalthresholdsfile | sed 's/.gz//' | sed 's/\.\///')
	
	#Remove .gz from tumorthresholdsfile.
        newtumorthresholdsfile=$(echo $tumorthresholdsfile | sed 's/.gz//' | sed 's/\.\///')

	#Grab everything but the head of the thresholds file.
        gunzip -c "$normalcorrectedthresholdsfile" > $newnormalthresholdsfile
	
	#Grab everything but the head of the thresholds file.
        gunzip -c "$tumorcorrectedthresholdsfile" > $newtumorthresholdsfile

	#Remove head from newthresholdsfile.
	cat $newnormalthresholdsfile | tail -n +2 > headless_$newnormalthresholdsfile

	#Remove head from newthresholdsfile.
        cat $newtumorthresholdsfile | tail -n +2 > headless_$newtumorthresholdsfile
	
	#Loop through thresholdsfile line by line (normal).
	while read line
	do
		#Set variables.
		normalstart=$(echo "$line" | cut -f 2)
		normalstop=$(echo "$line" | cut -f 3)
		normaltwenty=$(echo "$line" | cut -f 5)
		normalthirty=$(echo "$line" | cut -f 6)
		normalfifty=$(echo "$line" | cut -f 7)
		normalonehundred=$(echo "$line" | cut -f 8)
		normalonefifty=$(echo "$line" | cut -f 9)
		normaltwohundred=$(echo "$line" | cut -f 10)
		normaltwofifty=$(echo "$line" | cut -f 11)
		normalthreehundred=$(echo "$line" | cut -f 12)

		#Calculate demoninator.
		normaldenominator=$(($normalstop - $normalstart))	

		#Divide 20X, 30X, 50X, 100X, 150X, 200X, 250X, 300X by denominator.
		normaltwentyXdenom=$(echo "scale=4 ; $normaltwenty / $normaldenominator" | bc)
		normalthirtyXdenom=$(echo "scale=4 ; $normalthirty / $normaldenominator" | bc)
		normalfiftyXdenom=$(echo "scale=4 ; $normalfifty / $normaldenominator" | bc)	
		normalonehundredXdenom=$(echo "scale=4 ; $normalonehundred / $normaldenominator" | bc)
		normalonefiftyXdenom=$(echo "scale=4 ; $normalonefifty / $normaldenominator" | bc)
		normaltwohundredXdenom=$(echo "scale=4 ; $normaltwohundred / $normaldenominator" | bc)
		normaltwofiftyXdenom=$(echo "scale=4 ; $normaltwofifty / $normaldenominator" | bc)
		normalthreehundredXdenom=$(echo "scale=4 ; $normalthreehundred / $normaldenominator" | bc)


		#Print the link back with new columns added.
		echo "$line" | awk -v normaltwentyfinal="$normaltwentyXdenom" -v normalthirtyfinal="$normalthirtyXdenom" -v normalfiftyfinal="$normalfiftyXdenom" -v normalonehundredfinal="$normalonehundredXdenom" -v normalonefiftyfinal="$normalonefiftyXdenom" -v normaltwohundredfinal="$normaltwohundredXdenom" -v normaltwofiftyfinal="$normaltwofiftyXdenom" -v normalthreehundredfinal="$normalthreehundredXdenom" '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" normaltwentyfinal "\t" normalthirtyfinal "\t" normalfiftyfinal "\t" normalonehundredfinal "\t" normalonefiftyfinal "\t" normaltwohundredfinal "\t" normaltwofiftyfinal "\t" normalthreehundredfinal}' >> calculated_$newnormalthresholdsfile	
	done < headless_$newnormalthresholdsfile

	#Loop through thresholdsfile line by line (tumor).
        while read line
        do
		#Set variables.
                tumorstart=$(echo "$line" | cut -f 2)
                tumorstop=$(echo "$line" | cut -f 3)
                tumortwenty=$(echo "$line" | cut -f 5)
                tumorthirty=$(echo "$line" | cut -f 6)
                tumorfifty=$(echo "$line" | cut -f 7)
                tumoronehundred=$(echo "$line" | cut -f 8)
                tumoronefifty=$(echo "$line" | cut -f 9)
                tumortwohundred=$(echo "$line" | cut -f 10)
                tumortwofifty=$(echo "$line" | cut -f 11)
                tumorthreehundred=$(echo "$line" | cut -f 12)

                #Calculate demoninator.
                tumordenominator=$(($tumorstop - $tumorstart))

                #Divide 20X, 30X, 50X, 100X, 150X, 200X, 250X, 300X by denominator.
                tumortwentyXdenom=$(echo "scale=4 ; $tumortwenty / $tumordenominator" | bc)
                tumorthirtyXdenom=$(echo "scale=4 ; $tumorthirty / $tumordenominator" | bc)
                tumorfiftyXdenom=$(echo "scale=4 ; $tumorfifty / $tumordenominator" | bc)
                tumoronehundredXdenom=$(echo "scale=4 ; $tumoronehundred / $tumordenominator" | bc)
                tumoronefiftyXdenom=$(echo "scale=4 ; $tumoronefifty / $tumordenominator" | bc)
                tumortwohundredXdenom=$(echo "scale=4 ; $tumortwohundred / $tumordenominator" | bc)
                tumortwofiftyXdenom=$(echo "scale=4 ; $tumortwofifty / $tumordenominator" | bc)
                tumorthreehundredXdenom=$(echo "scale=4 ; $tumorthreehundred / $tumordenominator" | bc)


                #Print the link back with new columns added.
                echo "$line" | awk -v tumortwentyfinal="$tumortwentyXdenom" -v tumorthirtyfinal="$tumorthirtyXdenom" -v tumorfiftyfinal="$tumorfiftyXdenom" -v tumoronehundredfinal="$tumoronehundredXdenom" -v tumoronefiftyfinal="$tumoronefiftyXdenom" -v tumortwohundredfinal="$tumortwohundredXdenom" -v tumortwofiftyfinal="$tumortwofiftyXdenom" -v tumorthreehundredfinal="$tumorthreehundredXdenom" '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" tumortwentyfinal "\t" tumorthirtyfinal "\t" tumorfiftyfinal "\t" tumoronehundredfinal "\t" tumoronefiftyfinal "\t" tumortwohundredfinal "\t" tumortwofiftyfinal "\t" tumorthreehundredfinal}' >> calculated_$newtumorthresholdsfile
        done < headless_$newtumorthresholdsfile

	#Remove headless_$newnormalthresholdsfile.
	rm headless_$newnormalthresholdsfile
	rm $newnormalthresholdsfile

	#Remove headless_$newtumorthresholdsfile.
        rm headless_$newtumorthresholdsfile
        rm $newtumorthresholdsfile
done
