#! /bin/csh -f
set objdir='Data/Vos_SpitzerC14/wl_redo/'
#set objdir='Data/GEMINI2017B/'
#set spt=18.0
#get teff from Filippazzo relation
#set teff=`echo "4747. - (700.5*$spt) + (115.5*($spt^2)) - (11.91*($spt^3)) + (0.6318*($spt^4)) - (0.01606*($spt^5)) + (0.0001546*($spt^6))" | bc`
#echo "Teff is: $teff K"
#determine tlow and thigh to be teff -+ 500
#set tlow=`echo "$teff - 500.0" | bc`
#set thigh=`echo "$teff + 500.0" | bc`

echo " Finding files in " $objdir
set infiles=`find $objdir -name "2M1425-36_33_amoeba.fits"`
set nfiles=`echo $infiles | wc -w`
echo " "
echo " $nfiles files found ..."

@ i = 1
while( $i <= $nfiles )
#while( $i <= 1 )

#grep the prior results file
    echo "Gathering header info for $infiles[$i]"
    #for nirspec
    set target=`listhead $infiles[$i] | grep -i "TARGNAME=" | head -1 | awk '{print $2}' | sed s/\'//g`
    #for gnirs
    #set target=`listhead $infiles[$i] | grep -i "OBJECT  =" | head -1 | awk '{print $3}' | sed s/\'//g`
    set order=`listhead $infiles[$i] | grep -i "ORDERS  =" | head -1 | awk '{print $3}' | sed s/\'//g`
    echo "$target-Results/$target-dreamZS-btsettl-4-$order-report.txt"
    set sysunc=`grep -i "System uncer. =" $target-Results/$target-dreamZPT-btsettl-4-$order-report.txt | head -1 | awk '{print $4}' | sed s/\'//g`
    #set sysunc=`grep -i "System uncer. =" $target-Results/$target-dreamZS-btsettl-4-$order-report.txt | head -1 | awk '{print $4}' | sed s/\'//g`
    echo "Systematic Uncertainty:  $sysunc"
    echo "Running dreamzs on $infiles[$i] .... "
#    ./dreamzsBD8.py ./$infiles[$i] --niter 3000 --nchain 3 --pton --fsys $sysunc
    ./dreamzsBD8.py ./$infiles[$i] --niter 10000 --nchain 3 --pton --fsys $sysunc
    echo "creating plots"
    ./makePlots8.py ./$infiles[$i] --fsys
    ./makePlots8.py ./$infiles[$i]
    @ i++
end
