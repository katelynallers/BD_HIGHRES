#! /bin/csh -f
set objdir='Data/2M1126/order*'
set spt=16.5
#get teff from Filippazzo relation
set teff=`echo "4747. - (700.5*$spt) + (115.5*($spt^2)) - (11.91*($spt^3)) + (0.6318*($spt^4)) - (0.01606*($spt^5)) + (0.0001546*($spt^6))" | bc`
echo "Teff is: $teff K"
#determine tlow and thigh to be teff -+ 500
set tlow=`echo "$teff - 500.0" | bc`
set thigh=`echo "$teff + 500.0" | bc`

echo " Finding files in " $objdir
set infiles=`find $objdir -name "*order3*.fits"`
set nfiles=`echo $infiles | wc -w`
echo " "
echo " $nfiles files found ..."

@ i = 1
while( $i <= $nfiles )
     echo "Running gtmap3 on $infiles[$i] .... "
     ./gtmap3.py ./$infiles[$i] --tlow $tlow --thigh $thigh
     echo "Running dreamzs on $infiles[$i] .... "
     ./dreamzsBD8.py ./$infiles[$i] --niter 3000 --nchain 3 --tlow $tlow --thigh $thigh
     @ i++
end
