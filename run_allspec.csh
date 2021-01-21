#! /bin/csh -f
set objdir='Data/Vos_Oct2020'
#set objdir='SIMULATED_DATA/'
echo " Finding files in " $objdir
#specific for lists
set infiles=`find $objdir -name "0642+41_GNIRS.fits"`
set nfiles=`echo $infiles | wc -w`
echo " "
echo " $nfiles files found ..."

@ i = 1
while( $i <= $nfiles )
#get teff from Filippazzo relation
#    echo $infiles[$i]
#    echo "SpT is $spt[$i]"
#    set teff=`echo "4747. - (700.5*$spt[$i]) + (115.5*($spt[$i]^2)) - (11.91*($spt[$i]^3)) + (0.6318*($spt[$i]^4)) - (0.01606*($spt[$i]^5)) + (0.0001546*($spt[$i]^6))" | bc`
#    echo "Teff is: $teff K"
#determine tlow and thigh to be teff -+ 500
#    set tlow=`echo "$teff - 500.0" | bc`
#    set thigh=`echo "$teff + 500.0" | bc`
    echo "Running gtmap3 on $infiles[$i] .... "
    ./gtmap3.py ./$infiles[$i]
    echo "Running dreamzs on $infiles[$i] .... "
    ./dreamzsBD8.py ./$infiles[$i] --niter 3000 --nchain 3 --pton
    #./dreamzsBD8.py ./$infiles[$i] --niter 3000 --nchain 3 --fsys 0.005
    echo "creating plots"
    ./makePlots8.py ./$infiles[$i]
    @ i++
end
