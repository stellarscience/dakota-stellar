#!/bin/csh

#################################################################
#---- Named APPScomHC in the old APPSPACK code.
#################################################################

#set HERE = (/home/tgkolda/community/interface)
set HERE = (.)

####create the names of the files unique to this iterate's function evalulation
set out = A3-$3.out   #output file
set hed = A3-$3.hed   #head file
set wel = A3-$3.wel   #well file 
set mfn = A3-$3.mfn   #mfn file 

###STEP 1: Write the .mfn file
/bin/touch $HERE/$mfn

echo " LIST  26 $out" >> $HERE/$mfn     #output file here
echo " BAS    1 A3.bas" >> $HERE/$mfn
echo " BCF   11 A3.bcf" >> $HERE/$mfn 
echo " OC    10 A3.oc" >> $HERE/$mfn 
echo "DATA  30 $hed" >> $HERE/$mfn      #heads file goes here
echo " PCG   12 A3.pcg" >> $HERE/$mfn 
echo "WEL   13 $wel" >> $HERE/$mfn      #well file goes here
echo " RCH   20 A3.rch" >> $HERE/$mfn

####STEP 2: Calculate the function value by calling Katie's code
set input = input.$3
set output = output.$3
/bin/cp $1 $HERE/$input
/bin/touch $HERE/$output

cd $HERE
./evalHC $input $mfn $wel $hed $output >& /dev/null
#./evalHC $input $mfn $wel $hed $output 

####STEP 3: Clean up and pass the result back to HOPSPACK
/bin/rm -f $HERE/$input
/bin/rm -f $HERE/$mfn
/bin/rm -f $HERE/$wel
/bin/rm -f $HERE/$hed
/bin/rm -f $HERE/$out

if(! -z $HERE/$output)then
  /bin/cp $HERE/$output $2 
endif

/bin/rm -f $HERE/$output
