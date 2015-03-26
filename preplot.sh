clear
echo "#####################################################################"
echo "#                                                                   #"
echo "#                     PREPLOT CONVERTOR                             #" 
echo "#                                                                   #"
echo "#####################################################################"
echo
echo "Clearing data from preplot folder..."
rm ./output/preplot/*.plt
echo
cd output
preplotRun="/usr/local/tec360/bin/preplot"
for f in *.dat
do
   echo "Processing $f"
   $preplotRun $f 
done
mv *.plt ./preplot/ # This is the name of the directory you create in 'output'
echo "#####################################################################"
echo "#                                                                   #"
echo "#                             FIN.                                  #"
echo "#                                                                   #"
echo "#####################################################################"

