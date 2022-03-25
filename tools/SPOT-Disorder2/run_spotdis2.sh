#!/bin/bash

################################################################################
#UPDATE these to the correct installation directories or values
psiblastdir=''
spotcontactdir=/home/jack/Documents/Programs/Contact_Map_Predictors/SPOT-Contact-Helical-New
hhblitsdir=''
hhsuitedir=/home/jack/Documents/Programs/hh-suite
nrdir=/home/jack/Documents/Databases/Protein/nr/nr
uniprotdir=/home/jack/Documents/Databases/Protein/uniprot20_2013_03/uniprot20_2013_03/uniprot20_2013_03
ncpu=16
spider3dir=/home/jack/Documents/Programs/SPIDER3-numpy-server
spot1ddir=/home/jack/Documents/Programs/Secondary_Structure_Predictors/SPOT1D
spotdis2dir=/home/jack/Documents/Programs/Disorder_Predictors/SPOT-Disorder2
################################################################################

export HHLIB=$hhsuitedir
#Generate the input file list

#-------------------------------------------------------------------------------
# RUN ON ALL FASTA FILES IN INPUT
find $spotdis2dir/inputs -name "*.fasta" | sed "s:.fasta::" | sed "s:.*/::" > protlist.txt #COMMENT THIS LINE OUT IF YOU HAVE YOUR OWN PROTLIST.TXT!!
#-------------------------------------------------------------------------------
# RUN ONLY ON SPECIFIED FILES IN PROTLIST.txt
prots=`sed -e "s:^:$spotdis2dir/inputs/:" protlist.txt`
#-------------------------------------------------------------------------------
orig_dir=`pwd`
for i in $prots; do
	protname=`echo $i | sed "s:.*/\(.*\):\1:"`   
    echo "Protein: $protname"
    echo "Generating PSSM..."
	[ -f $i.pssm ] || ${psiblastdir}psiblast -db $nrdir -num_iterations 3 -num_alignments 1 -num_threads $ncpu -query $i.fasta -out  $i.bla -out_ascii_pssm $i.pssm #-out_pssm $i.chk
    echo "Generating HHblits outputs..."
	[ -f $i.hhm -o -f $i.a3m ] || ${hhblitsdir}hhblits -i $i.fasta -ohhm $i.hhm -oa3m $i.a3m -d $uniprotdir -v 0 -maxres 40000 -cpu $ncpu -Z 0
    echo "Generating SPIDER3 outputs..."
    cd inputs
    [ -f $i.spd33 -o -f $i.spotcon ] || ${spider3dir}/script/spider3_pred.py $protname --odir "$spotdis2dir/inputs/"
    cd $orig_dir
    if [[ ! -f $i.spotcon ]]; then
        cd $spotcontactdir
        echo "Generating CCMPRED outputs..."
        [ -f $i.mat -o -f $i.spotcon ] || sources/CCMpred/calCCMpred.sh $i.a3m
        [ -f $protname.mat -o -f $i.spotcon ] && mv $protname.mat $spotdis2dir/inputs/
        echo "Generating DCA outputs..."
        [ -f $i.di -o -f $i.spotcon ] || sources/DCA/calDI.sh $i.a3m
        [ -f $protname.di -o -f $i.spotcon ] && mv $protname.di $spotdis2dir/inputs/
        cd $orig_dir
    fi
done
#-------------------------------------------------------------------------------
cd $spotcontactdir
rm -rf $spotdis2dir/spotconlist
touch $spotdis2dir/spotconlist
for i in $prots; do 
    echo "Checking for SPOT-Contact files..."
    [[ -f $i.spotcon ]] || echo $i | sed "s:.*/\(.*\):\1:" >> $spotdis2dir/spotconlist
done
if [[ -s $spotdis2dir/spotconlist ]]; then
    echo "SPOT-Contact features collected! Running SPOT-Contact..."
    $spotcontactdir/scripts/run_all_models.py --input_list "$spotdis2dir/spotconlist" --gpu 0 --input_dir $spotdis2dir/
    for i in $prots; do 
	    protname=`echo $i | sed "s:.*/\(.*\):\1:"` 
        mv $spotcontactdir/outputs/$protname.spotcon $spotdis2dir/inputs/
    done
fi
#-------------------------------------------------------------------------------
cd $spot1ddir

rm -rf $spotdis2dir/spot1dlist
touch $spotdis2dir/spot1dlist
for i in $prots; do 
    echo "Checking for SPOT-1D files..."
    [[ -f $i.spot1d ]] || echo $i | sed "s:.*/\(.*\):\1:" >> $spotdis2dir/spot1dlist
done
if [[ -s $spotdis2dir/spot1dlist ]]; then
    echo "Running SPOT-1D..."
    $spot1ddir/run_all_models.py --input_list "$spotdis2dir/spot1dlist" --gpu 0 --input_dir $spotdis2dir/
    for i in $prots; do 
	    protname=`echo $i | sed "s:.*/\(.*\):\1:"` 
        mv $spot1ddir/outputs/$protname.spot1d $spotdis2dir/inputs/
    done
fi
#-------------------------------------------------------------------------------
cd $orig_dir
rm -rf $spotdis2dir/spotconlist
rm -rf $spotdis2dir/spot1dlist

echo "All features collected! Running SPOT-Disorder2..."
$spotdis2dir/run_all_models.py --input_list "protlist.txt" --gpu 0
#remove intermediate files- comment out if unnecessary
find $spotdis2dir/outputs -name '*\.[0-9]' | xargs rm
#find $spotdis2dir/outputs -name '*[0-9]' | xargs rm
