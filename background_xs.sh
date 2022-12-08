n=10000
run(){
    p=$1
    e=$2
    reaction=$3
    mkdir /data/sebouh/hepmc
    mkdir /data/sebouh/raw
    mkdir /data/sebouh/tuples
    hepmcfile=/data/sebouh/hepmc/${reaction}_raw_${1}x${2}_tmp.hepmc
    rawfile=/data/sebouh/raw/${reaction}_${1}x${2}.root
    outfile=/data/sebouh/tuples/${reaction}_${1}x${2}_tuples.root
    #rm $rawfile $hepmcfile $outfile
    mkdir tmp
    cmndfile=tmp/${reaction}_${1}x${2}.cmnd
    cat pythia8cards/${reaction}.cmnd | sed 's/Beams:eA *=.*/Beams:eA  = '$1'  ! proton energy/g' | sed 's/Beams:eB *= *.*/Beams:eB  = '$2' ! electron energy/g' | sed 's/Main:numberOfEvents *=.*/Main:numberOfEvents = '$n'/g' > $cmndfile

    main42 $cmndfile $hepmcfile
    #DelphesHepMC3 tcl/ATHENA_smear2.tcl $rawfile $hepmcfile
    #DelphesHepMC3 tcl/delphes_EICmatrixv2_3T.tcl $rawfile $hepmcfile
    #python3 python/eicroot2pandas.py $rawfile $outfile --n=${n} --match=old --reaction=$reaction &
    #outfile=/data/sebouh/tuples/${reaction}_${1}x${2}_tuples_highestEnergy.root
    #python3 python/eicroot2pandas.py $rawfile $outfile --n=${n} --match=highestEnergy --reaction=$reaction &
}

#for reaction in Photoproduction NC_DIS CC_DIS
for reaction in NC_DIS NC_DIS_positron Photoproduction Photoproduction_positron CC_DIS CC_DIS_positron
do
    
    #run 41 5 $reaction &
    #run 100 5 $reaction &
    #run 100 10 $reaction &
    echo "<"$reaction">"
    run 275 10 $reaction
    echo "</"$reaction">"
    #run 275 18 $reaction &
    
#    wait
#    hadd -f ${reaction}_all.root ${reaction}_tuples_*x*.root 
done
#wait
