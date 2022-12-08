n=1000000
run(){
    p=$1
    e=$2
    reaction=$3
    mkdir hepmc
    mkdir raw
    mkdir tuples
    hepmcfile=hepmc/${reaction}_raw_${1}x${2}.hepmc
    rawfile=tmp/${reaction}_${1}x${2}.root
    outfile=tuples/${reaction}_${1}x${2}_tuples.root
    rm $rawfile $hepmcfile $outfile
    mkdir tmp
    cmndfile=tmp/${reaction}_${1}x${2}.cmnd
    cat pythia8cards/${reaction}.cmnd | sed 's/Beams:eA *=.*/Beams:eA  = '$1'  ! proton energy/g' | sed 's/Beams:eB *= *.*/Beams:eB  = '$2' ! electron energy/g' | sed 's/Main:numberOfEvents *=.*/Main:numberOfEvents = '$n'/g' > $cmndfile

    main42 $cmndfile $hepmcfile
    DelphesHepMC3 tcl/ATHENA_smear2.tcl $rawfile $hepmcfile 
    python3 python/eicroot2pandas.py $rawfile $outfile --n=1000000 --match=old
}

#for reaction in Photoproduction NC_DIS CC_DIS
for reaction in NC_DIS
do
    
    run 41 5 $reaction &
    run 100 5 $reaction &
    run 100 10 $reaction &
    run 275 10 $reaction &
    run 275 18 $reaction &
    
#    wait
#    hadd -f ${reaction}_all.root ${reaction}_tuples_*x*.root 
done
wait
