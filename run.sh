#!/bin/bash

# You need to set these three variables
IMGPATH=/Users/yuncong/Documents/workspace/StairsPointCloud/staircase_new
OUTPUTPATH=/Users/yuncong/Documents/workspace/StairsPointCloud/staircase_new_result
PROGSPATH=/Users/yuncong/Documents/testStereoAlgs

IMGID=$1
ALG=$2
LEFTNAME=top${IMGID}_rect_rot_small
RIGHTNAME=bottom${IMGID}_rect_rot_small

LEFTIMG=$IMGPATH/$LEFTNAME.pgm
RIGHTIMG=$IMGPATH/$RIGHTNAME.pgm

#convert $IMGPATH/top${IMGID}_rect_rot.pgm -resize 25% $IMGPATH/top${IMGID}_rect_rot_small.pgm
#convert $IMGPATH/top${IMGID}_rect_rot.pgm -resize 25% $IMGPATH/top${IMGID}_rect_rot_small.pgm

if [ $# -ne 3 ]; then
    echo Usage: source run.sh algorithm image_id params_file
    return
fi

if [ "$ALG" == "elas" ]; then
    ELAS_BIN=$PROGSPATH/libelas/elas
    ELAS_PARAMSFILE=$3
    $ELAS_BIN $LEFTIMG $RIGHTIMG $ELAS_PARAMSFILE
    mv $IMGPATH/${LEFTNAME}_disp_elas.pgm $OUTPUTPATH/${LEFTNAME}_disp_elas.pgm
    mv $IMGPATH/${RIGHTNAME}_disp_elas.pgm $OUTPUTPATH/${RIGHTNAME}_disp_elas.pgm
elif [ "$ALG" == "mrf" ]; then
    MRF_BIN=$PROGSPATH/MRFStereo/mrfstereo/mrfstereo
    MRF_PARAMSFILE=$3
    MRF_PARAMS=`awk 'NR == 20' $MRF_PARAMSFILE`
    $MRF_BIN $MRF_PARAMS $LEFTIMG $RIGHTIMG $OUTPUTPATH/${LEFTNAME}_disp_mrf.pgm
elif [ "$ALG" == "sm" ]; then
    SM_BIN=$PROGSPATH/StereoMatch/StereoMatch
    SCRIPTPATH=$3
    OLDPATH=`pwd`
    cd $SCRIPTPATH
    mkdir stairs/results
    echo "input_file $LEFTIMG
input_file $RIGHTIMG" > stairs/data_in.txt
    $SM_BIN script exp_all.txt
    if [ -d $OUTPUTPATH/${LEFTNAME}_disp_sm ]; then
        mv stairs/results/* $OUTPUTPATH/${LEFTNAME}_disp_sm
        rm -r stairs/results
    else
        mv stairs/results $OUTPUTPATH/${LEFTNAME}_disp_sm
    fi
    cd $OLDPATH
else
    echo algorithm $ALG not recognized
    return
fi
