#/bin/bash

#number of processes

np=90

i=2

SOURCEDIR=$HOME/makedx_v20

#echo deploying makeda01
#cd $SOURCEDIR/makedx01/src
#make clean
#make makeda_he
#make makedx_v20
#$SOURCEDIR/makedx01
#mv makeda_he makedx01 
#cp makedx_v20 makedx01
#./clean.sh

while [ $i -le 9 ]; do
    echo deploying makedx0$i
#    mkdir makedx0$i
#    cd makedx0$i
#    mkdir src
#    mkdir masses
    cd $SOURCEDIR/makedx01/
#    cp controlparams $SOURCEDIR/makedx0$i/
#    cp makedx01.sh $SOURCEDIR/makedx0$i/makedx0$i.sh
    cp saverun.sh $SOURCEDIR/makedx0$i/
#    rm $SOURCEDIR/makedx0$i/src/*
#    cp src/* $SOURCEDIR/makedx0$i/src/
#    cp masses/* $SOURCEDIR/makedx0$i/masses/
#    cp clean.sh AUXIN5 EEOS* IEOS* inputprof SQOPAC $SOURCEDIR/makedx0$i/
#    cp clean.sh saverun.sh $SOURCEDIR/makedx0$i/
#    cd $SOURCEDIR/makedx0$i/src
#    make clean
#    make makeda_he
#    make makedx_v20
#    cd $SOURCEDIR/makedx0$i
#    sed -i -e "s/x01/x0$i/g" makedx0$i.sh
#    chmod +x makedx0$i.sh
#    mv makedx_v20 makedx0$i
#    mv makeda_he makedx0$i 
#    cp makedx_mesa makedx0$i
#    ./clean.sh
    let i=i+1
done

while [ $i -le $np ]; do
    echo deploying makedx$i
#    mkdir makedx0$i
#    cd makedx$i
#    mkdir src
#    mkdir masses
    cd $SOURCEDIR/makedx01/
#    cp controlparams $SOURCEDIR/makedx$i/
#    cp makedx01.sh $SOURCEDIR/makedx$i/makedx$i.sh
    cp saverun.sh $SOURCEDIR/makedx$i/
#    rm $SOURCEDIR/makedx$i/src/*
#    cp src/* $SOURCEDIR/makedx$i/src/
#    cp masses/* $SOURCEDIR/makedx0$i/masses/
#    clean.sh AUXIN5 EEOS* IEOS* inputprof SQOPAC $SOURCEDIR/makedx0$i/
#    cp clean.sh saverun.sh $SOURCEDIR/makedx$i/
#    cd $SOURCEDIR/makedx$i/src
#    make clean
#    make makeda_he
#    make makedx_v20
#    cd $SOURCEDIR/makedx$i
#    sed -i -e "s/x01/x$i/g" makedx$i.sh
#    chmod +x makedx$i.sh
#    mv makedx_v20 makedx$i
#    mv makeda_he makedx$i 
#   cp makedx_mesa makedx$i
#   ./clean.sh
    let i=i+1
done

#cd /home/axk55/src/
#pwd

