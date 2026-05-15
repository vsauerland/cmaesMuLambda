#!/bin/bash

numSessions=$(head -4 nIter.txt | tail -1 | awk '{print $1}')
maxIter=$(head -5 nIter.txt | tail -1 | awk '{print $1}')
startIter=$(tail -1 nIter.txt | awk '{print $1}')

for nIter in $(seq "$startIter" "$maxIter")
do
        gen=$((nIter+1))
        genDir="gen$gen"
        mkdir -p "$genDir"

        ./cmaes "$nIter" nIter.txt 1 > "$genDir/serial_${nIter}.log" || exit

        for nSession in $(seq 1 "$numSessions")
        do
                memberDir="$genDir/member$nSession"
                mkdir -p "$memberDir"

                paramFile="$genDir/parameters_${gen}_${nSession}.txt"

                while ! [ -e "$paramFile" ]
                do
                        sleep 1
                done

                cp "$paramFile" "$memberDir/parameters.txt"

                (
                        cd "$memberDir" || exit

                        ../../testfunctions \
                                ../../nIter.txt \
                                parameters.txt \
                                fitness.txt \
                                > "parallel_${gen}_${nSession}.log"
                ) &
        done

        for nSession in $(seq 1 "$numSessions")
        do
                memberDir="$genDir/member$nSession"
                fitnessFile="$genDir/fitness_${gen}_${nSession}.txt"

                while ! [ -s "$memberDir/fitness.txt" ]
                do
                        sleep 1
                done

                cp "$memberDir/fitness.txt" "$fitnessFile"
        done
done
