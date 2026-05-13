numSessions=$(head -4 nIter.txt | tail -1 | awk '{print $1}')
maxIter=$(head -5 nIter.txt | tail -1 | awk '{print $1}')
startIter=$(tail -1 nIter.txt | awk '{print $1}')
for nIter in $(seq $startIter $maxIter)
do

        ./cmaes $nIter nIter.txt gen$nIter > serial_$nIter.log
        if [ "$?" -ne 0 ]; then
                exit
        fi
        for nSession in $(seq 1 $numSessions)
        do
                ## wait until required parameter file for model simulation exists:
                while ! [ -e "parameters_"$((nIter+1))"_"$nSession".txt" ]
                do
                        sleep 1
                done
                ./testfunctions nIter.txt "gen$((nIter+1))" "parameters_"$((nIter+1))"_"$nSession".txt" "fitness_"$((nIter+1))"_"$nSession".txt" > "parallel_"$((nIter+1))"_"$nSession".log" &
        done
        ## wait until required fitness files for next optimization step exist:
        for nSession in $(seq 1 $numSessions)
        do
                while ! [ -e "fitness_"$((nIter+1))"_"$nSession".txt" ]
                do
                        sleep 1
                done
        done
done
