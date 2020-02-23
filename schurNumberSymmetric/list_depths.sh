#!/bin/sh

#  list_depths.sh
#  SchurNumber
#
#  Created by rubis on 01/06/2019.
#  Copyright © 2019 rubis. All rights reserved.

p=$1
nmin=$2
nmax=$3
dmin=$4
dmax=$5

function issumfree

{
j="0" #indice de partie
partition=$1
while read -ra set
do
    if [ $j -gt "0" ] && [ $j -le $p ]
    then
        for x in "${set[@]}"
        do
            for y in "${set[@]}"
            do
                if [ $x -le $y ]
                then
                    let "z = x + y"
                    for t in "${set[@]}"
                    do
                        if [ $z -eq $t ]
                        then
                            return "0"
                        fi
                    done
                fi
            done
        done
    fi
    let "j = j + 1"
done <<< "$partition"
return "1"
}

for n in `seq $nmin $nmax`;
do
    let "nmod3 = (n+1) % 3"
    if [ $nmod3 -ne "0" ]
    then
        d=$dmax
        depth="0"
        while [ $d -ge $dmin ] && [ $depth -eq "0" ]
        do
            res=`/Users/rubis/Library/Developer/Xcode/DerivedData/SchurNumber-ewodqxgllntcqmhionbjawrkvroz/Build/Products/Debug/schurNumberSymmetric -p -d $d $p $n`
#echo "$res"
            depth=`expr "$res" : '.*[^0-9]\([0-9]\{1,\}\)'`
#echo "Profondeur $depth"
            let "d = d - 1"
        done
        if [ $depth -ne "0" ]
        then
            issumfree "$res"
            bool=$?
            if [ $bool -eq "0" ]
            then
                echo "Partition avec somme trouvée\n$res"
                break
            fi
        fi
        echo "$n, $depth"
    fi
done
