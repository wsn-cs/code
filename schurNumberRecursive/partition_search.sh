#!/bin/sh

#  partition_search.sh
#  SchurNumber
#
#  Created by rubis on 16/02/2020.
#  Copyright © 2020 rubis. All rights reserved.

cmd1="/Users/rubis/Library/Developer/Xcode/DerivedData/SchurNumber-ewodqxgllntcqmhionbjawrkvroz/Build/Products/Release/schurNumberRecursive"
cmd2="/Users/rubis/Library/Developer/Xcode/DerivedData/SchurNumber-ewodqxgllntcqmhionbjawrkvroz/Build/Products/Release/SchurNumber"

constraint_set1="1 2 4 8 11 22 25 50 53 66 68 71 77"
constraint_set2="3 5 6 7 19 21 23 51 52 63 64 65 74 75 76"
constraint_set3="9 10 12 13 14 15 16 17 18 20 54 55 56 57 58 59 60 61 62"
constraint_set4="24 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49"
constraint_set5="67 69 70 72 73 78"

begin_set1="1 4"
begin_set2="2 3"
begin_set3="5 6 7 8 9 10 11 12 13"

#exec $cmd1 -a -p b 4 $constraint_set1 $constraint_set2 $constraint_set3 $constraint_set4 $begin_set1 $begin_set2 $begin_set3
exec $cmd2 -m w -a -p b 5 $constraint_set1 $constraint_set2 $constraint_set3 $constraint_set4 $constraint_set5
