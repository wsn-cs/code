#!/bin/sh

#  partition_search.sh
#  SchurNumber
#
#  Created by rubis on 16/02/2020.
#  Copyright © 2020 rubis. All rights reserved.

cmd1="/Users/rubis/Library/Developer/Xcode/DerivedData/SchurNumber-ewodqxgllntcqmhionbjawrkvroz/Build/Products/Release/schurNumberRecursive"
cmd2="/Users/rubis/Library/Developer/Xcode/DerivedData/SchurNumber-ewodqxgllntcqmhionbjawrkvroz/Build/Products/Release/SchurNumber"

constraint_set1='"'"1 2 4 8 11 22 25 50 63 69 135 140 150 155 178 183 193"'"'
constraint_set2='"'"3 5 6 7 19 21 23 51 52 53 64 65 66 137 138 139 151 152 153 180 181 182 194 195 196"'"'
constraint_set3='"'"9 10 12 13 14 15 16 17 18 20 54 55 56 57 58 59 60 61 62 141 142 143 144 145 146 147 148 149 184 185 186 187 188 189 190 191 192"'"'
constraint_set4='"'"24 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 154 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 179"'"'
constraint_set5='"'"67 68 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 136"'"'

begin_set1='"'"4 14 19 42 47 57"'"'
begin_set2='"'"1 2 3 15 16 17 44 45 46 58 59 60"'"'
begin_set3='"'"5 6 7 8 9 10 11 12 13 48 49 50 51 52 53 54 55 56"'"'
begin_set4='"'"18 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43"'"'

exec $cmd1 -a -p b 5 $constraint_set1 $constraint_set2 $constraint_set3 $constraint_set4 $constraint_set5 $begin_set1 $begin_set2 $begin_set3 $begin_set4
#exec $cmd2 -m w -a -p b 5 $constraint_set1 $constraint_set2 $constraint_set3 $constraint_set4 $constraint_set5
