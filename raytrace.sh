#!/bin/bash

for i in {1..50}
do
    time1=1400
    time2=$time1+10 
    file=Assets/project.lua

    sed -i '' "s/$time1/$time2/" $file

    ./a4 Assets/project.lua 
done