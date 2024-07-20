#!/bin/bash
# using formatted result of date command to make name of file
dzis=`date +%d%m%y`
napis=`echo today is $dzis`
echo $napis
ls /usr/bin -al > log.$dzis
