#!/bin/bash
# saving formatted date result to variable
dzis=$(date +%y%m%d)
ls /usr/bin -al > log.$dzis
