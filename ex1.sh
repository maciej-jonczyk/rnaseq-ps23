#!/bin/bash
# two syntaxes of doing command substitution
test=`date`
echo today is $test
test2=$(date)
echo another today is $test2
