#!/bin/bash
# using $[ expression ] syntax to make a computation - still only fot integers
var1=100
var2=50
var3=45
var4=$[$var1 * ($var2 - $var3)]
echo The final result is $var4
