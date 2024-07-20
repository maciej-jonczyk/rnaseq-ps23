#!/bin/bash
# using both variable from bc program and from outside (var1)
var1=10
var2=$(echo "scale=2; zmn1=1.1; zmn1 * 1.1 / $var1" | bc)
echo The final result is $var1 not $var2
