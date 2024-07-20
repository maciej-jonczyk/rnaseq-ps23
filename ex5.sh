#!/bin/bash
# saving computation done in bc program
# syntax
# "semicolon separated opts; expression" and redirection to bc
var1=`echo "scale=3; zmn1=4; zmn1 / 3" | bc`
var2=$(echo "scale=2; zmn1=1.1; zmn1 * 1.1" | bc)
echo The final result is $var1 not $var2
