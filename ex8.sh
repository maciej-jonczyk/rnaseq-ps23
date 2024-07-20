#!/bin/bash
# using variable value (must be integer) custom exit status number
# explicit value like "exit 5" can be used
# value must be <=255, otherwise modulo (remainder of division by 255 is used)
var1=10
var2=43
var3=33
var4=71

var5=$(bc << EOF
scale=4
a1 = ($var1 * $var2)
b1 = ($var3 * $var4)
b1 - a1
EOF
)

echo The final result is $var1 not $var5
exit $var5
