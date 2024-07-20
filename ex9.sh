#!/bin/bash
# if-else and using variable
testuser=adam
if grep $testuser /etc/passwd; then
	echo "It's ok"
	echo "Another ok"
	ls *8*
else
	echo "Not ok, the $testuser not exist"
fi

