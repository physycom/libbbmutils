#!/bin/bash

counter=0
passed=0
# Check if test log ".old" copies exist and create them if needed
old=$(find test -type f -name "*.log.old")
if [[ $old == "" ]]; then
	echo "First test call, creating old copies of output files."
	for test in test/*.log; do
		cp $test $test.old
	done
fi
# Test loop compare new test log with old ones
for test in test/*.log; do
	((counter++))
	df=$(diff $test $test.old)
	testname=$(basename $test)
	testname="${testname%.*}"
	echo -n "Test for $testname : "
	if [[ $df == "" ]]; then 
		echo "OK"
		((passed++))
	else
		echo "failed"
	fi
done
echo "Test passed      : $passed/$counter ( $(($passed/$counter*100)) % )"	
