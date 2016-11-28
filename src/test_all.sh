	#!/bin/bash

counter=0
passed=0
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
