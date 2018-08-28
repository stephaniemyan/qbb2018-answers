# calculates alignments per chromosome when given a .sam file

#!/bin/bash

FILE=$1

grep -v "^@" $FILE | grep -v 2110000 | datamash -s -g 3 count 3
