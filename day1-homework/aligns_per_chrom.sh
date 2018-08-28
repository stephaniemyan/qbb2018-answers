# calculates alignments per chromosome when given a .sam file

# this is an efficient-ish method with databash
#!/bin/bash

FILE=$1

grep -v "^@" $FILE | grep -v 2110000 | datamash -s -g 3 count 3

# a less efficient method uses sort, something like:
# sort -k 3 $FILE | uniq -f 2 -c | head
