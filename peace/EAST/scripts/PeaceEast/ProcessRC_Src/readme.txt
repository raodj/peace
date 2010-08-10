ProcessRC.cpp:
make all the reads into their reverse complemented form.

Input: estfile, mstfile, newestfile
Output: an ace info file with the name "$estfile.aceinfo". This file has the format "indexOfNode 1/-1". (1:forward direction, -1: reverse.) Each index has one line.

