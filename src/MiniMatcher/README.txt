code ist nun verf"ugbar in src/MiniMatcher

cd src/MiniMatcher

CXX=clang++ make

cd ${GTTL}/testsuite && make split_files_mn.x

Annahme: die Programme nt_mini_matcher.x sowie split_files_mn.x

sind in Deinem PATH.

Annahme: diese Datei sind im aktuellen Verzeichnis:

$ ls SRR2744682*.fastq.gz  unitigs.fa 
SRR2744682_1.fastq.gz  SRR2744682_2.fastq.gz  unitigs.fa

$ nt_mini_matcher_wrapper.py  SRR2744682 | sh -s > SRR2744682.tsv

match_file_cleaner.py SRR2744682.tsv

liefert dann die Ausgabe mit korrigierten Readnummern.
