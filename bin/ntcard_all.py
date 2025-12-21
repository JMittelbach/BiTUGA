#!/usr/bin/env python3

import sys, os, argparse

def parse_command_line(argv):
  p = argparse.ArgumentParser(description=('apply ntcard_mn.x to populus '
                                          'and ginkgo'))
  p.add_argument('kmer_length',type=int,
                  help='specify length of k-mers')
  return p.parse_args(argv)

args = parse_command_line(sys.argv[1:])

for species, suffix in [('populus_tremula', 'fastq'), ('ginkgo_biloba', 'fq')]:
  envvar = species.upper()
  directory = os.environ[envvar]
  if not os.path.isdir(directory):
    sys.stderr.write(f'{sys.argv[0]}: directory "{directory}" does not exist\n')
    exit(1)
  print(f'touch {species}.txt')
  for filename in os.listdir(directory):
    print(f'ntcard_mn.x --qgram_length {args.kmer_length} --binary '
          f'{directory}/{filename} >> {species}.txt')
