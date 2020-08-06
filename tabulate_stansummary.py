#!/usr/bin/python -tt
# the structure of this script comes from Google's Python Class
# http://code.google.com/edu/languages/google-python-class/
#
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

"""A Python program that reads in the results of stansummary and prints them out with nice formatting
  python tabulate_stansummary.py
"""

import sys

def main():
  if len(sys.argv) >= 2:
    name = sys.argv[1]
  else:
    print 'Error: filename needs to be supplied as an argument'

  import pandas as pd
  from tabulate import tabulate
  df = pd.read_csv(name, comment='#', index_col=0)
  df = df[['5%', '50%', '95%', 'N_Eff', 'N_Eff/s']]
  for name in df.index:
    if name in ['accept_stat__', 'treedepth__', 'stepsize__', 'n_leapfrog__', 'energy__', 'divergent__',
                'theta[1]', 'theta[2]', 'theta[3]', 'theta[4]', 'theta[5]']:
      df = df.drop([name], axis=0)
  df = df.dropna()
  df['N_Eff'] = df['N_Eff'].astype(int)
  print tabulate( df.round(3), headers='keys', tablefmt='psql' )

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
