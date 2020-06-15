#!/usr/bin/python -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

"""A tiny Python program to evaluate posterior quantiles from csv results
Try running this program from the command line like this:
  python quantiles.py
"""

import sys

def main():
  # Get the name from the command line, using 'World' as a fallback.
  if len(sys.argv) >= 2:
    name = sys.argv[1]
  else:
    print 'Error: filename needs to be supplied as an argument'

  import pandas as pd
  from tabulate import tabulate
  df = pd.read_csv(name)
  print tabulate( df.quantile([0.025, 0.5, 0.975]).transpose().round(3), headers='keys', tablefmt='psql' )

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
