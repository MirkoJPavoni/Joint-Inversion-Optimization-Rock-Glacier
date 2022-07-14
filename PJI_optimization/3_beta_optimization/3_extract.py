# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 08:55:41 2022

@author: Mirko
"""

import pandas as pd
import glob
file_list = glob.glob('*.txt')
print("Number of files in the list", len(file_list))
print(file_list)

with open('parameters.txt', 'w') as outfile:
    for fname in file_list:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

