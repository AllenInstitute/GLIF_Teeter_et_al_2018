'''Check to see if the contents of two csv files are the same.  Note that the files 
must be comma delineated, not tab delineated
'''

import csv
import os
relative_path=os.path.dirname(os.getcwd())
import dictdiffer
import csv

def read_csv(file):
    reader = csv.DictReader(open(file))

    result = {}
    for row in reader:
        key = row.pop('specimen_id')
        if key in result:
            # implement your duplicate row handling here
            result[key+'_2'] = row
        result[key] = row
    return result

file1=os.path.join(relative_path, 'SMfig_4','spikecut_standard_err.csv')
file2=os.path.join(relative_path, 'SMfig_4','out.csv')

d1=read_csv(file1)
d2=read_csv(file2)

the_same=list(dictdiffer.diff(d1,d2))
if the_same==[]:
    print 'your csv files are the same'
else:
    print the_same
    raise Exception('the csv files are not the same')