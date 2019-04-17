#!/usr/bin/env python

import sys
import numpy as np

things=[]
thingSet = set([])
counts = None
thingIndexDict = {}

for line in sys.stdin:
    lineThings = line.split()
    if counts is None:
        ncol = len(lineThings)
        colInds = range(ncol)
        counts = np.zeros(shape=(0,ncol))
    for i in colInds:
        if lineThings[i] in thingSet: counts[thingIndexDict[lineThings[i]],i] += 1
        else:
            thing = lineThings[i]
            things.append(thing)
            thingIndexDict[thing] = len(things) - 1
            thingSet.add(thing)
            counts = np.vstack([counts, np.zeros(ncol)])
            counts[-1,i] += 1

counts = counts.astype(int)

output = np.hstack([np.array([[t] for t in things]), counts.astype(str)])

output = output[output[:,0].argsort()]

sys.stdout.write("\n".join(["\t".join(list(row)) for row in output])+ "\n")

