#!/usr/bin/env python
import sys
import math

if len(sys.argv) < 2:
    print sys.argv[0]+': GRIDFILE'
    exit(0)

with open(sys.argv[1]) as f:
    points = map(lambda r: map(float, r.split(' ')), f)
    for p in points:
        x = p[0]
        y = p[1]
        theta = math.atan2(y,x)
        if (theta >= 0 and theta <= math.pi/6.0):
            print x, y
