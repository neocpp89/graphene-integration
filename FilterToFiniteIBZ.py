#!/usr/bin/env python
import sys
import math

if len(sys.argv) < 2:
    print sys.argv[0]+': GRIDFILE'
    exit(0)

L = math.sqrt(3.0);

sample_length = 5e-6; # meters
a = 2.46e-10; # meters

with open(sys.argv[1]) as f:
    points = map(lambda r: map(float, r.split(' ')), f)
    for p in points:
        x = p[0]
        y = p[1]
        magcp = math.sqrt(x*x + y*y);
        if (magcp == 0):
            theta = 0;
            rho = 0;
        else:
            theta = math.acos(math.cos(6.0*(math.acos(x / magcp) + math.pi / 2.0))) / 6.0;
            rho = L * magcp * math.cos(theta);

        cart_theta = math.atan2(y, x);
        if (cart_theta >= 0 and cart_theta <= math.pi/6.0 and rho > (a / sample_length)):
            print x, y
