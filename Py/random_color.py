#!/usr/bin/env python3

from argparse import ArgumentParser
from randomcolor import RandomColor

ap = ArgumentParser(description='Generate some colors')
ap.add_argument('-c','--count', type = int, 
        default = 1, help='number of colors to create')
args = vars(ap.parse_args())

rand_color = RandomColor()

params = {
        'hue':None,
        'luminosity':None,
        'count':args.get('count'),
        'format_':'hex'
        }

random_colors = rand_color.generate(**params)

print(random_colors)
