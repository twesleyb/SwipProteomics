#!/usr/bin/env python3
''' A thin wrapper around randomcolor.RandomColor.generate() '''

## IMPORTS
from argparse import ArgumentParser
from randomcolor import RandomColor

# Define a function that drops None type arguments.
def rmNone(args):
    rm = [key for key in args if args.get(key) is None]
    for key in rm: del args[key]
    return(args)


## INPUTS
ap = ArgumentParser(description='Generate a random color.')

# Count
ap.add_argument('-c','--count', default = 1, type = int, help='''
An integer which specifies the number of colors to generate.''')

# Hue
ap.add_argument('-u','--hue', type = str, default=None, help='''
Controls the hue of the generated color. You can pass a string
representing a color name: red, orange, yellow, green, blue, purple, pink and
monochrome are currently supported.''')

# Luminosity
ap.add_argument('-l','--luminosity', type = str, default=None, help='''
Controls the luminosity of the generated color. You can specify a
string containing bright, light, or dark.''')

# Seed
# FIXME: doesnt work
#ap.add_argument('-s','--seed', type = int, default=None, help='''
#An integer which when passed will cause randomColor to return
#the same color each time.''')

# Format
ap.add_argument('-f','--format', type = str, default=None, help='''
A string which specifies the format of the generated color. Possible
values are rgb, rgbArray, hsl, hslArray and hex (default).''')


# Parse input arguments.
args = vars(ap.parse_args())
args = rmNone(args)

# Generate n random colors.
rand_color = RandomColor()
print(rand_color.generate(**args))
