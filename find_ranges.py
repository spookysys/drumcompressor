#!env python3
import fileinput
import re
from glob import glob
import numpy as np

def grab_something(name):
	fnames = glob(r'export/struct/*.inc')
	vals = fileinput.input(fnames)
	vals = [
		re.sub(r'[^\-0-9.,]', r'', x)
		for x in vals
		if re.match(r'//.*// ' + name, x)
	]
	vals = [
		[
			float(y) 
			for y in x.split(",")
		]
		for x in vals
	]
	vals = np.array(vals)

	maxs = np.max(vals, 0)
	mins = np.min(vals, 0)

	print('max: ', maxs)
	print('min: ', mins)

print('Treble filter coefficients')
grab_something('treble filter floats')

print('Treble envelope coefficients')
grab_something('treble envelope floats')
