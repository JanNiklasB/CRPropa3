# This is a helper used by cmake/FindCRPropa.cmake as a fallback locator.
# Prints the path to CRPropa's swig_interface directory or install prefix. 
# Otherwise, it exits with a non-zero status if CRPropa is not importable.

import sys

try:
	import crpropa
except ImportError as e:
	sys.stderr.write("Could not import crpropa!\n")
	sys.exit(-1)

if len(sys.argv) < 2:
	sys.stderr.write("No argument was given, either use 'swig_interface' or 'install_prefix'.\n")
	sys.exit(-1)

if sys.argv[1] == 'swig_interface':
	sys.stdout.write(crpropa.getDataPath('swig_interface'))
elif sys.argv[1] == 'install_prefix':
	sys.stdout.write(crpropa.getInstallPrefix())
else:
	sys.stderr.write("No known argument, either use 'swig_interface' or 'install_prefix'.\n")
	sys.exit(-1)
