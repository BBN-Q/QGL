import os, sys, glob
from parse import *

qasm_path = "/home/gribeill/GitHub/openqasm/examples/*.qasm"

fails = []

files = glob.glob(qasm_path)

print(files)

for fn in files:
	p = QASM3Parser()
	with open(fn) as f: 
		try:
			p.build_tree(f.read())
		except:
			fails.append(fn)

Nt = len(files)
Nf = len(fails)

print(f"Failed on {Nf} out of {Nt} files!")
for f in fails:
	print(f)
