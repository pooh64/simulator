#!/usr/bin/env python3

import sys
import re

if len(sys.argv) != 2:
	sys.exit("wrong args")

sim_runtime = '''
_start:
movliu r2, 32767
bl r1, main
exit
__flush:
flush
br r1
__bkpt:
bkpt
br r1
'''
print(sim_runtime)

lines = []
with open(sys.argv[1]) as file:
    lines = file.readlines()
    lines = [line.lstrip().rstrip() for line in lines]

re_comment = re.compile(r"\s*(;.*)?$")
re_label = re.compile(r"^(.?\w*):$")

for l in lines:
	l = re.sub(re_comment, '', l)
	if len(l) == 0:
		continue
	if l[0] == '.':
		m = re.search(re_label, l)
		if not m:
			continue
	l = l.replace("\t", "\t ")
	print(l)
