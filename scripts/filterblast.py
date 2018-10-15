#!/usr/bin/env python

import sys
import os

inp = sys.argv[1]
oup = sys.argv[2]
scripts = sys.argv[3]

# need this to get local modules
sys.path.append(scripts)
from helpers import helpers as hp

# filter top hit
hp.tophitsfilter(inp, oup)
