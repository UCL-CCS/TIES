#!/usr/bin/env python3
"""
Visualise the bootstrapped replicas. We have 20 replicas altogether.
Show as a functin of a number of bootstrapped replicas how the values change
"""

import os
from pathlib import Path
from collections import OrderedDict
import glob
import time
import itertools
import random
import sys
import json
import re

import numpy as np
from numpy import genfromtxt
import scipy.stats as st
import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

print("hi")
