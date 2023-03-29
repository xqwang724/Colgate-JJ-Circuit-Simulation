import matplotlib.pyplot as plt
import math
import numpy as np
import os
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.tri as tri
from os import path
import time
import csv
import pandas as pd

def readOutputs(inputFile):
	dataset = []
	with open(inputFile) as inFile:
		rows = inFile.readlines()
		num = len(rows)
		for rawrow in rows[:]:
			row = rawrow.split()
			dataset.append(row)
		dataset = np.array(dataset).astype(float)
	return dataset, num

def processOutputs(inputFile,duration):
	data, N = readOutputs(inputFile)
	# print(y)
	#duration = 1e-10
	times = np.linspace(0,duration,N)

	y = []
	for i in range(len(data[0])):
		y.append([])
		for j in range(N):
			y[i].append(data[j][i])
	return times,y

if __name__ == "__main__":
	data, N = readOutputs("output_singlejunction_IV_0.02p.txt")
	# print(y)
	duration = 0.1e-9
	times = np.linspace(0,duration,N)

	y = []
	for i in range(len(data[0])):
		y.append([])
		for j in range(N):
			y[i].append(data[j][i])

	for i in range(len(data[0])):
		plt.plot(times, y[i])
		plt.show()








