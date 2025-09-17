import crpropa as crp
from tqdm import tqdm
from time import time
import numpy as np
from matplotlib import pyplot as plt
import pandas
import matplotlib

matplotlib.rcParams.update({'font.size': 23})
matplotlib.rcParams.update({'lines.linewidth': 2.5})
a = 15

def test(Nm, NTests = 1000, output=False):
	Brms = 1 * crp.muG
	lmin = 0.1 * crp.pc
	lmax = 100 * crp.pc

	Spectrum = crp.TurbulenceSpectrum(Brms, lmin, lmax)
	Turbulence = crp.PlaneWaveTurbulence(Spectrum, int(Nm), 42)
	RNGenerator = crp.Random(137)

	start = time()
	for i in tqdm(range(int(NTests)), disable= not output):
		vec = RNGenerator.randVector()*RNGenerator.rand()*crp.kpc
		Turbulence.getField(vec)
	stop = time()
	duration = stop - start
	if(output):
		print( "Duration of", NTests, "calls of Turbulence.getField:")
		print(" ", duration)
	return duration

def genData():
	resolution = 25
	countForAverage = 5
	NTests = 10000

	NumberOfModesArr = np.logspace(1, 7, resolution) 
	durations = np.array([])
	for Nm in tqdm(NumberOfModesArr):
		durtmp = 0
		for j in range(countForAverage):
			durtmp += test(Nm, NTests)
		durations = np.append(durations, durtmp/countForAverage)

	Data = pandas.DataFrame({"Nm" : NumberOfModesArr, "duration" : durations})
	Data.to_csv("PWTCUDADurations.csv")

def genDurationPlot():
	kinds = ["CPU", "CUDA", "SIMD", "SIMPLEGPU"]
	Files = [f"PWT{kind}Durations.csv" for kind in kinds]

	fig = plt.figure(figsize=(a, a*3.5/6))

	for i, file in enumerate(Files):
		data = pandas.read_csv(file)
		plt.plot(
			data["Nm"].values,
			data["duration"].values,
			"-o",
			label=f"Method = {kinds[i]}"
		)

	plt.xlabel("Number of Modes")
	plt.ylabel("Duration [s]")
	plt.grid()
	plt.loglog()
	plt.legend()

	fig.savefig("PWTDurationsPlot.pdf")

def genComparisonPlot():
	kinds = ["CPU", "CUDA", "SIMD", "SIMPLEGPU"]
	Files = [f"PWT{kind}Durations.csv" for kind in kinds]

	fig = plt.figure(figsize=(a, a*3.5/6))

	CUDA_DATA = pandas.read_csv(Files[1])

	for i, file in enumerate(Files):
		if "CUDA" in file:
			continue

		data = pandas.read_csv(file)
		plt.plot(
			data["Nm"].values,
			data["duration"].values/CUDA_DATA["duration"].values[:len(data["duration"].values)],
			"-o",
			label=f"Method = {kinds[i]}"
		)

	plt.xlabel("Number of Modes")
	plt.ylabel("Duration/CUDADuration")
	plt.grid()
	plt.loglog()
	plt.legend()

	fig.savefig("PWTComparisonPlot.pdf")

if __name__ == "__main__":
	genData()
	genDurationPlot()
	genComparisonPlot()