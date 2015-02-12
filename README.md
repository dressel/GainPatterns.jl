# GainPatterns

The purpose of this package is to allow the plotting of gain patterns. It is assumed that you have the gain pattern in some text file, and wish to plot it. Some basic gain pattern manipulation (like cross-correlation with similar patterns) will be allowed. It would be cool if gain patterns could be generated from simple antennas, but that is not really what this is for.

## Available Functions
`normalize(gains::Vector{Float64})` returns a new vector with normalized gains.

`normalize!(gains::Vector{Float64})` normalizes the provided vector of gains. Similar to `normalize`, but more memory conscious.

`crosscorrelate(ref, sample, angles)` performs a cross-correlation of the sample and a reference signal. The idea is that the reference is of length 360, and the sample can be of variable length. `angles` is a vector of the angles at which `sample` was sampled.

`plotgains(angles::Vector{Float64}, gains::Vector{Float64})` saves a plot of the given gains to a pdf file named temp.pdf in the current directory.

## Near-Term Plans
This is currently still very rough. Some things I want to add:
* Better/more documentation
* Plotting ability
* Maybe create a GainPattern type?
* Allow creation of said type through a text or csv file

## Future Plans
Future plans (like way down the road):
* Allow the use of PyPlot as well as PGFPlots
* Read in an image of a gain pattern and spit out the resulting values (a bit quixotic, but I suppose it is doable).
* Generate gain patterns for simple antennas. I do not currently have the expertise to make that happen.

## Build Status
This came default with the package when I made it. Not sure what to do with it or what it does.
[![Build Status](https://travis-ci.org/dressel/GainPatterns.jl.svg?branch=master)](https://travis-ci.org/dressel/GainPatterns.jl)
