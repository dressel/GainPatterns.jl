# GainPatterns

The purpose of this package is to allow the plotting of gain patterns. It is assumed that you have the gain pattern in some text file, and wish to plot it. Some basic gain pattern manipulation (like cross-correlation with similar patterns) will be allowed. It would be cool if gain patterns could be generated from simple antennas, but that is not really what this is for.

## Available Functions
`normalize(gains::Vector{Float64})` returns a new vector with normalized gains.

`normalize!(gains::Vector{Float64})` normalizes the provided vector of gains. Similar to `normalize`, but more memory conscious.

`crosscorrelate(ref, sample, angles)` performs a cross-correlation of the sample and a reference signal. The idea is that the reference is of length 360, and the sample can be of variable length. `angles` is a vector of the angles at which `sample` was sampled.

`plotgains(angles::Vector{Real}, gains::Vector{Real})` creates a PolarAxis object containing the plot of the gains vs the angles. You can save it using `save` from PGFPlots.jl.

`plotgains(angles::Vector{Real}, gains::Vector{Real}, ymin::Real)` is the same as above, except it allows you to specify a minimum y (radial) value.

`plotgains(angles::Vector{Real}, gains::Matrix{Real}, ymin::Real)` plots errors along the radial direction. Each row in gains corresponds to a specific angle. The values of the columns correspond to different samples taken at that angle.

`plotgains(angles::Vector{Real}, gains::Matrix{Real}, ymin::Real)` is the same as above, except it allows you to specify a minimum y (radial) value.

## Example Usage
```
angles = [0:360]
gains = float(cosd(angles))   # returns array of type any sometimes??
gp = GainPattern(angles, gains)
axis = plot(gp)
save("plot1.pdf", axis)
```
Now say you want to make the minimum y (radial) value -2, and not -1. You can specify this:
```
angles = [0:360]
gains = float(cosd(angles))   # returns array of type any sometimes??
gp = GainPattern(angles, gains)
axis = plot(gp, ymin=-2)
save("plot2.pdf", axis)
```

## Near-Term Plans
This is currently still very rough. Some things I want to add:
* ~~Create a GainPattern type~~ done!
* ~~Allow creation of said type through a text or csv file~~ done!
* Overhaul documentation and create Julia notebook with examples
* ~~make ymin work for plot~~ done!
* Create subtraction
* Create shift! and flip! functions (for angles...)
* Allow log-axis for gains
* Make plots show up immediately (this is a PGFPlots issue)
* Test `validgain` function and over-writing it

## Future Plans
Future plans (like way down the road):
* Allow the use of PyPlot as well as PGFPlots
* Read in an image of a gain pattern and spit out the resulting values (a bit quixotic, but I suppose it is doable).
* Generate gain patterns for simple antennas. I do not currently have the knowledge to make that happen.

## Build Status
This came default with the package when I made it. Not sure what to do with it or what it does.
[![Build Status](https://travis-ci.org/dressel/GainPatterns.jl.svg?branch=master)](https://travis-ci.org/dressel/GainPatterns.jl)
