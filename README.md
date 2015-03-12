# GainPatterns

The purpose of this package is to allow the plotting of gain patterns. It is assumed that you have the gain pattern in some text file, and wish to plot it. Some basic gain pattern manipulation (like cross-correlation with similar patterns) will be allowed. It would be cool if gain patterns could be generated from simple antennas, but that is not really what this is for.

## Available Functions
`normalize(gains::Vector{Float64})` returns a new vector with normalized gains.

`normalize!(gains::Vector{Float64})` normalizes the provided vector of gains. Similar to `normalize`, but more memory conscious.

`crosscorrelate(ref, sample, angles)` performs a cross-correlation of the sample and a reference signal. The idea is that the reference is of length 360, and the sample can be of variable length. `angles` is a vector of the angles at which `sample` was sampled.


## Example Usage
Check out the [examples](http://nbviewer.ipython.org/github/dressel/GainPatterns.jl/blob/master/doc/GainPatterns.ipynb).

```
angles = [0:359]
gains = float(cosd(angles))   # returns array of type any sometimes??
gp = GainPattern(angles, gains)
p = plot(gp)
save("plot1.pdf", p)
```
Now say you want to make the minimum y (radial) value -2, and not -1. You can specify this:
```
angles = [0:359]
gains = float(cosd(angles))   # returns array of type any sometimes??
gp = GainPattern(angles, gains)
p = plot(gp, ymin=-2)
save("plot2.pdf", p)
```
You can specify the maximum radial value as well, using the `ymax` parameter, similarly to how the `ymin` parameter is used.

Notice how the plot is always closed? We only have values for angles 0 to 359, but there is a line between 359 and zero. That is because GainPatterns automatically draws a line between the last point and the first point, to make the gain pattern look continuous over all angles.

Sometimes, this behavior is undesirable. For example, maybe you get more than one revolution of data; maybe you start measuring gains at 45 degrees and loop all the way around to 47 degrees. It might be silly to draw a line connecting the last and first points in such a case. For these situations (or any time you'd like), you can set the optional argument `lastleg` to false when plotting.
```
p = plot(gp, lastleg=false)
```
You can pass in a vector of gain patterns to plot them on a single axis.
```
p = plot([gp1, gp2])
```
If you are using `plot` to plot multiple gain patterns, you can specify the legend entries.
```
p = plot([gp1, gp2], legendentries=["Gain Pattern 1", "Gain Pattern 2"])
```
If the dimensions of the pattern array and the legend entries array don't match, it will probably fail. I must check this.

You can plot stuff and show the degrees with:
```
plot(gp, degrees=true)
```


## Near-Term Plans
This is currently still very rough. Some things I want to add:
* ~~Create a GainPattern type~~ done
* ~~Allow creation of said type through a text or csv file~~ done
* Overhaul documentation and create Julia notebook with examples
* ~~make ymin work for plot~~ done
* ~~Create subtraction~~ done
* Create shift! and flip! functions (for angles...)
* Allow log-axis for gains
* Make plots show up immediately (this is a PGFPlots issue)
* ~~Test `validgain` function and over-writing it~~ done
* Create sample function to sample from a gain pattern at a given angle
* If a samples vector has only one entry, eliminate it? or maybe not. Plotting doesn't seem to have an issue with this
* Group plots (you can't use tex grouplot package for this)
* ~~Sample an entire distribution and return an array (or fill one)~~done
* ~~Normalize function for gainpatterns~~done
* Allow for units on axes
* function that also prints samples  (not just mean gains)
* Check if file exists with open (csv +elsewhere)
* Handle file extension stuff better (csv +elsewhere)

## Future Plans
Future plans (like way down the road):
* Allow the use of PyPlot as well as PGFPlots
* Read in an image of a gain pattern and spit out the resulting values (a bit quixotic, but I suppose it is doable).
* Generate gain patterns for simple antennas. I do not currently have the knowledge to make that happen.

## Build Status
This came default with the package when I made it. Not sure what to do with it or what it does.
[![Build Status](https://travis-ci.org/dressel/GainPatterns.jl.svg?branch=master)](https://travis-ci.org/dressel/GainPatterns.jl)
