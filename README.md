# GainPatterns

This package allows the plotting and manipulation of gain patterns. 
A gain pattern is a collection of gains (or strengths) versus angle.
If you have a vector of angles, and a vector of gains taken at these angles, you can create a GainPattern.
You can then create publication-ready plots of these gain patterns using PGF.

I made this package because I was finding it tedious to create these plots.
PGF allows you to make polar plots, but you can't have negative values.
This means you have to shift all gains so the lowest one is zero, and then shift the axis labels back.
This package takes care of all that.

This package also allows you to manipulate gain patterns as well.
For example, you can normalize a gain pattern (subtract mean gain from every element, and divide by std deviation).
You can shift a gain pattern by some number of degrees.
You can add constant values to gain patterns.
You can also sample possible gains from a gain pattern.

This package does not generate gain patterns.
That is, given some antenna configuration, this will not calculate what the gain pattern should look like.
That would be useful, but is outside the scope of this package (for now).
At the moment, the package is mostly a plotting tool.

## Requirements
You need the Julia package PGFPlots to generate the plots.
PGFPlots requires you to have LaTeX set up on your computer.

## Creating GainPatterns
A GainPattern has three fields:
* `angles`, a vector of angles (in degrees) at which the gains are sampled at,
* `meangains`, a vector of the gains. `meangains[i]` is the gain measured at `angles[i]`,
* `samples`, a vector of the samples.

## Plotting
A brief overview is shown here, but check out the [examples](http://nbviewer.ipython.org/github/dressel/GainPatterns.jl/blob/master/doc/GainPatterns.ipynb).
For some reason, the notebook viewer has been acting oddly in Firefox.
If the axis labels don't show up in the examples, try looking at the exmples in another browser.

Plots can be created with:
```
plot(gp)			# plot single GainPattern
plot([gp1, gp2])	# plot multiple GainPatterns on same axis
```

`plot` creates a PolarAxis object (from PGFPlots package).
Currently, it doesn't show anything to the screen (does show in IJulia notebook though).
Use PGFPlots's `save` function to save the file as a PDF or tex file.

```
p = plot(gp)
save("plot.pdf", p)		# saves plot as pdf
save("plot.tex", p)		# saves plot as tex file
```

Optional arguments give you greater control over your plots:
* `ymin`
* `ymax`
* `lastleg` 
* `degrees` Set this to true if you want the angles to have degrees.
* `legendentries` Vector of strings. Length must match length of vector of GainPatterns to plot.
In Julia, optional arguments require you to include the argument name.
Order does not matter.
```
plot(gp, ymin=-100)
plot([gp1,gp2], legendentires=["plot1", "plot2"], degrees=true)
```



## Available Manipulations
On top of plotting, you can perform the following manipulations on a GainPattern.

`normalize!(gp::GainPattern)` modifies gp by normalizing its meangains and samples.

`rotate!(gp::GainPattern, degrees::Real)` adds degrees to gp's angles vector.

`sample(gp::GainPattern, angle::Real)` randomly samples from gp's samples vector at a given angle.
Currently, you must select an angle that exists in gp's angles vector.
In the future, we could do some sort of interpolation.

`sample(gp::GainPattern, angles::Vector{Real})`

`crosscorrelate(ref, sample, angles)` performs a cross-correlation of the sample and a reference signal. The idea is that the reference is of length 360, and the sample can be of variable length. `angles` is a vector of the angles at which `sample` was sampled.


## Near-Term Plans
This is currently still very rough. Some things I want to add:
* Overhaul documentation and create Julia notebook with examples
* ~~Create subtraction~~ done
* Create shift! and flip! functions (for angles...)
* Allow log-axis for gains
* Make plots show up immediately (this is a PGFPlots issue)
* ~~Test `validgain` function and overriding it~~ done
* ~~Create sample function to sample from a gain pattern at a given angle~~ done
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
