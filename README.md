# GainPatterns

This package allows manipulation of gain patterns.
A gain pattern is a collection of gains (or strengths) versus angle.
If you have a vector of angles, and a vector of gains taken at these angles, you can create a GainPattern.
You can then create publication-ready plots of these gain patterns using PGF with the GainPatternsTex.jl package.

This package also allows you to manipulate gain patterns as well.
For example, you can normalize a gain pattern (subtract mean gain from every element, and divide by std deviation).
You can shift a gain pattern by some number of degrees.
You can add constant values to gain patterns.
You can also sample possible gains from a gain pattern.

This package does not generate gain patterns.
That is, given some antenna configuration, this will not calculate what the gain pattern should look like.
That would be useful, but is outside the scope of this package (for now).
At the moment, the package is mostly a plotting tool.

If you would like to plot GainPatterns, check out the GainPatternsTex.jl package.

## Creating GainPatterns
A GainPattern has three fields:
* `angles`, a vector of angles (in degrees) at which the gains are sampled at,
* `meangains`, a vector of the gains. `meangains[i]` is the gain measured at `angles[i]`,
* `samples`, a vector of the samples.

You can create a gain pattern with the samples.
You can also create a gain pattern from a csv file, where the first column is the angles at which measurements were taken (in degrees, 0-359), and any subsequent columns are gain values taken at those degrees.
```
gp1 = GainPattern(angles, gains)
gp2 = GainPattern("file.csv")
```

A GainPattern can be saved to a csv file using `csv`.
First column is angles, second column is meangains.
```
csv(gp)						# saves file to temp.csv
csv(gp, "filename.csv")		# saves file to filename.csv
```

## Available Manipulations
On top of plotting, you can perform the following manipulations on a GainPattern.

`normalize!(gp::GainPattern)` modifies gp by normalizing its meangains and samples.

`rotate!(gp::GainPattern, degrees::Real)` adds degrees to gp's angles vector.

`sample(gp::GainPattern, angle::Real)` randomly samples from gp's samples vector at a given angle.
Currently, you must select an angle that exists in gp's angles vector.
In the future, we could do some sort of interpolation.

`sample(gp::GainPattern, angles::Vector{Real})` allows you to sample from a specified set of angles.
Each sample in the angles parameter must exist in gp.

`sample(gp::GainPattern)` samples from gp's set of angles.
A gain is sampled for every angle of gp.
This is equivalent to calling `sample(gp, gp.angles)`.

`crosscorrelate(ref, sample, angles)` performs a cross-correlation of the sample and a reference signal. The idea is that the reference is of length 360, and the sample can be of variable length. `angles` is a vector of the angles at which `sample` was sampled.

## Adding and Subtracting GainPatterns
You can also add and subtract gain patterns.
Addition and subtraction return new gain patterns.
For example, `gp1 + gp2` will create a new GainPattern.
`gp1` and `gp2` need to have the same number of angles for addition/subtraction.

If the length of the sample vectors are the same for both `gp1` and `gp2`, then the sample vectors are added/subtracted.
If `gp1` and `gp2` have sample vectors of different lengths, then only the mean gains are added/subtracted,
and the new gain pattern will be created just from this vector of added/subtracted mean gains.

You can also add a constant to a gain pattern.
This adds a constant to the meangain and every sample.

If you simply want to combine two gain patterns, use `addsamples!`, as shown below:
```
addsamples!(gp1, gp2)
```
This will add the append the saples of gp2 to the samples of gp1.

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
