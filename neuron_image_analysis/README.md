# README #

This repo contains various matlab scripts for manipulating images of neuronal activity.  It is under active development, but most of the scripts should work.

## Installation ##

To use these scripts follow these steps:

* Download and uncompress the source code
* Add the source code to the Matlab path.  You must add the root directory as well as the icons directory.

## Functionality ##

The most useful scripts at this point are:

### nia_playFlatMovie() ###
Displays a movie with tools to interactively adjust colormaps and visualize regions of interest

### nia_movie class ###
A class that stores image data with metadata to provide fast and flexible data manipulation

### nia_playMovie() ###
Displays a nia_movie object with tools similar to those for nia_playFlatMovie()
