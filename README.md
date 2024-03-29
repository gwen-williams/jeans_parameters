Calculate the density assuming spherical symmetry, the thermal velocity dispersion, the Jeans length, and the Jeans mass for any given mass, radius and temperature.

Currently, the script assumes $\mu$ = 2.8, and a core formation efficiency (CFE) of 13%. These can be changed on lines 174 and 175.

## Dependencies

Python 3.9

* argparse
* numpy
* astropy

## Usage

In iPython:

``` python

run jeans_parameters.py -h

``` 

Use the -h argument for help on input arguments.

``` python

run jeans_parameters.py

```

If no arguments are specified by the user, the script will run with default values of 100 Msun, 0.1 pc and 10 K.


``` python

run jeans_parameters.py -m 200 -r 0.2 -t 20

```

Will run with user defined arguments of 200 Msun, 0.2 pc, and 20 K.


