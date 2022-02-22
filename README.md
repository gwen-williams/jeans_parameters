Calculate the density assuming spherical symmetry, the thermal velocity dispersion, the Jeans length, and the Jeans mass for any given mass, radius and temperature.

## Dependencies

Python 3.9

* argparse
* numpy
* astropy

## Usage

In iPython:

`run jeans_parameters.py -h` 

Use the -h argument for help on input arguments.

`run jeans_parameters.py`

If no arguments are specified by the user, the script will run with default values of 100 Msun, 0.1 pc and 10 K.

`run jeans_parameters.py -m 200 -r 0.2 -t 20`

Will run with user defined arguments of 200 Msun, 0.2 pc, and 20 K.


