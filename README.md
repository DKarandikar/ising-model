# ising-model

Basic Ising model using monte-carlo in python

# Overview and Aims

Basic implementation of a monte-carlo method to the 2d Ising model. 

Phase change can be seen around the T=2.5 mark.

Aims: want to implement more graphs, such as heat capacity to see the phase shift from more angles. 

## Installation

To run on windows, if git is installed:

        git clone https://github.com/DKarandikar/ising-model.git
        cd ising-model
        pip install .

Otherwise use the github clone function directly and then just pip in the download location.


## Testing

From command line in directory run:

        nosetests
or 

        python setup.py test
        
Otherwise, install it and then in python type:

        import ising
        ising.test()

Result should be a matplotlib graph of the form:

<img src="images/testingexample.png" width="500">

## Usage

To use more generally import ising, then use:

        ising.standard()

or: 

        ising.gridplot(1)

or more generally:

        ising.ising_graphs(n_0, n_max, move_n, T_steps, T_range)

Where:

- n_0 is the number of iterations to run to equilibrate 
- n_max is the number of iterations to run per temperate step on top of n_0
- move_n is the number of times to attempt to flip a spin per step
- T_steps is the number of discrete temperate points to run for 
- T_range is a tuple (a,b) where a is the initial temperature and b is the final

Overall, the script should try to flip a spin (n_0 + n_max)* move_n * T_steps times per graph

## Changelog

### [0.4.2] - 29/01/2017

- Finished migrating to cython code
- Anything that doesn't require matplotlib should now be in the pyd file
- Significant speed increase as a result

### [0.4.1] - 28/01/2017

- Significant structural changes
- Setup.py files changed to actually work

### [0.4.0] - 28/01/2017

- Implemented some cython that is temporarily working 
- Mcmove method is approximately 10x faster

### [0.3.0] - 27/01/2017

- Gridplot now generates two plots and animates them both
- The second plot is always +5 temperature on the passed value
- Changed ising_graphs to display magnetization neater

### [0.2.0] - 27/01/2017

- Added method to generate animation, called gridplot

### [0.1.0] - 25/01/2017

- Migrated to a package format installed with pip or python install
- Minimal changes to main functionality
- Migrated test to a proper nosetests structure
- Setup should have relevant dependencies 

### [0.0.1] - 25/01/2017

- Initial version using a .py file in main directory
- Basic implementation established for only E-T graph