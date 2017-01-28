
import numpy
import timeit
import move
#import move_cy

def mcsetup(grid_size):
    '''Setups up the initial lattice and spins'''
    lattice = 2 * numpy.random.randint(2, size=(grid_size, grid_size)) - 1
    return lattice

exponentials = {4 : numpy.exp(-1*4*(1/2)), 8 : numpy.exp(-1*8*(1/2))}


lattice = mcsetup(100)

p= timeit.timeit('move.mcmove(lattice,200,exponential_low, exponential_high)',"import move; import numpy; lattice = 2 * numpy.random.randint(2, size=(100,100)) - 1 ; exponential_low = numpy.exp(-1*4*(1/2));  exponential_high = numpy.exp(-1*8*(1/2))", number=1000)

print(p)

p= timeit.timeit('move_cy.mcmove(lattice,200,exponential_low, exponential_high)',"import move_cy; import numpy; lattice = 2 * numpy.random.randint(2, size=(100,100)) - 1 ; exponential_low = numpy.exp(-1*4*(1/2));  exponential_high = numpy.exp(-1*8*(1/2))", number=1000)

print(p)

#print(move_cy.calc_energy(mcsetup(16),16))