import numpy
cimport numpy
import cython


from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport exp 

# cython: profile=True

@cython.cdivision(True) 


def mcsetup(int grid_size):
    '''Setups up the initial lattice and spins'''
    cdef int i
    cdef int j
    cdef numpy.ndarray[int, ndim=2, mode='c'] lattice = numpy.zeros((grid_size, grid_size), dtype=int)
    cdef int random 
    for i in range(grid_size):
        for j in range(grid_size):
            random = rand()%2
            if random==0:
                lattice[i,j] = 1
            else:
                lattice[i,j] = -1
    return lattice


def mcmove(numpy.ndarray[int, ndim=2, mode='c'] lattice, int move_n, float exp_low, float exp_high):
    '''Flip a spin if the energy change is beneficial'''
    cdef int pos_x
    cdef int pos_y
    cdef int new_e
    cdef int _
    cdef numpy.ndarray[int, ndim=2, mode='c'] output = lattice
    cdef int size_x = lattice.shape[0]
    cdef int size_y = lattice.shape[1]
    cdef double rand_float

    for _ in range(move_n):
        #pos_x = numpy.random.randint(0, lattice.shape[0])
        pos_x = int(rand() % size_x)
        #pos_y = numpy.random.randint(0, lattice.shape[1])
        pos_y = int(rand() % size_y)

        new_e = delta_e(output, pos_x, pos_y)

        rand_float = rand()
        rand_float /= RAND_MAX

        if new_e <= 0:
            output[pos_x, pos_y] *= -1
        elif new_e == 4 and exp_low > rand_float:
            output[pos_x, pos_y] *= -1
        elif new_e == 8 and exp_high > rand_float:
            output[pos_x, pos_y] *= -1

    return output

cdef int delta_e(numpy.ndarray[int, ndim=2, mode='c'] lattice, int pos_x, int pos_y):
    '''Calculates the new delta E for a given site'''
    cdef int size_x = lattice.shape[0]
    cdef int size_y = lattice.shape[1]
    cdef int energy = 2 * lattice[pos_x, pos_y] * (lattice[(pos_x-1)%size_x, pos_y] \
                                   + lattice[(pos_x+1)%size_x, pos_y] \
                                   + lattice[pos_x, (pos_y-1)%size_y] \
                                   + lattice[pos_x, (pos_y+1)%size_y])
    return energy

cdef int calc_magnet(numpy.ndarray[int, ndim=2, mode='c'] lattice):
    '''Sums over all spins to calculate magnetization'''
    cdef int magnet = numpy.sum(lattice)
    return magnet


cdef float calc_energy(numpy.ndarray[int, ndim=2, mode='c'] lattice, int grid_size):
    '''Calculates the average energy of the system'''
    cdef int result = 0
    cdef int i, j, site_e
    cdef int size_x = lattice.shape[0]
    cdef int size_y = lattice.shape[1]
    
    for i in range(size_x):
        for j in range(size_y):
            site_e = lattice[i, j] * -1  * (lattice[(i-1)%grid_size, j] \
                                   + lattice[(i+1)%grid_size, j] \
                                   + lattice[i, (j-1)%grid_size] \
                                   + lattice[i, (j+1)%grid_size])
            result += site_e
    cdef float p = result/4
    return p


def mcrun(float temp, int n_0, int n_max, int move_n, int grid_size):
    '''Runs an entire simulation for a given temperature and returns the energy average'''
    cdef numpy.ndarray[int, ndim=2, mode='c'] lattice = mcsetup(grid_size)

    exponential_low = numpy.exp(-1*4*(1/temp))
    exponential_high = numpy.exp(-1*8*(1/temp))

    for _ in range(n_0):
        lattice = mcmove(lattice, move_n, exponential_low, exponential_high)

    energy = 0
    magnet = 0

    for _ in range(n_max):
        lattice = mcmove(lattice, move_n, exponential_low, exponential_high)
        energy += calc_energy(lattice, grid_size)
        magnet += calc_magnet(lattice)

    return (energy/n_max, magnet/n_max)
