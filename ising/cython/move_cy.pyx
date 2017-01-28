import numpy
cimport numpy
import cython


from libc.stdlib cimport rand, RAND_MAX

# cython: profile=True

@cython.cdivision(True) 

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