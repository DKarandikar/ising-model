import numpy


def mcmove(lattice, move_n, exp_low, exp_high):
    '''Flip a spin if the energy change is beneficial'''

    for _ in range(move_n):
        pos_x = numpy.random.randint(0, lattice.shape[0])
        pos_y = numpy.random.randint(0, lattice.shape[1])

        new_e = delta_e(lattice, (pos_x, pos_y))

        if new_e <= 0:
            lattice[pos_x, pos_y] *= -1
        elif new_e == 4 and exp_low > numpy.random.rand():
            lattice[pos_x, pos_y] *= -1
        elif new_e == 8 and exp_high > numpy.random.rand():
            lattice[pos_x, pos_y] *= -1

    return lattice

def delta_e(lattice, site): 
    '''Calculates the new delta E for a given site'''
    pos_x = site[0]
    pos_y = site[1]
    size_x = lattice.shape[0]
    size_y = lattice.shape[1]
    energy = 2 * lattice[pos_x, pos_y] * (lattice[(pos_x-1)%size_x, pos_y] \
                                   + lattice[(pos_x+1)%size_x, pos_y] \
                                   + lattice[pos_x, (pos_y-1)%size_y] \
                                   + lattice[pos_x, (pos_y+1)%size_y])
    return energy