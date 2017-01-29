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

def mcsetup(grid_size):
    '''Setups up the initial lattice and spins'''
    lattice = 2 * numpy.random.randint(2, size=(grid_size, grid_size)) - 1
    return lattice


def calc_energy(lattice, grid_size):
    '''Calculates the average energy of the system'''
    result = 0
    for i in range(lattice.shape[0]):
        for j in range(lattice.shape[1]):
            site_e = lattice[i, j] * -1  * (lattice[(i-1)%grid_size, j] \
                                   + lattice[(i+1)%grid_size, j] \
                                   + lattice[i, (j-1)%grid_size] \
                                   + lattice[i, (j+1)%grid_size])
            result += site_e

    return result/(4)

def calc_magnet(lattice):
    '''Sums over all spins to calculate magnetization'''
    magnet = numpy.sum(lattice)
    return magnet


def mcrun(temp, n_0, n_max, move_n, grid_size):
    '''Runs an entire simulation for a given temperature and returns the energy average'''
    lattice = mcsetup(grid_size)

    exponential_low = numpy.exp(-1*4*(1/temp))
    exponential_high = numpy.exp(-1*8*(1/temp))

    for _ in range(n_0):
        lattice = move_cy.mcmove(lattice, move_n, exponential_low, exponential_high)

    energy = 0
    magnet = 0

    for _ in range(n_max):
        lattice = move_cy.mcmove(lattice, move_n, exponential_low, exponential_high)
        energy += move_cy.calc_energy(lattice, grid_size)
        magnet += calc_magnet(lattice)

    return (energy/n_max, magnet/n_max)