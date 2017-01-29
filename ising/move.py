import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation


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
        lattice = mcmove(lattice, move_n, exponential_low, exponential_high)

    energy = 0
    magnet = 0

    for _ in range(n_max):
        lattice = mcmove(lattice, move_n, exponential_low, exponential_high)
        energy += calc_energy(lattice, grid_size)
        magnet += calc_magnet(lattice)

    return (energy/n_max, magnet/n_max)

    
def ising_graphs(n_0, n_max, move_n, temp_steps, temp_range, testing=False):
    '''Runs multiple temp simulations and then produces relevant graphs'''

    grid_size = 16

    temperatures = numpy.linspace(temp_range[0], temp_range[1], temp_steps)
    energies = numpy.zeros(temp_steps)
    magnetizations = numpy.zeros(temp_steps)

    for k, temp in enumerate(temperatures):
        values = mcrun(temp, n_0, n_max, move_n, grid_size)
        energies[k] = values[0]/(grid_size**2)
        magnetizations[k] = abs(values[1]/(grid_size**2))
        print(k/temp_steps)

    figure = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')

    dummy_sp = figure.add_subplot(1, 2, 1)
    plt.plot(temperatures, energies, 'o', color="#A60628", label=' Energy')
    plt.xlabel(r'$Temperature [J/k_{B}]$', fontsize=20)
    plt.ylabel(r'$Energy [J]$', fontsize=20)

    dummy_sp = figure.add_subplot(1, 2, 2)
    plt.plot(temperatures, magnetizations, 'o', color="#A60628", label=' Magnetization')
    plt.xlabel(r'$Temperature [J/k_{B}]$', fontsize=20)
    plt.ylabel(r'$Magnetization [\mu]$', fontsize=20)

    if not testing:
        plt.show()

    print("done")

ising_graphs(100, 100, 50, 50, (1, 4), False)