'''
First implementation of Ising model
Start date: 24 January 2017
Author: DKarandikar
'''

import numpy
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation


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


def mcmove(lattice, beta, move_n, exponentials):
    '''Flip a spin if the energy change is beneficial'''

    for _ in range(move_n):
        pos_x = numpy.random.randint(0, lattice.shape[0])
        pos_y = numpy.random.randint(0, lattice.shape[1])

        new_e = delta_e(lattice, (pos_x, pos_y))

        if new_e <= 0:
            lattice[pos_x, pos_y] *= -1
        elif exponentials[new_e] > numpy.random.rand():
            lattice[pos_x, pos_y] *= -1

    return lattice


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
    magnet = numpy.sum(lattice)
    return magnet


def mcrun(size, temp, n_0, n_max, move_n, grid_size):
    '''Runs an entire simulation for a given temperature and returns the energy average'''
    lattice = mcsetup(size)

    exponentials = {4 : numpy.exp(-1*4*(1/temp)), 8 : numpy.exp(-1*8*(1/temp))}

    for _ in range(n_0):
        lattice = mcmove(lattice, 1/temp, move_n, exponentials)

    energy = 0
    magnet = 0

    for _ in range(n_max):
        lattice = mcmove(lattice, 1/temp, move_n, exponentials)
        energy += calc_energy(lattice, grid_size)
        magnet += calc_magnet(lattice)

    return (energy/n_max, magnet/n_max)

def ising_graphs(n_0, n_max, move_n, temp_steps, temp_range, testing=False):
    '''Runs multiple temp simulations and then produces an energy-temperature graph'''

    grid_size = 16
    temperatures = numpy.linspace(temp_range[0], temp_range[1], temp_steps)

    energies = numpy.zeros(temp_steps)
    magnetizations = numpy.zeros(temp_steps)

    for k, temp in enumerate(temperatures):
        values = mcrun(grid_size, temp, n_0, n_max, move_n, grid_size)
        energies[k] = values[0]/(grid_size**2)
        magnetizations[k] = abs(values[1]/(grid_size**2))
        print(k/temp_steps)

    figure = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')

    subplot = figure.add_subplot(2, 1, 1)
    plt.plot(temperatures, energies, 'o', color="#A60628", label=' Energy')
    plt.xlabel("Temperature ", fontsize=20)
    plt.ylabel("Energy ", fontsize=20)

    subplot = figure.add_subplot(2, 1, 2)
    plt.plot(temperatures, magnetizations, 'o', color="#A60628", label=' Magnetization')
    plt.xlabel("Temperature ", fontsize=20)
    plt.ylabel("Magnetization ", fontsize=20)

    if not testing:
        plt.show()

    print("done")

def test(testing=False):
    '''Runs a quick test to see that e_t is working'''
    ising_graphs(100, 100, 50, 50, (1, 4), testing)

def time_test():
    '''Performs multiple tests and averages the time'''
    times = 0
    number_tests = 10
    for _ in range(number_tests):
        start = time.time()
        test(True)
        times += (time.time() - start)
    print('Took ' + '%.2f' % (times/number_tests) + ' seconds')


def standard():
    '''Runs a standard density simulation'''
    ising_graphs(400, 400, 300, 100, (1, 4))




def gridplot(temp):
    ''' Displays an animation of the 2d ising model at temperate temp, in units kb=1'''

    def update(i, lattice, beta, exponentials):
        lattice = mcmove(lattice, beta, 150, exponentials)
        mat.set_data(lattice)
        return lattice

    grid_size = 100

    exponentials = {4 : numpy.exp(-1*4*(1/temp)), 8 : numpy.exp(-1*8*(1/temp))}
    beta = 1/temp

    lattice = mcsetup(grid_size)

    fig, axis = plt.subplots()
    mat = axis.matshow(lattice)
    plt.colorbar(mat)
    ani = animation.FuncAnimation(fig, update, fargs=(lattice, beta, exponentials),
                                  interval=50, save_count=50)

    plt.show()
