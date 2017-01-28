'''
First implementation of Ising model
Start date: 24 January 2017
Author: DKarandikar
'''
import time

import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import move_cy


def mcsetup(grid_size):
    '''Setups up the initial lattice and spins'''
    lattice = 2 * numpy.random.randint(2, size=(grid_size, grid_size)) - 1
    return lattice



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


def gridplot(temp):
    ''' Displays an animation of the 2d ising model at temperate temp, in units kb=1'''

    grid_size = 100

    exponential_low = numpy.exp(-1*4*(1/temp))
    exponential_high = numpy.exp(-1*8*(1/temp))
    exponential_low2 = numpy.exp(-1*4*(1/(temp+5)))
    exponential_high2 = numpy.exp(-1*8*(1/(temp+5)))

    lattice = mcsetup(grid_size)
    lattice2 = mcsetup(grid_size)

    def data_gen():
        '''Returns the lattice to the animation update'''
        while True:
            yield lattice, lattice2

    def update(data):
        ''' Updates the animation by flipping (potentially) 150 spins every tick '''
        lattice, lattice2 = data

        lattice = mcmove(lattice, 150, exponential_low, exponential_high)
        matrix1.set_data(lattice)

        lattice2 = mcmove(lattice2, 150, exponential_low2, exponential_high2)
        matrix2.set_data(lattice2)

        return None

    fig = plt.figure(figsize=(18, 10), dpi=80)
    axis1 = plt.subplot(121)
    axis2 = plt.subplot(122)

    axis1.axis('off')
    axis1.set_title("Ising model at T=" + str(temp))
    axis2.axis('off')
    axis2.set_title("Ising model at T=" + str(temp+5))

    matrix1 = axis1.matshow(lattice, cmap='spring')
    matrix2 = axis2.matshow(lattice2, cmap='spring')

    dummy_ani = animation.FuncAnimation(fig, update, data_gen,
                                        interval=50, save_count=50)

    plt.show()



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
