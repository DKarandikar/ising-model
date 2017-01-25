'''
First implementation of Ising model
Start date: 24 January 2017
Author: Daniel Karandikar
'''

import numpy
import matplotlib.pyplot as plt


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


def mcmove(lattice, beta, move_n):
    '''Flip a spin'''
    for _ in range(move_n):
        pos_x = numpy.random.randint(0, lattice.shape[0])
        pos_y = numpy.random.randint(0, lattice.shape[1])

        new_e = delta_e(lattice, (pos_x, pos_y))

        if new_e <= 0:
            lattice[pos_x, pos_y] *= -1
        elif numpy.exp(-1 * new_e * beta) > numpy.random.rand():
            lattice[pos_x, pos_y] *= -1

    return lattice


def mcsetup(grid_size):
    '''Setups up the initial lattice and spins'''
    lattice = 2 * numpy.random.randint(2, size=(grid_size, grid_size)) - 1
    return lattice

def calc_energy(lattice):
    '''Calculates the average energy of the system'''
    result = 0
    size_x = lattice.shape[0]
    size_y = lattice.shape[1]
    for i in range(lattice.shape[0]):
        for j in range(lattice.shape[1]):
            site_e = lattice[i, j] * -1  * (lattice[(i-1)%size_x, j] \
                                   + lattice[(i+1)%size_x, j] \
                                   + lattice[i, (j-1)%size_y] \
                                   + lattice[i, (j+1)%size_y])
            result += site_e

    return result/(4*size_x*size_y)


def mcrun(size, temp, n_0, n_max, move_n):
    '''Runs an entire simulation for a given temperature and returns the energy average'''
    lattice = mcsetup(size)
    for _ in range(n_0):
        lattice = mcmove(lattice, 1/temp, move_n)

    energy = 0

    for _ in range(n_max):
        lattice = mcmove(lattice, 1/temp, move_n)
        energy += calc_energy(lattice)

    return energy/n_max

def e_t_graph(n_0, n_max, move_n, temp_steps):
    '''Runs multiple temp simulations and then produces an energy-temperature graph'''

    grid_size = 16
    temperatures = numpy.linspace(1, 4, temp_steps)
    energies = numpy.zeros(temp_steps)

    for k in range(len(temperatures)):
        energies[k] = mcrun(grid_size, temperatures[k], n_0, n_max, move_n)
        print(k/100)

    #figure = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')

    plt.plot(temperatures, energies, 'o', color="#A60628", label=' Energy')
    plt.xlabel("Temperature (T)", fontsize=20)
    plt.ylabel("Energy ", fontsize=20)
    plt.show()

    print("done")

def configplot():
    '''Plots multiple subplots'''

    n = 16
    f= plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')

    lattice = mcsetup(n)
    for i in range(5001):
        lattice = mcmove(lattice, 1/3)
        if i == 1:
            subplot(f, lattice, i, n, 2)
        if i == 10:
            subplot(f, lattice, i, n, 3)
        if i == 100:
            subplot(f, lattice, i, n, 4)
        if i == 100:
            subplot(f, lattice, i, n, 5)
        if i == 5000:
            subplot(f, lattice, i, n, 6)
        #print(i)

    plt.show()

def subplot(f, lattice, i, n, n_):
    '''Appends a subplot to f'''
    X, Y = numpy.meshgrid(range(n), range(n))
    sp = f.add_subplot(3, 3, n_ )
    plt.pcolormesh(X, Y, lattice)



### Main Section ###

N_0 = 400
N_MAX = 400
MOVE_N = 300
T_STEPS = 100

testing = "y"

if testing == "n":
    e_t_graph(N_0, N_MAX, MOVE_N, T_STEPS)
elif testing == "y":
    e_t_graph(100, 100, 50, 50)


#configplot()
