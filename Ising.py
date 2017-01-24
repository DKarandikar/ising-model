'''
First implementation of Ising model
Start date: 24 January 2017
Author: Daniel Karandikar
'''

import numpy
import matplotlib.pyplot as plt

N_0 = 400
N_MAX = 400
T_STEPS = 50

def delta_e(lattice, site):
    '''Calculates the new delta E for a given site'''
    x = site[0]
    y = site[1]
    size_x = lattice.shape[0]
    size_y = lattice.shape[1]
    energy = 2 * lattice[x, y] * (lattice[(x-1)%size_x, y] \
                                   + lattice[(x+1)%size_x, y] \
                                   + lattice[x, (y-1)%size_y] \
                                   + lattice[x, (y+1)%size_y])
    return energy


def mcmove(lattice, beta):
    '''Flip a spin'''
    size_x = lattice.shape[0]
    size_y = lattice.shape[1]
    for i in range(size_x):
        for j in range(size_y):
            x = numpy.random.randint(0, lattice.shape[0])
            y = numpy.random.randint(0, lattice.shape[1])

            new_e = delta_e(lattice, (x, y))

            if new_e <= 0:
                lattice[x, y] *= -1
            elif numpy.exp(-1 * new_e * beta) > numpy.random.rand():
                lattice[x, y] *= -1

    return lattice


def mcsetup(N):
    '''Setups up the initial lattice and spins'''
    lattice = 2 * numpy.random.randint(2, size=(N, N)) - 1
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


def mcrun(size, temp):
    lattice = mcsetup(size)
    for i in range(N_0):
        lattice = mcmove(lattice, 1/temp)

    energy = 0

    for i in range(N_MAX):
        lattice = mcmove(lattice, 1/temp)
        energy += calc_energy(lattice)

    return energy/N_MAX

def e_t_graph():
    temperatures = numpy.linspace(1, 4, T_STEPS)
    energies = numpy.zeros(T_STEPS)

    for k in range(len(temperatures)):
        energies[k] = mcrun(16, temperatures[k])
        print(k/100)

    figure = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')

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

e_t_graph()

#configplot()
