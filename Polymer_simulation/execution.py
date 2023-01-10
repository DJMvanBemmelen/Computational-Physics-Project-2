from class_polymer import Polymer
import numpy as np
import matplotlib.pyplot as plt
import random as rd
import math as m
import scipy.optimize as sc
import os
import time


def generation(L_max, N, c_plus, name_dir):
    """
    Generates the polymers
    
    Parameters
    ----------
    L_max : int
        Maximum length of the polymers
    N : int
        Number of iniialy generated polymers
    c_plus : float
        Upperbound value for the PERM algortihm (if c_plus = 0 the PERM method is not used)
    name_dir : string
        Sets name of the directory for the outputted data

    Returns
    -------
    polymer_ls : list
        List of pointers to the generated Polymer objects
    """

    # time simulation start
    print('\n Start build progress:')
    start = time.time()

    # Initialise N polymers 
    polymer_ls = []
    for k in range(N):
        polymer = Polymer()
        polymer_ls.append(polymer)

    # Lists for polymer grow data
    if c_plus != 0:
        prune_ls = []
        enrich_ls = []
        nr_polymers_after_ls = []
        prune_count = 0
        enrich_count = 0
    nr_polymers_before_ls = []

    # Grow the polymers til L_max
    for L in range(L_max):
        nr_polymers_before = 0
        for polymer in polymer_ls:
            polymer.grow()
            if polymer.weight[-1] != 0:
                # Count how many polymers are left
                nr_polymers_before += 1

        nr_polymers_before_ls.append(nr_polymers_before)

        # Implement the PERM method 
        if c_plus != 0:
            c_minus = 0.1*c_plus
            polymer_ls, prune_count, enrich_count, nr_polymers = PERM(polymer_ls, c_plus, c_minus, prune_count, enrich_count)

            prune_ls.append(prune_count)
            enrich_ls.append(enrich_count)
            nr_polymers_after_ls.append(nr_polymers)
        
        # Display progress
        if L % (L_max/4) == 0:
                percentage = L/(L_max)*100
                print("Grow percentage ", percentage)

    # Determining runtime
    finish = time.time()
    runtime = round(finish - start, 2)
    print('\n Finished polymer generation in Runtime:', runtime, 'seconds\n')

    # Save polymer grow data plot
    if c_plus != 0:
        np.save(os.path.join(name_dir, 'prune_ls'), prune_ls)
        np.save(os.path.join(name_dir, 'enrich_ls'), enrich_ls)
        np.save(os.path.join(name_dir, 'nr_polymers_after_ls'), nr_polymers_after_ls)
    
    np.save(os.path.join(name_dir, 'nr_polymers_before_ls'), nr_polymers_before_ls)

    with open(os.path.join(name_dir, 'input.txt'), "w") as file:
        file.write("===== The input data of the simulation ===== \n")
        file.write("Maximum polymer length (L_max):\n")
        file.write(f"{L_max}\n")
        file.write("Number of initialy generated polymers:\n")
        file.write(f"{N}\n")
        file.write("Upperbound value for the PERM algortihm (if c_plus = 0 the PERM method is not used):\n")
        file.write(f"{c_plus}\n")

    if c_plus != 0:
        with open(os.path.join(name_dir, 'input.txt'), "a") as file:
            file.write("Total number of prunes:\n")
            file.write(f"{prune_count}\n")
            file.write("Total number of enrichments:\n")
            file.write(f"{enrich_count}\n")

    return polymer_ls


def generate_free(L_max, N, c_plus, name_dir):
    """
    Generates the polymers as free random walks. 
    
    Parameters
    ----------
    L_max : int
        Maximum length of the polymers
    N : int
        Number of iniialy generated polymers
    c_plus : float
        Upperbound value for the PERM algortihm (if c_plus = 0 the PERM method is not used)
    name_dir : string
        Sets name of the directory for the outputted data

    Returns
    -------
    polymer_ls : list
        List of pointers to the generated Polymer objects
    """
    # time simulation start
    print('\n Start build progress:')
    start = time.time()

    # Lists for polymer grow data
    if c_plus != 0:
        prune_ls = []
        enrich_ls = []
        nr_polymers_after_ls = []
        prune_count = 0
        enrich_count = 0
    nr_polymers_before_ls = []

    # Initialise N polymers 
    polymer_ls = []
    for k in range(N):
        polymer = Polymer()
        polymer_ls.append(polymer)
    
    # Grow polymers to length L_max
    for L in range(L_max):
        nr_polymers_before = 0
        for polymer in polymer_ls:
            polymer.grow_free()
            if polymer.weight[-1] != 0:
                # Count how many polymers are left
                nr_polymers_before += 1

        nr_polymers_before_ls.append(nr_polymers_before)
        # Implement the PERM method 
        if c_plus != 0:
            c_minus = 0.1*c_plus
            polymer_ls, prune_count, enrich_count, nr_polymers = PERM(polymer_ls, c_plus, c_minus, prune_count, enrich_count)

            prune_ls.append(prune_count)
            enrich_ls.append(enrich_count)
            nr_polymers_after_ls.append(nr_polymers)
        # Display progress
        if L % (L_max/4) == 0:
                percentage = L/(L_max)*100
                print("Grow percentage ", percentage)
    
    # Determining runtime
    finish = time.time()
    runtime = round(finish - start, 2)
    print('\n Finished polymer generation in Runtime:', runtime, 'seconds\n')
    
    # Save polymer grow data plot
    if c_plus != 0:
        np.save(os.path.join(name_dir, 'prune_ls'), prune_ls)
        np.save(os.path.join(name_dir, 'enrich_ls'), enrich_ls)
        np.save(os.path.join(name_dir, 'nr_polymers_after_ls'), nr_polymers_after_ls)
    
    np.save(os.path.join(name_dir, 'nr_polymers_before_ls'), nr_polymers_before_ls)

    with open(os.path.join(name_dir, 'input.txt'), "w") as file:
        file.write("===== The input data of the simulation ===== \n")
        file.write("Maximum polymer length (L_max):\n")
        file.write(f"{L_max}\n")
        file.write("Number of initialy generated polymers:\n")
        file.write(f"{N}\n")
        file.write("Upperbound value for the PERM algortihm (if c_plus = 0 the PERM method is not used):\n")
        file.write(f"{c_plus}\n")

    if c_plus != 0:
        with open(os.path.join(name_dir, 'input.txt'), "a") as file:
            file.write("Total number of prunes:\n")
            file.write(f"{prune_count}\n")
            file.write("Total number of enrichments:\n")
            file.write(f"{enrich_count}\n")

    return polymer_ls


def observables(L_max, polymer_ls, n_bootstrap_samples, name_dir):
    """
    Computes the observables and the error of the observables using eigther the bootstrapping or the approximate analytical method

    Parameters
    ----------
    L_max : int
        Maximum length of the polymers
    n_bootstrap_samples : int
        Number of bootstrap samples used in the bootstrapping method (if 0 the aprox analytical method is used)
    polymer_ls : list
        List of pointers to all Polymer objects
    name_dir : string
        Sets name of the directory for the outputted data

    """

    # time simulation start
    print('\n Start observable calculation progress:')
    start = time.time()
    
    end_to_end = np.zeros(L_max-1, dtype="float")
    end_to_end_error = np.zeros(L_max-1, dtype="float")
    r_gyration = np.zeros(L_max-1, dtype="float")
    r_gyration_error = np.zeros(L_max-1, dtype="float")
    
    if n_bootstrap_samples > 0:
        BS_samples = bootstrap_samples(polymer_ls, n_bootstrap_samples)
        end_to_end_BSerror = np.zeros(L_max-1, dtype="float")
        r_gyration_BSerror = np.zeros(L_max-1, dtype="float")

    for L in range(1, L_max):
        end_to_end[L - 1] = average(polymer_ls, 'end-to-end', L)
        end_to_end_error[L - 1] = approx_analytic_error(polymer_ls, L, 'end-to-end', end_to_end[L - 1])
        r_gyration[L - 1] = average(polymer_ls, 'r_gyration', L)
        r_gyration_error[L - 1] = approx_analytic_error(polymer_ls, L, 'r_gyration', r_gyration[L - 1])
        
        if n_bootstrap_samples > 0:
            end_to_end_BSerror[L - 1] = bootstrap_error(BS_samples, n_bootstrap_samples, L, 'end-to-end')
            r_gyration_BSerror[L - 1] = bootstrap_error(BS_samples, n_bootstrap_samples, L, 'r_gyration')
    
        # Display progress
        if L % (L_max/4) == 0:
                percentage = L/(L_max)*100
                print("Observable percentage ", percentage)
    

    # Determining runtime
    finish = time.time()
    runtime = round(finish - start, 2)
    print('\n Finished observable calculation in Runtime:', runtime, 'seconds\n')

    np.save(os.path.join(name_dir, 'end_to_end'), end_to_end)
    np.save(os.path.join(name_dir, 'end_to_end_error'), end_to_end_error)
    np.save(os.path.join(name_dir, 'r_gyration'), r_gyration)
    np.save(os.path.join(name_dir, 'r_gyration_error'), r_gyration_error)

    if n_bootstrap_samples > 0:
        np.save(os.path.join(name_dir, 'end_to_end_BSerror'), end_to_end_BSerror)
        np.save(os.path.join(name_dir, 'r_gyration_BSerror'), r_gyration_BSerror)
        with open(os.path.join(name_dir, 'input.txt'), "a") as file:
            file.write("Number of bootstrap samples used in the bootstrapping method (if 0 the aprox analytical method is used):\n")
            file.write(f"{n_bootstrap_samples}\n")



def average(polymers, observable, L):
    """
    Computes the weighted average of the requested observable over the requested number of datapoints starting at the first lattice point.

    Parameters
    ----------
    polymers : list
        List of pointers to all Polymer objects
    observable : string
        Indicates of which observable the average should be calculated. Allowed values: 'end-to-end', 'r_gyration'
    L : int
        The number of points to be included in the calculation
    
    Return
    ------
    avr : float

    """
    observable_ls = np.zeros(len(polymers), dtype=np.longdouble)
    weights = np.zeros(len(polymers), dtype=np.longdouble)

    # get observable data
    if observable == 'end-to-end':
        for i in range(len(polymers)):
            observable_ls[i] = polymers[i].end_to_end(L)
            weights[i] = polymers[i].weight[L - 1]

    elif observable == 'r_gyration':
        for i in range(len(polymers)):
            observable_ls[i] = polymers[i].r_gyration(L)
            weights[i] = polymers[i].weight[L - 1]
    else:
        print('Misspelled "observable"')
        return 0

    # calculate weighted average
    numerator = np.sum(weights*observable_ls)
    denominator = np.sum(weights)

    # avoid division by 0
    if denominator == 0:
        return 0

    avg = numerator/denominator


    return avg

def approx_analytic_error(polymers, L, observable, avg):
    """
    Computes an approximation of the error on the weighted average of an observable using the formula given in the lecture notes.

    polymers : list
        List containing the original sample of polymer objects
    L : int
        Length at which the observable average was computed
    observable : str
        The observable of which to compute the error. Only accepted values are: 'end-to-end', 'r_gyration'
    avg : numpy.longdouble
        The weighted average of the observable of the original sample
    """
    if observable != 'end-to-end' and observable != 'r_gyration':
        print('Misspelled "observable"')
        return 0

    N = len(polymers)
    sample_size = N
    numerator = 0
    denominator = 0

    # construct quantity from eq (5) on https://compphys.quantumtinkerer.tudelft.nl/proj2-polymers
    for i in range(N):
        polymer = polymers[i]

        # weight = 0 means the polymer didn't reach this length and is not part of the sample
        if polymer.weight[L - 1] == 0:
            sample_size -= 1

        # convert datatype for numpy calculation
        weight = np.longdouble(polymer.weight[L - 1])

        if observable == 'end-to-end':
            numerator += (weight**2) * (polymer.end_to_end(L) - avg)**2
        else:
            numerator += (weight**2) * (polymer.r_gyration(L) - avg)**2
        
        denominator += weight

    numerator = sample_size * numerator
    denominator = (sample_size-1) * denominator**2

    # avoid division by 0
    if denominator == 0:
        return 0

    error = np.sqrt(numerator/denominator)

    return error

def bootstrap_samples(polymers, n_bootstrap_samples):
    """
    Generates new samples for the error calculation with the bootstrap method.

    polymers : list
        List containing the original sample of polymer objects
    n_bootstrap_samples : int
        Number of "bootstrap" samples
    bootstrap_samples : 2D list
        List containing n_bootstrap_samples number of bootstrap samples lists containing N randomly chosen polymer objects
    """

    N = len(polymers)
    
    # Making n bootstrap samples 
    bootstrap_samples = []
    for n in range(n_bootstrap_samples):
        # Choosing N polymers randomly from "polymers" (S) to form S_i
        polymer_ls_bs = []
        for i in range(N):
            polymer_ls_bs.append(rd.choice(polymers))
        
        bootstrap_samples.append(polymer_ls_bs)

    return bootstrap_samples

def bootstrap_error(bootstrap_samples, n_bootstrap_samples, L, observable):
    '''
    Calculates the error 
    bootstrap_samples : 2D list
        List containing n_bootstrap_samples number of bootstrap samples lists containing N randomly chosen polymer objects
    n_bootstrap_samples : int
        Number of "bootstrap" samples
    L : int
        Length at which the observable average was computed
    observable : str
        The observable of which to compute the error. Only accepted values are: 'end-to-end', 'r_gyration'
    error : float
        Bootstrap error
    '''

    if observable != 'end-to-end' and observable != 'r_gyration':
        print('Misspelled "observable"')
        return 0
    
    # Calculating average for every bootstrap sample using function average
    bootstrap_observable = np.zeros(n_bootstrap_samples, dtype=np.longdouble)
    for n in range(n_bootstrap_samples):
        bootstrap_observable[n] = average(bootstrap_samples[n], observable, L)

    # Calculating error using the estimated standard deviation
    E_X2 = 1/n_bootstrap_samples * np.sum(np.power(bootstrap_observable, 2))
    EX_2 = m.pow(1/n_bootstrap_samples * np.sum(bootstrap_observable), 2)
    error = m.sqrt(E_X2 - EX_2)

    return error


def average_weight(polymer_ls):
    """
    Calculates the average weight of a polymer of length L
    polymer_ls : list
        Contains all Polymer objects
    """

    total_weight = 0
    number_polymers = 0
    for polymer in polymer_ls: 
        if polymer.weight[-1] != 0:
            total_weight += polymer.weight[-1]
            number_polymers += 1
    
    if number_polymers == 0:
        average_weight = 0
    else: 
        average_weight = total_weight / number_polymers

    return average_weight


def PERM(polymer_ls, c_plus, c_minus, prune_count, enrich_count):
    """
    Implementation of the Perm Algorithm
    polymer : list
        Contains all Polymer objects
    c_plus : float
            Predetermined value for c_plus where c_min is founc via c_plus/10
    """
    new_polymer_ls = []
    W_average = average_weight(polymer_ls)
    nr_polymers = 0
    for polymer in polymer_ls:
        if polymer.weight[-1] < c_minus * W_average and polymer.weight[-1] > 0:
            polymer.prune()
            prune_count += 1

        elif polymer.weight[-1] > c_plus * W_average:
            polymer.enrich_weight()
            new_polymer = Polymer()
            new_polymer.enrich_copy(polymer.positions, polymer.weight)
            new_polymer_ls.append(new_polymer)
            enrich_count += 1
            nr_polymers += 1 

        if polymer.weight[-1] != 0:
            nr_polymers += 1 

    polymer_ls = polymer_ls + new_polymer_ls

    return polymer_ls, prune_count, enrich_count, nr_polymers


def make_dir(name):
    """
    Makes directory where all the simulation data is put in
    
    """

    name_folder = "Polymer Simulation Data"
    name_dir = "../" + name_folder + "/" + name

    while os.path.exists(name_dir):
        name_dir = name_dir + " copy"

    os.makedirs(name_dir)

    return name_dir