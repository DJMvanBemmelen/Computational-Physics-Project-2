# from concurrent.futures.process import _ExecutorManagerThread
# from pyrsistent import b
# from sympy import Q
from class_polymer import Polymer
import execution as ex

from matplotlib import pyplot as plt
import numpy as np
import time



def main():
    
    L_max = 150
    N = 10000
    n_bootstrap_samples = 0

    # no perm
    c_plus = 3
    name = "L_100_N_10000_cplus_3_free"
    name_dir = ex.make_dir(name)
    polymer_ls = ex.generate_free(L_max, N, c_plus, name_dir)
    ex.observables(L_max, polymer_ls, n_bootstrap_samples, name_dir)

    print('\nSimulation 1 finished\n')

    # # ## Draw Polymers ###
    # # ex.draw_polymer(polymer_ls[0], len(polymer_ls[0].positions))

    # # ## Plot growth direction distribution ###
    # # ex.plot_directions(directions)



if __name__ == '__main__':
    main()