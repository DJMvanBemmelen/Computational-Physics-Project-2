import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
import os


def main():
    # Name of the simulation you want to visualise
    name = 'L_100_N_10000_cplus_3_free'
    PERM_on = True


    # Create path to directory
    name_folder = "Polymer Simulation Data"
    name_dir = "../" + name_folder + "/" + name

    # Reading out main input data
    with open(name_dir + "/input.txt", "r") as input_file:
        lines = input_file.readlines()
        L_max = int(lines[2])
        N = int(lines[4])


    plot_growdata(name, name_dir, PERM_on)
    
    end_to_end = data_reader(name_dir, "end_to_end")
    end_to_end_error = data_reader(name_dir, "end_to_end_error")
    end_to_end_BSerror = data_reader(name_dir, "end_to_end_BSerror")

    title = f'End-to-end of a Free walk'
    plot_fit(L_max, end_to_end, end_to_end_error, title, "end-to-end", "approx analytical", name_dir, PERM_on)
    title = f'End-to-end using the Bootstrapping error'
    plot_fit(L_max, end_to_end, end_to_end_BSerror, title, "end-to-end", "bootstrapping", name_dir, PERM_on)
    
    r_gyration = data_reader(name_dir, "r_gyration")
    r_gyration_error = data_reader(name_dir, "r_gyration_error")
    r_gyration_BSerror = data_reader(name_dir, "r_gyration_BSerror")
    
    title = f'Radius of gyration of a Free walk'
    plot_fit(L_max, r_gyration, r_gyration_error, title, "r_gyration", "approx analytical", name_dir, PERM_on)
    title = f'Radius of gyration using the Bootstrapping error'
    plot_fit(L_max, r_gyration, r_gyration_BSerror, title, "r_gyration", "bootstrapping", name_dir, PERM_on)
    
    return 



def data_reader(name_dir, data):
    """
    Read the data from the .npy file and returns as array for plotting

    Parameters
    ----------
    name_dir : string
        Name of the simulation repository that contains the data.npy files
    data : string
        Name of the data to be returned, without .npy

    Returns
    -------
    data_array : np.ndarray
        The requested data as array
    """
    data_array = np.load(os.path.join(name_dir, data + '.npy'))

    return data_array


def plot_fit(L_max, observable_avr, observable_error, title, observable, error_type, name_dir, PERM_on):
    """
    Plots the inputted data.
    
    Parameters
    ----------
    L_max : int
        Maximum length of the polymers
    observable_avr : ndarray
        Observable data 
    observable_error : ndarray
        Observable error 
    title : string
        The title of the plot
    observable : str
        The observable of which to compute the fit. Only accepted values are: 'end-to-end', 'r_gyration'
    error_type :  str
    """
    x_data = np.arange(1,L_max,1)
    x_label = 'polymer length [subunits]'
    y_label = r'$\langle r^{2} \rangle$'

    fit_y, fit_para, fit_res = expo_fit(observable_avr, observable_error, L_max, observable, PERM_on)


    print(fit_res)
    plt.errorbar(x_data, observable_avr, label='Polymer data', color = "C0")
    plt.fill_between(x_data, observable_avr - observable_error, observable_avr + observable_error, color = "C0", alpha=.4)
    plt.plot(x_data, fit_y, label=r'$\langle r^{2} \rangle = A L^{2\nu}$', color="C1")
    plt.title(title, fontsize=15)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.legend()
    plt.savefig(os.path.join(name_dir, title))
    plt.clf()  

    with open(os.path.join(name_dir, 'input.txt'), "a") as file:
        file.write("Fit parameter A for observable " + observable + " using error_type " + error_type + " \n")
        file.write(f"{fit_para[0]}\n")
        file.write("Fit parameter nu for observable " + observable + " using error_type " + error_type + " \n")
        file.write(f"{fit_para[1]}\n")
        file.write("Residuals for observable " + observable + " using error_type " + error_type + " \n")
        file.write(f"{fit_res}\n")



def expo_fit(observable_avr, observable_error, L_max, observable, PERM_ON):
    """
    Fit the exponential function A * L^nu to the given observable data
    observable_avr : ndarray
        Observable data 
    observable_error : ndarray
        Observable error 
    L_max : int
        Maximum length of the polymers
    observable : string
        Indicates of which observable the average should be calculated. Allowed values: 'end-to-end', 'r_gyration'
    """

    # get observable data
    if observable == 'end-to-end':
        xdata = np.arange(1, L_max, 1)
        ydata = observable_avr
        p0 = [1, 3/4]
        sigma = observable_error
        
        if PERM_ON == False:
            xdata = np.arange(2, L_max, 1)
            ydata = np.delete(observable_avr,0)
            p0 = [1, 3/4]
            sigma = np.delete(observable_error,0)

    elif observable == 'r_gyration':
        xdata = np.arange(2, L_max, 1)
        ydata = np.delete(observable_avr,0)
        p0 = [1, 3/4]
        sigma = np.delete(observable_error,0)
        
        if PERM_ON == False:
            xdata = np.arange(3, L_max, 1)
            ydata = np.delete(observable_avr, [0, 1])
            p0 = [1, 3/4]
            sigma = np.delete(observable_error,[0, 1])

    else:
        print('Misspelled "observable"')
        return 0

    # Exponential function to fit to 
    def f(x, A, nu):
        fx = A * x ** (2 * nu)
        return fx

    # popt, pcov = sc.curve_fit(f, xdata, ydata, p0, sigma)
    popt, pcov = sc.curve_fit(f, xdata, ydata, p0)
    fit_y = f(xdata, popt[0], popt[1])

    if observable == 'r_gyration' and PERM_ON == False:
        fit_y = np.insert(fit_y, [0, 1] , [observable_avr[0], observable_avr[1]])

    if observable == 'end-to-end' and PERM_ON == False:
        fit_y = np.insert(fit_y, 0, observable_avr[0])

    if observable == 'r_gyration' and PERM_ON == True:
        fit_y = np.insert(fit_y, 0, observable_avr[0])

    return fit_y, popt, pcov


def plot_growdata(name, name_dir, PERM_ON):
    """
    Plotting
    """

    if PERM_ON == True:
        # Reading out main input data
        with open(name_dir + "/input.txt", "r") as input_file:
            lines = input_file.readlines()
            L_max = int(lines[2])
            c_plus = float(lines[6])
            prune_count = int(lines[8])
            enrich_count = int(lines[10])

        prune_ls = data_reader(name_dir, "prune_ls")
        enrich_ls = data_reader(name_dir, "enrich_ls")
        nr_polymers_after_ls = data_reader(name_dir, "nr_polymers_after_ls")
        nr_polymers_before_ls = data_reader(name_dir, "nr_polymers_before_ls")

        print("The number of polymers with length 150: ")
        print(nr_polymers_after_ls[-1])

    if PERM_ON == False:
        with open(name_dir + "/input.txt", "r") as input_file:
            lines = input_file.readlines()
            L_max = int(lines[2])

        nr_polymers_before_ls = data_reader(name_dir, "nr_polymers_before_ls")

        print("The number of polymers with length 150: ")
        print(nr_polymers_before_ls[-1])

    if PERM_ON == True:
        print('Total number of prunes: ', prune_count)
        print('Total number of enrichments: ', enrich_count)

        plt.plot(range(0,L_max), prune_ls, label="Nr. prunes")
        plt.plot(range(0,L_max), enrich_ls, label="Nr. enrichments")
        plt.plot(range(0,L_max), nr_polymers_after_ls, label="Nr. polymers")


    plt.plot(range(0,L_max), nr_polymers_before_ls, label="Nr. polymers before")
    plt.xlabel("polymer length [subunits]")
    plt.ylabel("N")
    plt.legend()
    plt.title(f'Grow data')
    plt.savefig(os.path.join(name_dir, "Grow_data" + name))
    plt.clf()


if __name__ == '__main__':
    main()