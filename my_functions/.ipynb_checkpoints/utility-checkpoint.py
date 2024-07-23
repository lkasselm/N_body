# general utility functions
import numpy as np

def f_of_x(f, X, x, interpolation = None):
    """
    given an array X and an array f, where f is supposed to be a function of X,
    find f(x) assuming that it is in the bounds of X

    Inputs:
    f (array)
    X (array)
    x (float)
    """
    X = np.nan_to_num(X, nan=0) # set nan entries to zero
    idx = np.argmin(np.abs(X-x)) # get index of X, whose value is closest to x
    if interpolation is None:
        return f[idx]
    elif interpolation == 'linear':
        # check edge cases
        if x < X[idx]:
            if idx == 0:
                return f[idx]
            else:
                return f[idx -1] + (x-X[idx-1])*(f[idx] - f[idx-1])/(X[idx] - X[idx-1])
        else:
            # same with 1 idx higher
            if idx == len(X)-1:
                return f[idx]
            else:
                return f[idx] + (x-X[idx])*(f[idx+1] - f[idx])/(X[idx+1] - X[idx])
    else:
        raise NameError('Unknown interpolation method')

def print_in_scientific_notation(input):
    """
    Inputs:
    * A floating point number
    or:
    * A tuple consisting of strings and floating point numbers

    Prints the numbers in scientific (base-10) notation
    """
    if type(input) is float:
        print(f"{number:.2e}")
    elif type(input) is tuple:
        out = ''
        for i in input:
            if type(i) is float:
                out = out + f"{i:.2e}"
            elif type(i) is str:
                out = out + i
            else:
                raise TypeError('Unsupported type', type(i))
    else:
        raise TypeError('Unsupported type', type(input))
    print(out)
    
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
    
def append_to_new_line(filename, text):
    # Open the file in append mode
    with open(filename, 'a') as file:
        # Add a newline character before the text if the file isn't empty
        file.write('\n' + text)