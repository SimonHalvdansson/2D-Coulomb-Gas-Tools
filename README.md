# Lemniscate-Polynomial-Data

Each zip file contains a .pckl file, see e.g. https://docs.python.org/3/library/pickle.html for more on this file format. In short it stores a Python object as a file. For useage see the runner.py file. The short version is that setting `ps = load('ps600')` makes 'ps' a 600 times 600 matrix with coefficients so the 40th coefficient for p_300 is `ps[300][40]`.
