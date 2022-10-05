# SteadyState_BH
Code using [AMUSE](https://amuse.readthedocs.io/en/latest/) which looks at simulating intermediate mass black hole (IMBH) clusters, testing their stability in an attempt to investigate whether IMBH mergers could help explain the formation of super massive black holes.

To run the code, run interface.py while being in the SteadyState_BH/Code directory.

When using ALICE supercomputer, change number\_of\_workers on line 38 in evol.py, uncomment out lines 304 - 310 in evol.py and change file name on line 44 in data_init.py.
