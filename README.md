# NTT Performance Test

This is a test application for NTT-related algorithms developed by Víctor Gayoso Martínez and David García Lleyda.

USAGE: 

java -jar NTTPerformanceTest.jar [+|-] [N] [I]

  +: negacyclic convolution
  
  -: regular convolution
  
  N: number of coefficients
  
  I: number of iterations

The variable "debug" allows to print out the content of all the polynomials involved in the calculations.  

A sample output for the input "java -jar NTTPerformanceTest.jar - 8 1" is the following:


===== PARAMETERS =====

N: 8

N_INV: 10753

Q: 12289

OMEGA: 4043

OMEGA_INV: 5146

x^(N-1)

-----------------------------------------------

Polynomial a(x):                 [12150, 3263, 7054, 3105, 8575, 8324, 5304, 4583]

Polynomial b(x):                 [3743, 2619, 766, 11047, 279, 3661, 6096, 5555]

Regular NTT(a(x)):               [3202, 9594, 11447, 8719, 1519, 387, 5287, 7889]

Regular NTT(b(x))):              [9188, 2832, 6149, 11346, 291, 4743, 460, 7224]

Regular NTT(a(x)·b(x)):          [110, 11518, 8500, 11613, 11914, 4480, 11087, 6043]

Direct multiplication a(x)·b(x): [6622, 3145, 6033, 11772, 4353, 1122, 1293, 2637]

Plain NTT a(x)·b(x):             [6622, 3145, 6033, 11772, 4353, 1122, 1293, 2637]

Recursive NTT a(x)·b(x):         [6622, 3145, 6033, 11772, 4353, 1122, 1293, 2637]

-----------------------------------------------

END OF COMPUTATION: 1 iteration/s completed

Timing values provided in nanoseconds

-----------------------------------------------

Time NTT Plain Init:              16965500

Time NTT Recursive Init:          34400

Total time Direct:                7800

Total time NTT Plain:             94000

Total time NTT Recursive:         16700


Mean time Direct:                 7800

Mean time NTT Plain:              94000

Mean time NTT Recursive:          16700


Mean time Direct:                 7800

Mean time NTT Plain + Init:       17059500

Mean time NTT Recursive + Init:   51100


Mean time Direct:                 7800

Mean time 1 NTT Plain + Init:     17059500

Mean time 1 NTT Recursive + Init: 51100
