#!/usr/bin/python

import sys,os
import numpy as np
from read_from_file import read_from_file

class reference_prf:
    def __init__(self, version='1_0'):
        if version == '1_0':

            self.gas_defaults = {
                'N2': [7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,   7.810e-01,
                       7.810e-01,   7.810e-01,   7.810e-01],
                'NH3':[1.940e-16,   1.940e-16,   1.940e-16,   1.940e-16,   1.940e-16,
                       1.940e-16,   1.940e-16,   1.940e-16,   1.940e-16,   1.941e-16,
                       1.967e-16,   2.046e-16,   2.356e-16,   2.506e-16,   2.680e-16,
                       2.854e-16,   3.156e-16,   3.458e-16,   3.852e-16,   4.246e-16,
                       4.787e-16,   5.568e-16,   6.589e-16,   7.629e-16,   9.387e-16,
                       1.177e-15,   1.531e-15,   2.549e-15,   5.739e-15,   2.846e-14,
                       5.213e-13,   3.548e-12,   7.032e-12,   8.621e-12,   1.002e-11,
                       1.137e-11,   1.270e-11,   1.428e-11,   1.748e-11,   2.380e-11,
                       3.129e-11,   5.004e-11,   8.721e-11],
                'OCS':[1.472e-14,   1.472e-14,   1.472e-14,   1.472e-14,   1.472e-14,
                       1.472e-14,   1.472e-14,   1.422e-14,   1.372e-14,   1.321e-14,
                       1.275e-14,   1.232e-14,   1.189e-14,   1.171e-14,   1.155e-14,
                       1.139e-14,   1.192e-14,   1.244e-14,   7.536e-14,   1.383e-13,
                       3.088e-13,   1.686e-12,   5.195e-12,   1.271e-11,   3.741e-11,
                       8.618e-11,   1.629e-10,   2.413e-10,   3.088e-10,   3.652e-10,
                       4.021e-10,   4.243e-10,   4.308e-10,   4.311e-10,   4.310e-10,
                       4.310e-10,   4.310e-10,   4.310e-10,   4.310e-10,   4.310e-10,
                       4.310e-10,   4.310e-10,   4.310e-10],
                'SO2':[3.880e-11,   3.880e-11,   3.880e-11,   3.880e-11,   3.880e-11,
                       3.880e-11,   3.880e-11,   3.879e-11,   3.880e-11,   3.884e-11,
                       3.822e-11,   3.693e-11,   3.300e-11,   3.108e-11,   2.801e-11,
                       2.494e-11,   2.051e-11,   1.608e-11,   1.513e-11,   1.417e-11,
                       1.534e-11,   1.733e-11,   1.969e-11,   2.245e-11,   2.681e-11,
                       3.157e-11,   3.715e-11,   4.498e-11,   5.454e-11,   6.813e-11,
                       9.056e-11,   1.324e-10,   1.975e-10,   2.271e-10,   2.493e-10,
                       2.666e-10,   2.794e-10,   2.850e-10,   2.856e-10,   2.895e-10,
                       2.915e-10,   2.909e-10,   2.910e-10],
                'HI': [3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,   3.000e-12,
                       3.000e-12,   3.000e-12,   3.000e-12],
                'HONO': [1.000e-14,   1.000e-14,   1.000e-14,   1.000e-14,   1.000e-14,
                         1.000e-14,   1.000e-14,   3.095e-14,   1.043e-13,   2.564e-13,
                         8.372e-13,   1.767e-12,   4.118e-12,   5.330e-12,   7.275e-12,
                         9.220e-12,   1.016e-11,   1.110e-11,   9.660e-12,   8.220e-12,
                         6.300e-12,   4.620e-12,   3.070e-12,   1.760e-12,   9.500e-13,
                         5.110e-13,   2.780e-13,   1.680e-13,   1.150e-13,   8.960e-14,
                         1.030e-13,   1.000e-13,   9.500e-14,   1.040e-13,   1.130e-13,
                         1.770e-13,   2.410e-13,   3.220e-13,   4.030e-13,   4.330e-13,
                         4.630e-13,   4.630e-13,   4.630e-13],
                'CH3F': [5.000e-14,   5.000e-14,   5.000e-14,   5.000e-14,   5.000e-14,
                         5.000e-14,   5.000e-14,   5.000e-14,   5.000e-14,   5.000e-14,
                         5.000e-14,   5.000e-14,   5.000e-14,   5.000e-14,   5.005e-14,
                         5.010e-14,   5.065e-14,   5.120e-14,   5.410e-14,   5.700e-14,
                         6.100e-14,   6.700e-14,   7.650e-14,   9.230e-14,   1.120e-13,
                         1.500e-13,   1.990e-13,   2.420e-13,   2.870e-13,   3.350e-13,
                         3.750e-13,   3.870e-13,   3.970e-13,   4.000e-13,   4.000e-13,
                         4.000e-13,   4.000e-13,   4.000e-13,   4.000e-13,   4.000e-13,
                         4.000e-13,   4.000e-13,   4.000e-13],
                'CF4' :[3.000e-11,   3.000e-11,   3.000e-11,   3.000e-11,   3.000e-11,
                        3.000e-11,   3.000e-11,   3.000e-11,   3.023e-11,   3.080e-11,
                        3.492e-11,   3.989e-11,   4.641e-11,   5.000e-11,   5.000e-11,
                        5.000e-11,   5.350e-11,   5.700e-11,   6.100e-11,   6.500e-11,
                        6.500e-11,   6.500e-11,   6.500e-11,   6.500e-11,   6.500e-11,
                        6.500e-11,   6.500e-11,   6.750e-11,   7.000e-11,   7.000e-11,
                        7.000e-11,   7.000e-11,   7.000e-11,   7.000e-11,   7.000e-11,
                        7.000e-11,   7.000e-11,   7.000e-11,   7.000e-11,   7.000e-11,
                        7.000e-11,   7.000e-11,   7.000e-11],
                'CCL2F2':[2.899e-21,   1.844e-20,   3.835e-19,   1.220e-18,   4.347e-18,
                         1.540e-17,   4.002e-17,   8.649e-17,   2.034e-16,   5.494e-16,
                         1.694e-15,   6.200e-15,   2.714e-14,   4.978e-14,   8.815e-14,
                         1.514e-13,   2.664e-13,   5.237e-13,   1.154e-12,   2.675e-12,
                         5.954e-12,   1.186e-11,   2.151e-11,   3.725e-11,   6.003e-11,
                         9.000e-11,   1.307e-10,   1.897e-10,   2.649e-10,   3.460e-10,
                         4.024e-10,   4.358e-10,   4.672e-10,   4.753e-10,   4.788e-10,
                         4.801e-10,   4.805e-10,   4.806e-10,   4.806e-10,   4.804e-10,
                         4.802e-10,   4.800e-10,   4.798e-10],
                'CCL3F':[2.061e-21,   2.064e-21,   2.070e-21,   2.073e-21,   2.074e-21,
                         2.075e-21,   2.078e-21,   2.080e-21,   2.083e-21,   2.085e-21,
                         2.087e-21,   2.090e-21,   2.094e-21,   2.101e-21,   2.148e-21,
                         2.490e-21,   5.989e-21,   5.975e-20,   9.694e-19,   1.749e-17,
                         3.742e-16,   6.012e-15,   5.818e-14,   4.438e-13,   2.068e-12,
                         7.100e-12,   1.912e-11,   4.419e-11,   8.411e-11,   1.357e-10,
                         1.756e-10,   2.005e-10,   2.242e-10,   2.303e-10,   2.329e-10,
                         2.339e-10,   2.342e-10,   2.343e-10,   2.344e-10,   2.344e-10,
                         2.343e-10,   2.342e-10,   2.341e-10],
                'COCLF':[3.750e-17,   3.750e-17,   3.750e-17,   3.750e-17,   3.750e-17,
                         3.750e-17,   3.750e-17,   3.750e-17,   7.000e-17,   9.500e-17,
                         1.000e-16,   1.000e-16,   9.400e-16,   1.780e-15,   7.240e-15,
                         1.270e-14,   6.285e-14,   1.130e-13,   5.815e-13,   1.050e-12,
                         2.900e-12,   7.000e-12,   1.370e-11,   2.170e-11,   2.810e-11,
                         2.950e-11,   2.580e-11,   1.890e-11,   1.110e-11,   4.770e-12,
                         1.370e-12,   4.970e-13,   3.570e-13,   3.140e-13,   2.710e-13,
                         2.380e-13,   2.040e-13,   1.770e-13,   1.500e-13,   1.280e-13,
                         1.060e-13,   8.880e-14,   7.175e-14],
                'CHF2CL':[1.947e-13,   1.029e-12,   6.029e-12,   1.112e-11,   1.933e-11,
                          2.707e-11,   3.131e-11,   3.386e-11,   3.552e-11,   3.661e-11,
                          3.733e-11,   3.788e-11,   3.883e-11,   3.959e-11,   4.069e-11,
                          4.223e-11,   4.420e-11,   4.654e-11,   4.917e-11,   5.206e-11,
                          5.525e-11,   5.877e-11,   6.268e-11,   6.724e-11,   7.219e-11,
                          7.747e-11,   8.389e-11,   9.264e-11,   1.035e-10,   1.158e-10,
                          1.249e-10,   1.307e-10,   1.375e-10,   1.396e-10,   1.406e-10,
                          1.411e-10,   1.412e-10,   1.413e-10,   1.413e-10,   1.413e-10,
                          1.414e-10,   1.416e-10,   1.418e-10],
                'COCL2':[1.000e-16,   1.000e-16,   1.000e-16,   1.000e-16,   1.000e-16,
                         1.000e-16,   1.000e-16,   1.000e-16,   1.000e-16,   1.000e-16,
                         1.000e-16,   1.000e-16,   1.000e-16,   1.000e-16,   1.000e-16,
                         1.000e-16,   4.100e-15,   8.100e-15,   1.190e-13,   2.300e-13,
                         9.340e-13,   3.090e-12,   7.780e-12,   1.510e-11,   2.260e-11,
                         2.630e-11,   2.440e-11,   1.850e-11,   1.110e-11,   4.860e-12,
                         1.410e-12,   5.130e-13,   3.690e-13,   3.240e-13,   2.800e-13,
                         2.460e-13,   2.110e-13,   1.830e-13,   1.550e-13,   1.320e-13,
                         1.100e-13,   9.200e-14,   7.456e-14],
                'CH3I': [1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,   1.000e-17,
                         1.000e-17,   1.000e-17,   1.000e-17],
                'H2S':[1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,   1.500e-11,
                       1.500e-11,   1.500e-11,   1.500e-11],
                'CHCL2F':[4.520e-13,   4.790e-13,   5.060e-13,   5.330e-13,   5.600e-13,
                         5.870e-13,   6.130e-13,   6.400e-13,   6.670e-13,   6.940e-13,
                         7.210e-13,   7.320e-13,   7.420e-13,   7.530e-13,   7.640e-13,
                         7.740e-13,   7.850e-13,   1.590e-12,   2.420e-12,   3.270e-12,
                         4.140e-12,   5.030e-12,   5.950e-12,   6.880e-12,   7.840e-12,
                         8.820e-12,   9.820e-12,   1.080e-11,   1.190e-11,   1.290e-11,
                         1.370e-11,   1.400e-11,   1.420e-11,   1.430e-11,   1.440e-11,
                         1.440e-11,   1.450e-11,   1.460e-11,   1.470e-11,   1.480e-11,
                         1.480e-11,   1.450e-11,   1.460e-11],
                'SF6':[4.706e-13,   4.706e-13,   4.706e-13,   4.706e-13,   4.706e-13,
                       4.706e-13,   4.706e-13,   4.874e-13,   5.068e-13,   5.302e-13,
                       5.674e-13,   6.005e-13,   6.333e-13,   6.642e-13,   6.830e-13,
                       7.017e-13,   7.008e-13,   7.000e-13,   7.020e-13,   7.039e-13,
                       7.917e-13,   7.992e-13,   8.004e-13,   7.962e-13,   8.861e-13,
                       8.999e-13,   8.963e-13,   9.867e-13,   1.003e-12,   1.094e-12,
                       1.104e-12,   1.191e-12,   1.286e-12,   1.370e-12,   1.423e-12,
                       1.460e-12,   1.504e-12,   1.500e-12,   1.499e-12,   1.500e-12,
                       1.500e-12,   1.500e-12,   1.500e-12],
                'F134A':[1.249e-12,   1.249e-12,   1.249e-12,   1.249e-12,   1.249e-12,
                         1.249e-12,   1.249e-12,   1.270e-12,   1.288e-12,   1.302e-12,
                         1.311e-12,   1.316e-12,   1.316e-12,   1.316e-12,   1.311e-12,
                         1.305e-12,   1.300e-12,   1.294e-12,   1.288e-12,   1.282e-12,
                         1.294e-12,   1.371e-12,   1.393e-12,   1.427e-12,   1.572e-12,
                         1.773e-12,   1.996e-12,   2.297e-12,   2.676e-12,   3.524e-12,
                         4.505e-12,   5.442e-12,   6.824e-12,   7.606e-12,   8.153e-12,
                         8.522e-12,   8.766e-12,   8.839e-12,   9.031e-12,   9.186e-12,
                         9.342e-12,   9.504e-12,   9.707e-12],
                'F142B':[1.392e-12,   1.392e-12,   1.392e-12,   1.392e-12,   1.392e-12,
                         1.392e-12,   1.392e-12,   1.415e-12,   1.435e-12,   1.451e-12,
                         1.461e-12,   1.467e-12,   1.467e-12,   1.467e-12,   1.461e-12,
                         1.454e-12,   1.448e-12,   1.442e-12,   1.435e-12,   1.429e-12,
                         1.442e-12,   1.528e-12,   1.553e-12,   1.591e-12,   1.752e-12,
                         1.976e-12,   2.225e-12,   2.561e-12,   2.983e-12,   3.928e-12,
                         5.021e-12,   6.065e-12,   7.606e-12,   8.478e-12,   9.087e-12,
                         9.498e-12,   9.770e-12,   9.852e-12,   1.007e-11,   1.024e-11,
                         1.041e-11,   1.059e-11,   1.082e-11],
                'F141B':[1.448e-12,   1.448e-12,   1.448e-12,   1.448e-12,   1.448e-12,
                         1.448e-12,   1.448e-12,   1.472e-12,   1.493e-12,   1.510e-12,
                         1.520e-12,   1.526e-12,   1.526e-12,   1.526e-12,   1.520e-12,
                         1.513e-12,   1.507e-12,   1.500e-12,   1.493e-12,   1.487e-12,
                         1.500e-12,   1.590e-12,   1.616e-12,   1.655e-12,   1.823e-12,
                         2.056e-12,   2.315e-12,   2.664e-12,   3.103e-12,   4.086e-12,
                         5.223e-12,   6.310e-12,   7.913e-12,   8.820e-12,   9.454e-12,
                         9.882e-12,   1.016e-11,   1.025e-11,   1.047e-11,   1.065e-11,
                         1.083e-11,   1.102e-11,   1.125e-11],
                'CHF3':[2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,   2.750E-11,
                        2.750E-11,   2.750E-11,   2.750E-11]}

            self.gases = ['H2O',     'CO2',     'O3',     'N2O',     'CO',      #1-5
                          'CH4',     'O2',      'NO',     'SO2',     'NO2',     #6-10 
                          'NH3',     'HNO3',    'OH',     'HF',      'HCL',     #11-15
                          'HBR',     'HI',      'CLO',    'OCS',     'CH2O',    #16-20 
                          'HOCL',    'HO2',     'H2O2',   'HONO',    'HO2NO2',  #21-25
                          'N2O5',    'CLONO2',  'HCN',    'CH3F',    'CH3CL',   #26-30
                          'CF4',     'CCL2F2',  'CCL3F',  'CH3CCL3', 'CCL4',    #31-35
                          'COF2',    'COCLF',   'C2H6',   'C2H4',    'C2H2',    #36-40 
                          'N2',      'CHF2CL',  'COCL2',  'CH3BR',   'CH3I',    #41-45
                          'HCOOH',   'H2S',     'CHCL2F', 'O2CIA',   'SF6',     #46-50
                          'NF3',     'N2CIA',   'OTHER',  'OTHER',   'PH3',     #51-55
                          'OTHER',   'OTHER',   'OCLO',   'F134A',   'C3H8',    #56-60
                          'F142B',   'CFC113',  'F141B',  'CH3OH',   'OTHER',   #61-65
                          'OTHER',   'PAN',     'CH3CHO' ,'CH3CN',   'CHF3',   #66-70
                          'CH3COOH', 'C5H8',    'MVK',    'MACR',    'C3H6',    #71-75
                          'C4H8',    'OTHER',   'OTHER',  'OTHER',   'OTHER',   #76-80
                          'OTHER',   'OTHER',   'OTHER',  'OTHER',   'OTHER',   
                          'OTHER',   'OTHER',   'OTHER',  'OTHER',   'OTHER',    
                          'OTHER',   'OTHER',   'OTHER',  'OTHER',   'OTHER',   
                          'OTHER',   'OTHER',   'OTHER',  'OTHER']


    def create_ref_from_WACCM(self, direc,default=None):
        # inserts WACCM output in reference. gasnames musst be defined already
        # interpolates to the altitude grid already defined in the class

        if (default != None):
            defref = reference_prf()
            defref.read_reference_prf(default)
        else:
            defref = -1
            
            
        self.nr_gas = len(self.gases)
        self.gas_nr = np.zeros(self.nr_gas)
        self.gasname = []
        self.notes = []
        alts, vmr, note = self.load_waccmfile('%s/%s.txt'%(direc,'T'))
        self.nr_layers = len(alts)
        self.t = np.array(vmr)
        self.z = np.array(alts)
        alts, vmr, note = self.load_waccmfile('%s/%s.txt'%(direc,'P'))
        self.p = np.array(vmr)
        self.vmr = np.zeros((self.nr_gas, self.nr_layers))
        for nr in range(0,self.nr_gas):
            waccm_file = '%s/%s.txt'%(direc,self.gases[nr])
            self.gas_nr[nr] = nr + 1
            if self.gases[nr] == 'CH2O':
                self.gases[nr] = 'H2CO'
            self.gasname.append(self.gases[nr])
            self.notes.append('created by Create_REF_from_WACCM')
            if os.path.exists(waccm_file):
                print (waccm_file)
                alts, vmr, note = self.load_waccmfile(waccm_file)
                self.vmr[nr,0:self.nr_layers] = vmr
                self.notes[-1] = note
            else:
                if (defref != -1):
                    self.vmr[nr,0:self.nr_layers] = defref.get_gas_from_refprofile(self.gas_defaults[nr], z=self.z)
                    self.notes[-1] = 'copied from {}'.format(default)
                else:
                    if self.gases[nr] in self.gas_defaults:
                        self.vmr[nr,:] = self.gas_defaults[self.gases[nr]]
                    else:
                        self.vmr[nr,:] = 0.0
                    self.notes[-1] = 'default set from Create_REF_from_WACCM.py'

            self.gas_nr[nr] = nr+1

        # Fill in defaults
        for nr in range(0,self.nr_gas):
            if self.gases[nr] == 'O2CIA':
                nr2 = self.gasname.index('O2')
                self.vmr[nr] = self.vmr[nr2]
                self.notes[nr] = 'set equal to O2 profile'
            if self.gases[nr] == 'N2CIA':
                nr2 = self.gasname.index('N2')
                self.vmr[nr] = self.vmr[nr2]
                self.notes[nr] = 'set equal to N2 profile'
                    

                    
    def load_waccmfile(self, filename):
        # loads the second block formatted for reference.prf from the respective file
        wfile = read_from_file(filename)
        wfile.skipline()
        nr = wfile.next(1).pop()
        wfile.next(nr*5)
        nr = wfile.next(1).pop()
        wfile.skipline()
        alt = np.array(wfile.next(nr))
        line = wfile.get_line().strip()
        gas,tmp,note = line.partition(' ')
        vmr = wfile.next(nr)
        return(alt, vmr, note)

    def write_reference_prf(self, prffile, nopt=False):
        
        fid = open(prffile, 'w')
        line = '%5d%6d%6d\n' % (1,self.nr_layers,self.nr_gas)
        fid.write(line)
        line = '    ALTITUDE\n'
        fid.write(line)
        for n in range(0,self.nr_layers,5):
            line = ''.join(map (lambda x:'%11.3e,'%x, self.z[n:np.min((n+5,self.nr_layers))]))
            fid.write(' '+line+'\n')
        if not nopt:
            line = '    PRESSURE\n'
            fid.write(line)
            for n in range(0,self.nr_layers,5):
                line = ''.join(map (lambda x:'%11.3e,'%x, self.p[n:np.min((n+5,self.nr_layers))]))
                fid.write(' '+line+'\n')
            line = '    TEMPERATURE\n'
            fid.write(line)
            for n in range(0,self.nr_layers,5):
                line = ''.join(map (lambda x:'%11.3e,'%x, self.t[n:np.min((n+5,self.nr_layers))]))
                fid.write(' '+line+'\n')

        for gas in range(0,self.nr_gas):
            line = '%5d%8s %s' %(self.gas_nr[gas],self.gasname[gas],self.notes[gas])
            fid.write(line+'\n')
            for n in range(0,self.nr_layers,5):
                line = ''.join(map (lambda x:'%11.3e,'%x, 
                                        self.vmr[gas,n:np.min((n+5,self.nr_layers))]))
                fid.write(' '+line+'\n')


        fid.close()

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print ('Usage: Create_REF_from_WACCM.py <path to WACCM profiles> <optional: result.prf>')
        exit()
    if len(sys.argv) < 3:
        ref_file = 'waccm.prf'
    else:
        ref_file = sys.argv[2]
        
    waccm_dir = sys.argv[1]
    sr = reference_prf()
    sr.create_ref_from_WACCM(waccm_dir)
    sr.write_reference_prf(ref_file)

else:
    print ('This is a script to be run from commandline only')
