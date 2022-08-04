#H2 partition function from (Popovas & Jørgensen 2016)

def PartitionFunction(Temp):
    from scipy import interpolate
    import numpy as np
    T = np.array([500,600,700,800,900,1000,1100,1200,1300,1400,1500,2000,2500,3000,3500,4000,4500,5000,6000,7000])
    Q = np.array([3.128,3.725,4.324,4.928,5.536,6.151, 6.774, 7.406, 8.051, 8.709, 9.382,13.015,17.181,21.964,27.428,33.632,40.634,48.492,66.988,89.448 ])
    tck = interpolate.splrep(T, Q, k=3)
    PartFunc = interpolate.splev(Temp, tck, der=0)
    return PartFunc