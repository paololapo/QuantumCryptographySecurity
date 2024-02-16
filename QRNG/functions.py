import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
from scipy import stats

def group_events(mat_file):
    """Group events by channel and compute the time differences"""
    #Create a Pandas DataFrame
    time = mat_file["timetag"][0]
    channel = mat_file["channel"][0]
    df = pd.DataFrame({"channel" : channel, "time" : time})
    print("Analyzing " + str(len(time)) + " events...")

    #Select only the events between channels (1, 2) and (1, 3)
    ch2 = df[df["channel"] != 3]
    ch3 = df[df["channel"] != 2]

    #Compute the differences row with the previous
    ch2 = ch2.diff()
    ch3 = ch3.diff()

    return (ch2, ch3)


def find_coincidences(events, cut):
    """Find the coincidences between the two channels"""
    ch2, ch3 = events[0], events[1]

    #Select only the coincidences (clicks in different channels)
    ch2 = ch2[ch2["channel"]!=0]
    ch3 = ch3[ch3["channel"]!=0]

    #Cut the tail
    ch2 = ch2[ch2["time"] < cut]
    ch3 = ch3[ch3["time"] < cut]

    return (ch2, ch3)


def filter_coincidences(coincidences, prec):
    """Align the time tagger and return the coincidences"""
    ch2, ch3 = coincidences[0], coincidences[1]

    #Gaussian fit (storing only the averages)
    m2 = stats.norm.fit(ch2.time)[0]
    m3 = stats.norm.fit(ch3.time)[0]
    
    #Coincidences within precision
    coinc2 = ch2[np.abs((ch2["time"]-m2)) < prec]
    coinc3 = ch3[np.abs((ch3["time"]-m3)) < prec]
    print("Found " + str(len(coinc2)) + " coincidences on channel 2.")
    print("Found " + str(len(coinc3)) + " coincidences on channel 3.")

    return (coinc2, coinc3)


def find_probabilities_from_coincidences(filter_coincidences):
    """Compute the probabilities from the coincidences"""
    coinc2, coinc3 = filter_coincidences[0], filter_coincidences[1]
    
    #Count coincidences
    n2 = len(coinc2)
    n3 = len(coinc3)
    N = n2 + n3

    #Return (frequentist) probabilities
    p2 = n2/N
    p3 = n3/N
    print("Probabilities: " + str((p2, p3)))
    return p2, p3


def find_probabilities_from_mat_file(mat_file, prec, cut=10000):
    """Compute the probabilities directly from the raw file"""
    events = group_events(mat_file)
    coinc = find_coincidences(events, cut)
    coinc = filter_coincidences(coinc, prec)
    prob = find_probabilities_from_coincidences(coinc)
    return prob


def min_entropy(prob_x):
    return -np.log2(max(prob_x))

def max_entropy(prob_x):
    prob_x = np.array(list(prob_x))
    return 2*np.log2(np.sum(prob_x**0.5))


def stokes_parameters(px, py=[0, 0], pz=[0, 0]):
    s1 = px[0] - px[1]
    s2 = py[0] - py[1]
    s3 = pz[0] - pz[1]
    return s1, s2, s3


def density_matrix(strokes_parameters):
    s1, s2, s3 = strokes_parameters[0], strokes_parameters[1], strokes_parameters[2]
    rho_real = 0.5*np.array([[1+s3, s1], [s1, 1-s3]])
    rho_imag = 0.5*np.array([[0, -s2], [s2, 0]])
    return rho_real+1j*rho_imag


def f_rho(stokes_parameters):
    s1, s2, s3 = stokes_parameters
    k = s1-s2*1j
    mod2 = k.imag**2 + k.real**2
    par = (1 + np.sqrt(1-mod2))/2
    return -np.log2(par)


def security_parameter(l, h):
    par = 0.5*np.sqrt(2**(l-h))
    return par


    