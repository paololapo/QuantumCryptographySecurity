import numpy as np
import math


# Binary entropy
def h2(p):
    assert p > 0, "While computing h2(p), p must be greater than 0"
    return -p*np.log2(p)-(1-p)*np.log2(1-p)


def read_file(path, max_keys=None):
    """Read the file and return the keys"""
    block_lengths = []
    keys = []

    # Get the bytes
    bytes = np.fromfile(path, dtype="uint8")

    i = 0
    j = 0
    while (i < len(bytes)):
        # Get the length of the key block
        N = int.from_bytes(bytes[i:i+8], byteorder="big")
        block_lengths.append(N)
        i += 8

        # Get the key
        N_bytes = bytes[i:i+N]
        key = np.concatenate([np.unpackbits(b, bitorder="little") for b in N_bytes])
        keys.append(key)
        i += N
        
        j += 1
        if max_keys != None and j >= max_keys: break
    
    return block_lengths, keys


def key_decoding(key, mode):
    """Decode the key (given convention)"""
    assert mode == "AB" or mode == "D", "Mode issue"

    # Split the key into blocks
    k = np.split(key, len(key)//2)

    # Define the convertion table
    def state_convertion(pair, mode):
        if mode == "AB":
            if (pair == [0, 0]).all(): return "H"
            if (pair == [0, 1]).all(): return "V"
            if (pair == [1, 0]).all(): return "D"
            if (pair == [1, 1]).all(): return "A"
        
        if mode == "D":
            if (pair == [0, 0]).all(): return "S"
            if (pair == [0, 1]).all(): return "L"
            if (pair == [1, 0]).all(): return "Unused"
            if (pair == [1, 1]).all(): return "Unused"

    # Decode the key
    decoded = [state_convertion(pair, mode) for pair in k]

    return decoded


def basis_reconciliation(key_A, key_B, decoy_states, nZ_target=None):
    """Sift the raw keys"""
    assert len(key_A) == len(key_B) == len(decoy_states), "size problem"

    # Get the intensities and the indexing
    #intensities = list(set(decoy_states))
    intensities = ['S', 'L']
    intensity_index = {element: index for index, element in enumerate(intensities)}

    X_A = [[] for _ in intensities]
    Z_A = [[] for _ in intensities]
    X_B = [[] for _ in intensities]
    Z_B = [[] for _ in intensities]

    pulses = 0

    for i in range(len(key_A)):
        if nZ_target != None and (len(Z_A[0]) + len(Z_A[1])) >= nZ_target:
            break
        
        measure_A = key_A[i]
        measure_B = key_B[i]
        k = decoy_states[i]

        # Basis reconciliation
        if measure_A in ["H", "V"] and measure_B in ["H", "V"]:
            Z_A[intensity_index[k]].append(measure_A)
            Z_B[intensity_index[k]].append(measure_B)

        if measure_A in ["A", "D"] and measure_B in ["A", "D"]:
            X_A[intensity_index[k]].append(measure_A)
            X_B[intensity_index[k]].append(measure_B)
        
        pulses += 1
    
    return X_A, Z_A, X_B, Z_B, pulses


def errors_estimation(b_A, b_B):
    """Count the mismatches for every key"""
    m_b = [len([_ for i in range(len(b_A[k])) if b_A[k][i] != b_B[k][i]]) for k in range(len(b_A))]
    return m_b


def secret_key_length(parameters, eps1, eps2, eps_sec, eps_cor, lambda_EC):  
    """Compute all the quantities for the secret key length estimation"""
    assert len(parameters) == 6, "Expected 6 parameters"
    
    k_list = parameters[0]
    p_k_list = parameters[1]
    nZ_k_list = parameters[2]
    mZ_k_list = parameters[3]
    nX_k_list = parameters[4]
    mX_k_list = parameters[5]
    
    assert len(k_list) == len(p_k_list) == len(nZ_k_list) == len(mZ_k_list) == len(mX_k_list) == 2, "lists should have 2 elements"
    assert k_list[0] > k_list[1], "mu_1 is expected to be greater than mu_2"

    # Helper function
    def delta(a, e):
        return np.sqrt(0.5*a*np.log(1/e))


    # Helper parameters of the block lengths
    def nZ_k_pm(k, p_k, nZ_k, nZ, eps1):
        first_term = np.exp(k)/p_k
        #first_term = 1
        sqrt = delta(nZ, eps1)

        n_p = first_term*(nZ_k + sqrt)
        n_m = first_term*(nZ_k - sqrt)

        #print("nZ_k_pm", n_p, n_m)
        return n_p, n_m
    

    # Helper parameters of the error number
    def mZ_k_pm(k, p_k, mZ_k, mZ, eps2):
        first_term = np.exp(k)/p_k
        #first_term = 1
        sqrt = delta(mZ, eps2)

        m_p = first_term*(mZ_k + sqrt)
        m_m = first_term*(mZ_k - sqrt)

        #print("mZ_k_pm", m_p, m_m)
        return m_p, m_m
    

    # Probability of a n-photons state  
    def tau_n(n, k_list, p_k_list):
        r = 0
        for k, p_k in zip(k_list, p_k_list):
            r += np.exp(-k)*(k**n)*p_k/math.factorial(n)

        return r
    

    # Lower bound to the vacuum events
    def sZ_0_low(k_list, p_k_list, nZ_k_list, eps1):
        mu_1 = k_list[0]
        mu_2 = k_list[1]
        
        first_term = tau_n(0, k_list, p_k_list)/(mu_1-mu_2)
        nZ_2_m = nZ_k_pm(mu_2, p_k_list[1], nZ_k_list[1], np.sum(nZ_k_list), eps1)[1]
        nZ_1_p = nZ_k_pm(mu_1, p_k_list[0], nZ_k_list[0], np.sum(nZ_k_list), eps1)[0]
        second_term = mu_1*nZ_2_m - mu_2*nZ_1_p

        #print("sZ_0_low", first_term*second_term)
        result = np.max([first_term*second_term, 0])
        #if result == 0: print("Set sZ_0_low = 0")
        return result
    

    # Upper bound to the vacuum events (*)
    def sZ_0_up(nZ_k_list, mZ_k_list, eps1):
        # First bound
        first_bound = 2*(np.sum(mZ_k_list) + delta(np.sum(nZ_k_list), eps1))

        return first_bound
    
    
    # Lower bound to the single-photon events
    def sZ_1_low(k_list, p_k_list, nZ_k_list, mZ_k_list, eps1):
        mu1 = k_list[0]
        mu2 = k_list[1]

        term1 = tau_n(1, k_list, p_k_list)*mu1/mu2/(mu1-mu2)
        term2 = nZ_k_pm(mu2, p_k_list[1], nZ_k_list[1], np.sum(nZ_k_list), eps1)[1]
        term3 = (mu2**2)*nZ_k_pm(mu1, p_k_list[0], nZ_k_list[0], np.sum(nZ_k_list), eps1)[0]/mu1**2
        term4 = (mu1**2-mu2**2)*sZ_0_up(nZ_k_list, mZ_k_list, eps1)/mu1**2/tau_n(0, k_list, p_k_list)

        #print("sZ_1_low", term1*(term2-term3-term4))
        return np.max([term1*(term2-term3-term4), 0])
        

    # Helper function
    def gamma(a, b, c, d):
        x = (c+d)*(1-b)*b/(c*d*np.log(2))
        y = np.log2(441*(c+d)/c/d/(1-b)/b/a**2)
        return np.sqrt(x*y)
    

    # Upper bound to the error in X due to single-photons
    def vX_1_up(k_list, p_k_list, mX_k_list):
        first_term = tau_n(1, k_list, p_k_list)/(k_list[0]-k_list[1])
        mX_mu1_p = mZ_k_pm(k_list[0], p_k_list[0], mX_k_list[0], np.sum(mX_k_list), eps2)[0]
        mX_mu2_m = mZ_k_pm(k_list[1], p_k_list[1], mX_k_list[1], np.sum(mX_k_list), eps2)[1]

        #print("vX_1_up", first_term*(mX_mu1_p-mX_mu2_m))
        return first_term*(mX_mu1_p-mX_mu2_m)
    

    # Phase error rate in the Z basis
    def phi_x_up(k_list, p_k_list, nZ_k_list, mZ_k_list, nX_k_list, mX_k_list, eps1, eps_sec):
        aZ = sZ_1_low(k_list, p_k_list, nZ_k_list, mZ_k_list, eps1)
        aX = sZ_1_low(k_list, p_k_list, nX_k_list, mX_k_list, eps1)
        term1 = vX_1_up(k_list, p_k_list, mX_k_list)/aX
        
        #print("phi_x_up", term1 + gamma(eps_sec, term1, aZ, aX))
        return term1 + gamma(eps_sec, term1, aZ, aX)
    

    phi = phi_x_up(k_list, p_k_list, nZ_k_list, mZ_k_list, nX_k_list, mX_k_list, eps1, eps_sec)
    if not np.isnan(phi) and phi>0:
        t1 = sZ_0_low(k_list, p_k_list, nZ_k_list, eps1)
        t2 = sZ_1_low(k_list, p_k_list, nZ_k_list, mZ_k_list, eps1)
        t3 = 1-h2(phi)
        t4 = lambda_EC + 6*np.log2(19/eps_sec) - np.log2(2/eps_cor)
        l = t1 + t2*t3 - t4
        return  l
    else:
        return -1
        

def get_keys(keys_A, keys_B, keys_D, nZ_max):
    """ Get a pair of keys of length nZ_max """
    used_intervals = []
    used_pulses = []
    X_As = []
    Z_As = []
    X_Bs = []
    Z_Bs = []

    j = 0
    while(j < len(keys_A)):

        # Keys for each decoy
        X_A = [[], []]
        Z_A = [[], []]
        X_B = [[], []]
        Z_B = [[], []]
        pulses = 0

        i = 0
        nZ = 0
        while(nZ < nZ_max):
            # Out of range
            if (j+i) >= len(keys_A):
                break

            # Get the partial keys
            nZ_target = nZ_max - nZ
            X_A_p, Z_A_p, X_B_p, Z_B_p, pulses_p = basis_reconciliation(keys_A[j+i], keys_B[j+i], keys_D[j+i], nZ_target)
            pulses += pulses_p

            # Collect the partial keys
            X_A[0] += X_A_p[0]
            X_A[1] += X_A_p[1]
            Z_A[0] += Z_A_p[0]
            Z_A[1] += Z_A_p[1]
            X_B[0] += X_B_p[0]
            X_B[1] += X_B_p[1]
            Z_B[0] += Z_B_p[0]
            Z_B[1] += Z_B_p[1]


            nZ += len(Z_A_p[0])+len(Z_A_p[1])
            i += 1
        
        used_intervals.append((j, j+i))
        used_pulses.append(pulses)
        j += i + 1

        X_As.append(X_A)
        Z_As.append(Z_A)
        X_Bs.append(X_B)
        Z_Bs.append(Z_B)

    return X_As, Z_As, X_Bs, Z_Bs, used_intervals, used_pulses


def get_parameters(X_A, Z_A, X_B, Z_B):
    """Return the useful parameters"""
    # Errors per basis per decoy
    mZ_k_list = errors_estimation(Z_A, Z_B)
    mX_k_list = errors_estimation(X_A, X_B)

    # Number of bits per basis per decoy
    nZ_k_list = [len(Z_k) for Z_k in Z_A]
    nX_k_list = [len(X_k) for X_k in X_A]

    # Decoys
    k_list = [0.6, 0.1818]
    p_k_list = [0.7, 0.3]

    parameters = [k_list, p_k_list, nZ_k_list, mZ_k_list, nX_k_list, mX_k_list]
    #print("Parameters: ", parameters)
    return parameters
    

def do_analyis(keys_A, keys_B, keys_D, nZ_max):
    """Collective function to do the whole analysis"""
    # Get the keys of length nZ_max 
    X_As, Z_As, X_Bs, Z_Bs, used_intervals, used_pulses = get_keys(keys_A, keys_B, keys_D, nZ_max)
    print("Got " + str(len(X_As)) + " keys!")
    
    # Define what to store
    qZ = np.array([])
    qX = np.array([])
    N = np.array([])
    secure_lengths = np.array([])
    used_intervals = np.array(used_intervals)
    used_pulses = np.array(used_pulses)
    
    #Iterate the analysys
    for X_A, Z_A, X_B, Z_B in zip(X_As, Z_As, X_Bs, Z_Bs):
        # Get the parameter
        parameters = get_parameters(X_A, Z_A, X_B, Z_B)
        
        # Compute the security key length
        lambda_EC = 1.16*h2(np.sum(parameters[3])/np.sum(parameters[2]))
        l = secret_key_length(parameters, eps1=1e-9/19, eps2=1e-9/19, eps_sec=1e-9, eps_cor=1e-15, lambda_EC=lambda_EC)
        
        
        qZ = np.append(qZ, np.sum(parameters[3])/np.sum(parameters[2]))
        qX = np.append(qX, np.sum(parameters[5])/np.sum(parameters[4]))
        N = np.append(N, np.sum(parameters[2]) + np.sum(parameters[4]))
        secure_lengths = np.append(l, secure_lengths)
        if l>0: print(l/(np.sum(parameters[2]) + np.sum(parameters[4])))

    #Clean the results
    qZ = qZ[secure_lengths > 0]
    qX = qX[secure_lengths > 0]
    N = N[secure_lengths > 0]
    used_intervals = used_intervals[secure_lengths > 0]
    used_pulses = used_pulses[secure_lengths > 0]
    secure_lengths = secure_lengths[secure_lengths > 0]

    return qZ, qX, N, secure_lengths, used_intervals, used_pulses

