import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
from scipy.stats import binom
import seaborn as sns


def generate_key(length):
    return np.random.randint(0, 2, (1, length))[0]

def add_errors(a, error_prob):
    error_mask = np.random.choice(2, size=a.shape, p=[1.0-error_prob, error_prob])
    return np.where(error_mask, ~a+2, a)

def h2(p):
    return -p*np.log2(p)-(1-p)*np.log2(1-p)

def prob_detected_error(qber, block_length):
    prob = 0
    for i in range(1, block_length+1, 2):
        prob += binom(block_length, qber).pmf(i)
    return prob



def EC_via_syndrome_decoding(x_block1, y_block1, H):
    """Find and correct errors via syndrome decoding"""
    done = False

    #Compute the sydromes
    s_x = np.dot(H, x_block1)%2
    s_y = np.dot(H, y_block1)%2
    
    #Compute the syndrome difference
    S = np.bitwise_xor(s_x, s_y)
    
    if np.sum(S) > 0:
        #Spot the error
        err_index = np.packbits(S, bitorder="little")[0] - 1

        #Correct the error
        y_block1[err_index] = 1 - y_block1[err_index]
        done = True

    return x_block1, y_block1, done



def Winnow(n, k, qber, H):
    """Error detection"""
    #Generate the key and add some errors
    n = int(n)
    x = generate_key(n)
    y = add_errors(x, qber)

    #Length check
    elements = n//(2**k)*(2**k)
    x = x[:elements]
    y = y[:elements]
    errs = np.array([i for i in range(len(x)) if x[i] != y[i]])
    
    #Split the sifted keys into block of length 2**k
    x_blocks = np.array_split(x, n//(2**k))
    y_blocks = np.array_split(y, n//(2**k))

    #Compute the parity of each block
    x_parity = np.array([np.sum(x_i)%2 for x_i in x_blocks])
    y_parity = np.array([np.sum(y_i)%2 for y_i in y_blocks])

    #Compare the parity
    errors_indices = np.array([i for i in range(len(x_parity)) if x_parity[i] != y_parity[i]])

    """Error correction"""
    #Perform error correction
    for error_index in errors_indices:

        #Discard the last bit (privacy maintencance)
        x_block1 = x_blocks[error_index][:-1].copy()
        y_block1 = y_blocks[error_index][:-1].copy()

        #Syndrome_decoding
        x_block1, y_block1, done = EC_via_syndrome_decoding(x_block1, y_block1, H)

        #Remerge the last bit to the block
        if done:
            #If EC then the error has already been found
            y_blocks[error_index][:-1] = y_block1
        else:
            #If not, the error was in the last bit
            y_blocks[error_index][-1] = 1 - y_blocks[error_index][-1]

    #Keys
    x = np.concatenate(x_blocks)
    y = np.concatenate(y_blocks)
    errs = np.array([i for i in range(len(x)) if x[i] != y[i]])
    
    #Compute how many bits I revealed to Eve
    disclosed_bits = len(x_parity)
    disclosed_bits += len(errors_indices)*k

    return disclosed_bits, len(errs)


def construct_H(k):
    n = int(2**k-1)
    H = 3*np.ones(shape=(k, n))
    for i in range(k):
        for j in range(n):
            H[i, j] = np.floor((j+1)/2**i)%2
    return np.int64(H)
        


def binary_search(x_blocks, y_blocks, disclosed_bits):
    """Binary search: find and correct errors via iterative parity comparison"""
    #Compute the parity of each block
    x_parity = np.array([np.sum(x_i)%2 for x_i in x_blocks])
    y_parity = np.array([np.sum(y_i)%2 for y_i in y_blocks])
    
    #Compare the parity and find the blocks with errors
    errors_indices = np.array([i for i in range(len(x_parity)) if x_parity[i] != y_parity[i]])
    
    # Base case: If there are no errors, return
    if len(errors_indices) == 0:
        return x_blocks, y_blocks, disclosed_bits
    
    #Iterate over the blocks with errors and split them into two subblocks
    for error_index in errors_indices:
        x_block = x_blocks[error_index]
        y_block = y_blocks[error_index]
        
        #Check if the error has been found now
        if len(x_block) == 1:
            #Correct the error
            y_block = 1 - y_block
            y_blocks[error_index] = y_block
            
        else: 
            #Split the block with an error
            x_subblocks = np.array_split(x_block, 2)
            y_subblocks = np.array_split(y_block, 2)
            
            #Update the number of disclosed bits
            disclosed_bits += 1

            # Recursively call the function on the subblocks
            corrected_x_subblocks, corrected_y_subblocks, disclosed_bits = binary_search(x_subblocks, y_subblocks, disclosed_bits)

            # Update the corrected subblocks in the original x_blocks and y_blocks
            x_blocks[error_index] = np.concatenate(corrected_x_subblocks)
            y_blocks[error_index] = np.concatenate(corrected_y_subblocks)
    
    return x_blocks, y_blocks, disclosed_bits
    

def Cascade(n, qber, k1=None, it=2):
    #Compute k1, initialize disclosed_bits
    if k1==None: k1 = int(np.ceil(0.73/qber))
    n = int(n)
    disclosed_bits = 0

    #Generate the key and add some errors
    x = generate_key(n)
    y = add_errors(x, qber)

    # Define the single iteration
    def one_iteration(x, y, k1, disclosed_bits):
        """First step"""
        #Split the sifted keys into block of length (approx) k1
        x_blocks = np.array_split(x, n//k1)
        y_blocks = np.array_split(y, n//k1)

        #Perform binary search and error correction
        disclosed_bits += len(y_blocks)
        result = binary_search(x_blocks, y_blocks, disclosed_bits)
        disclosed_bits = result[2]
        
        #Update parameter for further steps
        x = np.concatenate(result[0])
        y = np.concatenate(result[1])
        k2 = 2*k1    

        """Second step"""
        #Random permutation of the key
        shuffle_indices = np.random.permutation(np.arange(n))
        x_s = x[shuffle_indices]
        y_s = y[shuffle_indices]

        #Split the new keys into block of length (approx) k2
        x_blocks1 = np.array_split(x_s, n//k2)
        y_blocks1 = np.array_split(y_s, n//k2)

        #Perform (again) binary search and error correction
        disclosed_bits += len(y_blocks1)
        result = binary_search(x_blocks1, y_blocks1, disclosed_bits)
        disclosed_bits = result[2]
        x_s = np.concatenate(result[0])
        y_s = np.concatenate(result[1])

        #Inverse transformation
        inverse_indices = np.argsort(shuffle_indices)
        x1 = x_s[inverse_indices]
        y1 = y_s[inverse_indices]

        #Find the errors corrected in the second step
        errs = [i for i in range(n) if y[i] != y1[i]]

        #Re-split again with k1
        x_blocks = np.array_split(x1, n//k1)
        y_blocks = np.array_split(y1, n//k1)

        for err in errs:
            #Find the correspondent blocks of the first step
            cumulative_lengths = np.cumsum([len(b) for b in x_blocks])
            block_index = np.searchsorted(cumulative_lengths, err, side='right')

            #Perform binary search only in that block
            x_block = x_blocks[block_index]
            y_block = y_blocks[block_index]
            disclosed_bits += 1
            result = binary_search([x_block], [y_block], disclosed_bits)
            disclosed_bits = result[2]
            
            x_blocks[block_index] = result[0][0]
            y_blocks[block_index] = result[1][0]

        x = np.concatenate(x_blocks)
        y = np.concatenate(y_blocks)

        return x, y, disclosed_bits
    
    for _ in range(it):
        x, y, disclosed_bits = one_iteration(x, y, k1, disclosed_bits)
        k1 = 2*k1
    
    errs = len([i for i in range(len(x)) if x[i] != y[i]])

    return disclosed_bits, errs