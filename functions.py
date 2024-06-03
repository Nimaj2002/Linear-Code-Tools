import numpy as np
from itertools import product, combinations


def mod_inv(a, p):
    """ Returns the modular multiplicative inverse of a under modulo p. """
    for x in range(1, p):
        if (a * x) % p == 1:
            return x
    raise ValueError(f"No modular inverse for {a} under modulo {p}")


def transform_to_standard_form(G, base):
    G = np.array(G, dtype=int)
    k, n = G.shape

    # Gaussian elimination to form identity matrix I_k on the left
    for i in range(k):
        # Find pivot
        for j in range(i, n):
            if G[i][j] % base != 0:
                break
        else:
            raise ValueError("Matrix cannot be transformed to standard form")

        # Swap columns if necessary
        if j != i:
            G[:, [i, j]] = G[:, [j, i]]

        # Normalize pivot row
        pivot_inv = mod_inv(G[i][i] % base, base)
        G[i] = (G[i] * pivot_inv) % base

        # Eliminate other entries in pivot column
        for r in range(k):
            if r != i:
                G[r] = (G[r] - G[r][i] * G[i]) % base

    # Swap columns to get the identity matrix in the leftmost position
    for i in range(k):
        if G[i][i] % base != 1:
            for j in range(i + 1, n):
                if G[i][j] % base == 1:
                    G[:, [i, j]] = G[:, [j, i]]
                    break

    return G


def parity_check_matrix_generator(G, b):
    k, n = G.shape
    submatrix = G[0:k, k:n]
    submatrix = submatrix.T
    submatrix = submatrix * -1
    submatrix = submatrix + b
    submatrix = submatrix % b

    I = np.eye(n-k)
    H = np.concatenate((submatrix, I), axis=1)
    H = H.astype(int)
    return H


def generate_codewords(G, base):
    k, n = G.shape
    # Generate all possible input vectors (2^k combinations)
    input_vectors = list(product(range(base), repeat=k))

    codewords = []
    for input_vector in input_vectors:
        input_vector = np.array(input_vector)
        codeword = np.dot(input_vector, G) % base
        codewords.append(codeword)

    return np.array(codewords)


def coset_leaders(generator_matrix, base):
    k, n = generator_matrix.shape
    identity_matrix = np.eye(n - k, dtype=int)
    coset_leaders = []

    for coset_representative in identity_matrix:
        coset_leader = coset_representative @ generator_matrix
        coset_leader = coset_leader + base
        coset_leader = coset_leader % base
        coset_leader = coset_leader.astype(int)
        # coset_leader = np.mod(np.dot(coset_representative, generator_matrix.T), 2)
        coset_leaders.append(coset_leader)

    return coset_leaders


def calculate_syndrome(vector, H, base):
    return np.dot(vector, H.T) % base


def find_coset_leaders(G, H, base):
    k, n = G.shape

    all_vectors = list(product(range(base), repeat=n))
    syndromes = {}

    for vector in all_vectors:
        vector = np.array(vector)
        syndrome = tuple(calculate_syndrome(vector, H, base))
        weight = np.sum(vector != 0)

        if syndrome not in syndromes or weight < syndromes[syndrome][1]:
            syndromes[syndrome] = (vector, weight)

    coset_leaders = {k: v[0] for k, v in syndromes.items()}

    return coset_leaders


def generate_syndrome_table(H, base):
    # Number of bits in the code
    n = H.shape[1]

    # Generate all possible error vectors
    error_vectors = [
        np.array(list(format(i, f'0{n}b')), dtype=int) for i in range(base**n)]
    # Calculate syndromes for each error vector
    syndrome_dict = {}
    for e in error_vectors:
        e = np.array(e)
        syndrome = tuple(e @ H.T % base)
        if syndrome not in syndrome_dict:
            syndrome_dict[syndrome] = [e]
        else:
            syndrome_dict[syndrome].append(e)

    return syndrome_dict


def arrays_with_lowest_sum(array_of_arrays):
    if not array_of_arrays:
        return []

    # Calculate the sum of each sub-array
    sums = [sum(arr) for arr in array_of_arrays]

    # Find the minimum sum
    min_sum = min(sums)

    # Collect all sub-arrays that have the minimum sum
    result = [list(arr)
              for arr, s in zip(array_of_arrays, sums) if s == min_sum]

    return result


def hamming_distance(x, y):
    """Calculate the Hamming distance between two binary vectors."""
    return np.sum(x != y)


def calculate_minimum_distance(codewords):
    """Calculate the minimum distance of the code defined by the generator matrix G."""
    min_distance = np.inf
    for codeword1, codeword2 in combinations(codewords, 2):
        distance = hamming_distance(codeword1, codeword2)
        if distance < min_distance:
            min_distance = distance
    return min_distance
