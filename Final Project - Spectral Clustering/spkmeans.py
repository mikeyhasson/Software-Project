
import capi_spk as spk
import numpy as np
import sys
from sklearn.preprocessing import normalize

def kmeans_pp(data, k):
    n, d = data.shape
    centroids_indices = np.empty(shape=(k,), dtype=int)
    centroids_indices[0] = np.random.randint(0, n-1)
    centroids = np.ndarray(shape=(k, d))
    centroids[0] = data[centroids_indices[0]]

    z = 1
    while z < k:
        probs = np.ndarray(shape=(n, ))

        for i in range(0, n):
            D_i = sum((data[i]-centroids[0]) ** 2)
            for j in range(1, z):
                calc = sum((data[i] - centroids[j]) ** 2)
                if calc < D_i:
                    D_i = calc

            probs[i] = D_i

        probs = probs / sum(probs)
        centroids_indices[z] = np.random.choice(n, 1, p=probs)
        centroids[z] = data[centroids_indices[z]]

        z += 1

    return centroids, centroids_indices # the initial centroids for the kmeans algorithm


def file_to_Array(fd):
    almost_vectors = fd.readlines()
    n = len(almost_vectors)
    d = len(almost_vectors[0].split(","))
    data_points = np.ndarray(shape=(n, d))

    # converting the rows into vectors
    for i in range(0, n):
        stringed_row = almost_vectors[i].split(',')
        data_points[i] = [float(val) for val in stringed_row]

    return data_points


def print_double(x, sep):
    if (x < 0) and (x > -0.00005):
        print("0.0000", end=sep)
    else:
        print("{:.4f}".format(x), end=sep)


def print_matrix(matrix):
    m = matrix.shape[0];
    n = matrix.shape[1];

    for i in range(0, m-1):
        for j in range(0, n-1):
            print_double(matrix[i][j], sep=',')
        
        print_double(matrix[i][n-1], sep='\n')

    for j in range(0, n-1):
        print_double(matrix[m-1][j], sep=',')

    print_double(matrix[m-1, n-1], sep='')


def print_Vector(vec, is_int=False):
    n = vec.shape[0]

    if is_int is False:
        for i in range(0, n-1):
            print_double(vec[i], sep=",")

        print_double(vec[n-1], sep="\n")
    else:
        for i in range(0, n-1):
            print(vec[i], end=',')

        print(vec[n-1], end="\n")


def mat_to_list(mat):
    n, d = mat.shape
    list = [1.1] * (n*d)

    for i in range(0, n):
        for j in range(0, d):
            list[i*d + j] = mat[i, j]

    return list

def vector_to_list(vec):
    n = vec.shape[0]
    list = [1.1] * n

    for i in range(0, n):
        list[i] = vec[i]

    return list


def list_to_matrix(list, n, m):

    mat = np.ndarray(shape=(n, m))

    for i in range(0, n):
        for j in range(0, m):
            mat[i, j] = list[i * m + j]

    return mat


def list_to_vector(list, size=None):

    if size is None:
        size = len(list)

    vec = np.ndarray(shape=(size, ))

    for i in range(0, size):
        vec[i] = list[i]

    return vec


def print_input_err_exit():
    print("Invalid Input!", end='')
    exit(1)


def print_err_exit():
    print("An Error Has Occured", end='')
    exit(1)


if __name__ == "__main__":  # The Python Interface Implementation
    np.random.seed(0)  #IMPORTANT

    argv = sys.argv
    argc = len(argv)

    if argc != 4:
        print_input_err_exit()

    goal = argv[2]
    filename = argv[3]

    fd = open(filename, 'r')
    data = file_to_Array(fd)
    fd.close()
    data_list = mat_to_list(data)
    n, d = data.shape

    if goal == "wam":
        wam = spk.wam(data_list, n, d)
        wam = list_to_matrix(wam, n, n)  # converts a 2d list to a np.ndarray
        print_matrix(wam)

        exit(0)

    if goal == "ddg":
        ddg = spk.ddg(data_list, n, d)
        ddg = list_to_matrix(ddg, n, n)  # converts a 2d array to an np.ndarray
        print_matrix(ddg)

        exit(0)

    if goal == "lnorm":
        lnorm = spk.lnorm(data_list, n, d)
        lnorm = list_to_matrix(lnorm, n, n)  # converts a 2d array to an np.ndarray
        print_matrix(lnorm)

        exit(0)

    if goal == "jacobi":

        eig_Values_list, eig_Vectors_list = spk.jacobi(data_list, n)

        eig_Values = list_to_vector(eig_Values_list)
        eig_Vectors = list_to_matrix(eig_Vectors_list, n, n)

        print_Vector(eig_Values)
        print_matrix(eig_Vectors.T)

        exit(0)

    if goal == "spk":

        try:
            k = int(argv[1])
        except ValueError:
            print_input_err_exit()

        if (k < 0) or (k > n):
            print_input_err_exit()

        lnorm = spk.lnorm(data_list, n, d)

        eig_Values_list, eig_Vectors_list = spk.sorted_jacobi(lnorm, n)

        eig_Vectors = list_to_matrix(eig_Vectors_list, n, n)
        eig_Values = list_to_vector(eig_Values_list, n)
        
        if k == 0:
            # find the suitable k
            k = np.argmax(np.diff(eig_Values[:n//2 + 1])) + 1  # np.argmax takes the minimum at argmax as wanted.
        
        U = eig_Vectors[:, :k]
        T = normalize(U, axis=1, norm='l2')

        init_centroids, init_centroids_indices = kmeans_pp(T, k)

        final_centroids_list = spk.kmeans(mat_to_list(T), 
            vector_to_list(init_centroids_indices), n, k, k)
        
        final_centroids = list_to_matrix(final_centroids_list, k, k)

        #output the initial centroids indices, and then the final centroids calculated.
        print_Vector(init_centroids_indices, is_int=True)
        print_matrix(final_centroids)

        exit(0)

    # else:
    print_input_err_exit()
