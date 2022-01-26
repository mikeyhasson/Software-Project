import mykmeanssp as Kmeans
import numpy as np;
import pandas as pd;
import sys;


def kmeans_pp(data, n, k, d):
    np.random.seed(0);

    centroids_indices = np.empty(shape=(k, 1), dtype=int);
    centroids_indices[0][0] = np.random.randint(0, n-1);
    centroids = np.ndarray(shape=(k, d));
    centroids[0] = data[centroids_indices[0]];
    indexes_array = np.arange(0, n);
    probs = np.ones(shape=(n, )) * sys.float_info.max

    z = 1;
    while z < k:
        for i in range(0, n):
            calc = sum((data[i] - centroids[z-1]) ** 2);
            if calc <= probs[i]:
                probs[i] = calc;

        centroids_indices[z][0] = np.random.choice(indexes_array, 1, p=probs/sum(probs))[0];
        centroids[z] = data[centroids_indices[z][0]];

        z += 1;

    return centroids, centroids_indices; # the initial centroids for the kmeans algorithm;


def scan_files(fd1, fd2):
    ''' returns the collected input out of the given files '''

    df1 = pd.read_csv(fd1, header=None);
    df2 = pd.read_csv(fd2, header=None);
    # files are supplied in csv format;

    d1 = df1.shape[1];
    d2 = df2.shape[1];
    df1.columns = np.arange(0, d1);
    df2.columns = np.arange(0, d2);
    df1 = df1.sort_values(0);
    df2 = df2.sort_values(0);
    df = pd.merge(df1, df2, how='inner', on=0);  # inner join on the first column which is the index of every row;
    df = df.sort_values(0);  # sorting the data points based on the first column;
    df.drop(columns=0, axis=1, inplace=True);
    
    return df.to_numpy(); # want to work with this type of object


def pack_input(n, k, d, num_of_iter, data_points, initial_centroids):
    ''' Packs the parameters to a valid input list of the C API function and returns it.
        n, k, d, num_of_iter - are integers
        data_points - the given data points - a pandas data frame
        initial_centroids - the calculated initial centroid - a panda data frame
    '''

    list = [0.1] * (4 + n*d + k*d);  # a list of size 4+(n+k)*d

    list[0] = n;
    list[1] = k;
    list[2] = d;
    list[3] = num_of_iter;

    for i in range(0, n):
        for j in range(0, d):
            list[4 + i*d + j] = data_points[i, j];

    for i in range(0, k):
        for j in range(0, d):
            list[4 + n*d + i*d + j] = initial_centroids[i, j];

    return list;


def unpack_output(output, k, d):
    result = np.empty(shape=(k, d));

    for i in range(0, k):
        for j in range(0, d):
            result[i][j] = output[d*i + j];

    return result;


def print_matrix(matrix, end='\n'):
    m = matrix.shape[0];
    n = matrix.shape[1];

    for i in range(0, m-1):
        for j in range(0, n-1):
            print(matrix[i][j], end=',');  # without '\n'
        print(matrix[i][n-1]); # with '\n'

    for j in range(0, n-1):
        print(matrix[m-1][j], end=',');  # without '\n'
    print(matrix[m-1][n-1], end=end); # with the desired ender;


def output_result(centroids, indices):
    print_matrix(indices);  # a (k, 1) np.ndarray
    print_matrix(centroids);  # a (k, d) np.ndarray


if __name__ == "__main__":

    args = sys.argv;
    num_of_args = len(args) - 1;

    k = int(args[1]);

    if num_of_args == 4:
        max_iter = int(args[2]);
        file_name_1 = args[3];
        file_name_2 = args[4];
    elif num_of_args == 3:
        max_iter = 300;
        file_name_1 = args[2];
        file_name_2 = args[3];
    else:
        print("Invalid number of arguments!");
        sys.exit();

    f1 = open(file_name_1, 'r');
    f2 = open(file_name_2, 'r');

    data_points = scan_files(f1, f2);

    n = data_points.shape[0];  # the number of data points;
    d = data_points.shape[1];  # the dimension of every data point;

    if (k <= 0) or (k >= n) or (max_iter < 0):
        print("Invalid Arguments!")
        sys.exit()

    f1.close();
    f2.close();

    initial_centroids, initial_centroids_indices = kmeans_pp(data_points, n, k, d);
    initial_centroids_indices = initial_centroids_indices.reshape(1, -1);

    c_api_func_input = pack_input(n, k, d, max_iter,data_points, initial_centroids);

    final_centroids = unpack_output(Kmeans.kmeans(c_api_func_input), k, d);

    final_centroids = np.round(final_centroids, decimals=4);

    output_result(final_centroids, initial_centroids_indices);






