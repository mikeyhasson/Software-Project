import sys  # forum discussion mentioned it's allowed to import sys


def main():
    args = sys.argv
    global k
    global d

    if len(args) > 3 or len(args) < 2:  # only 2 arguments are allowed, k and max_iter
        raise Exception("Number of arguments received is invalid")

    k = int(args[1])
    max_iter = 200

    if len(args) == 3:
        max_iter = int(args[2])

    # if k or max_iter are not in the right format, meaning integers, python will throw exception.

    f = sys.stdin  # input file
    datapoints = f.read().strip().split("\n")
    datapoints = [[float(y) for y in x.split(',')] for x in datapoints]
    d = len(datapoints[0])
    centroids = [datapoints[i][:] for i in range(k)]  # setting first k datapoints as centroids.

    if k >= len(datapoints) or k <= 0:  # checking if k is valid
        raise Exception("k is invalid")

    for i in range(max_iter):
        clusters = [[] for j in range(k)]
        assign_to_clusters(centroids, clusters, datapoints)
        if not update_centroids(centroids, clusters):  # if no change was made, terminate
            break
    s = ""
    for data in centroids:
        for i in range(d - 1):
            s += "{0:.4f}".format(data[i]) + ","
        s += "{0:.4f}".format(data[d - 1]) + "\n"
    s= s[:-1]
    print(s)


def assign_to_clusters(centroids, clusters, datapoints):
    for data in datapoints:
        min_val = None
        min_index = -1
        for i in range(k):
            subtracted_and_squared = [(x - y) ** 2 for x, y in zip(centroids[i], data)]
            cur_val = sum(subtracted_and_squared)
            if min_val is None or min_val > cur_val:
                min_val = cur_val
                min_index = i

        clusters[min_index].append(data)


def update_centroids(centroids, clusters):
    has_changed = False
    for i in range(k):
        cur_cluster = clusters[i]
        cluster_count = len(cur_cluster)
        for j in range(d):
            cur_mean_cluster_line = sum([data[j] for data in cur_cluster]) / cluster_count
            if cur_mean_cluster_line != centroids[i][j]:
                has_changed = True

            centroids[i][j] = cur_mean_cluster_line
    return has_changed


if __name__ == "__main__":
    main()
