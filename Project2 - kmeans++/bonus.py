import matplotlib
from sklearn import datasets
from sklearn import cluster
import matplotlib.pyplot as plt
import math


def find_elbow(interias):
    n_ks = len(interias)
    result = 1
    prev_incline = interias[1] - interias[0]
    min_angle = -1

    for i in range(2, n_ks):
        curr_incline = interias[i] - interias[i-1]
        m1 = curr_incline
        m2 = prev_incline
        curr_angle = math.atan((m1-m2)/(1+m1*m2))
        prev_incline = curr_incline

        if curr_angle < min_angle:
            min_angle = curr_angle
            result = i

    return result


if __name__ == "__main__":
    # Importing the data set:
    data = datasets.load_iris().data

    # Calculating the interias for the desired k values:
    interias = list()

    for k in range(1, 11):
        model = cluster.KMeans(n_clusters=k, random_state=0)
        model.fit(data)
        interias.append(-model.score(data))

    # Finding the elbow:
    k = find_elbow(interias)

    # Plotting the result:
    plt.plot(list(range(1, 11)), interias)
    plt.xlabel("K")
    plt.ylabel("Average Dispersion")
    plt.annotate(text="Elbow Point", xy=(k, interias[k]), xytext=(k+0.01,  interias[k] + 0.01),
                 arrowprops=dict(arrowstyle='<-', connectionstyle='arc3,rad=0.5', color='red'))

    plt.savefig('elbow.png')
