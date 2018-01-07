import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.legend_handler import HandlerLine2D
from random import *

N = 4000 #  global population
t = 4 # time spent in group
T = 1000 # number of generation
K = 0.1 # death rate
R_small = 4.0 # total resource, R for small group of 4
R_large = 50.0 # R for large group of 40
Gc = 0.018 # cooperative growth rate
Gs = 0.02 # selfish growth rate
Cc = 0.1 # cooperative consumption rate
Cs = 0.2 # selfish consumption rate

# lplot = []
geno0, geno1, geno2, geno3 = [], [], [], []

def consume(x):
    y = (23.0 * x - 20)/18
    # y = ((23.0*x*x)-(184*x)+2960)/648
    return y

def genoResource(ni, nj, Gi, Ci, Gj, Cj, R): # group resource the genotype received
    ri = (float(ni)*Gi*Ci*R) / ((float(nj)*Gj*Cj) + (float(ni)*Gi*Ci))
    return ri

def genoNumber(ni, ri, Ci, K): # number of individuals in a group with the genotype
    n = float(ni) + (ri/Ci) - (K*float(ni))
    return n

#1# Initialise the migrant pool with N individuals.
# 4 possible genotypes
genoSize = N/4.0
genotypes = [genoSize, genoSize, genoSize, genoSize] # total = N
# 0 cooperative + small
# 1 cooperative + large
# 2 selfish + small
# 3 selfish + large

genoPlot = [0, 0, 0, 0] # for plotting
geno0.append(0.25)
geno1.append(0.25)
geno2.append(0.25)
geno3.append(0.25)

# for y in range(1):
# for y in range(T):
for y in range(30):

#2# Group formation (aggregation): Assign individuals in the migrant pool to groups
# small size groups, size 4
    print "generation " + str(y+1)
    print "genotypes= " + str(genotypes)
    smallGroups = [] # previous all small groups
    smallGroup = []
    for _ in range(4):
        smallGroup.append(0)
    cont = 1
    while (cont == 1) and (genotypes[0] > 0) and (genotypes[2] > 0):
        for j in range(4):
            flip = random()
            if (flip < 0.5):
                if (genotypes[0] > 0):
                    # print "coop"
                    smallGroup[j] = 0 # put coop into group
                    genotypes[0] = genotypes[0] - 1
                    # print "genotypes[0]= " + str(genotypes[0])
                else:
                    cont = 0
                    break
            else:
                if (genotypes[2] > 0):
                    # print "selfish"
                    smallGroup[j] = 2 # put selfish into group
                    genotypes[2] = genotypes[2] - 1
                    # print "genotypes[2]= " + str(genotypes[2])
                else:
                    cont = 0
                    break
        # print "genotypes " + str(genotypes)
        # print smallGroup
        smallGroups = np.append(smallGroups, smallGroup) # append into all small groups

    smallGroups = np.reshape(smallGroups, (-1, 4))
    if (j != 3):
        smallGroups = smallGroups[:-1]
    # print smallGroups
    print "genotypes " + str(genotypes)

    # large size groups, size 40
    largeGroups = [] # all large groups
    largeGroup = []
    for _ in range(40):
        largeGroup.append(0)
    cont = 1
    while (cont == 1) and (genotypes[1] > 0) and (genotypes[3] > 0):
        for j in range(40):
            flip = random()
            if (flip < 0.5):
                if (genotypes[1] > 0):
                    # print "coop"
                    largeGroup[j] = 1 # put coop into group
                    genotypes[1] = genotypes[1] - 1
                    # print "genotypes[1]= " + str(genotypes[1])
                else:
                    cont = 0
                    break
            else:
                if (genotypes[3] > 0):
                    # print "selfish"
                    largeGroup[j] = 3 # put selfish into group
                    genotypes[3] = genotypes[3] - 1
                    # print "genotypes[3]= " + str(genotypes[3])
                else:
                    cont = 0
                    break
        # print largeGroup
        largeGroups = np.append(largeGroups, largeGroup) # append into all large groups

    largeGroups = np.reshape(largeGroups, (-1, 40))
    if (j != 39):
        largeGroups = largeGroups[:-1]
    # print largeGroups
    # print "genotypes " + str(genotypes)

    #3# Reproduction: Perform reproduction within groups for t time-steps
    genotypes = [0, 0, 0, 0]
    # reproduce small groups
    for j in range(smallGroups.shape[0]):
        g0 = 0 # number of genotype 0
        g2 = 0 # number of genotype 2
        for i in range(4):
            if smallGroups[j][i] == 0:
                g0 += 1
            else:
                g2 += 1
        # print "g0 " + str(g0)
        # print "g2 " + str(g2)
        # r0 = float(g0)/(g0+g2)
        # r2 = float(g2)/(g0+g2)
        # print "ratio of g0 " + str(r0)
        # print "ratio of g2 " + str(r2)
        totalI = R_small
        for x in range(t):
            # compute resource consumed
            print "size " + str(totalI)
            Rn = consume(totalI)
            print "total resource " + str(Rn)
            r0 = genoResource(g0, g2, Gc, Cc, Gs, Cs, Rn)
            r2 = Rn - r0
            # print "resource consumed by g0 " + str(r0)
            # print "resource consumed by g2 " + str(r2)

            # reproduce genotypes
            g0 = genoNumber(g0, r0, Cc, K) # g0
            g2 = genoNumber(g2, r2, Cs, K) # g2
            # print "number of g0 " + str(g0)
            # print "number of g2 " + str(g2)
            # r0 = float(g0)/(g0+g2)
            # r2 = float(g2)/(g0+g2)
            # print "ratio of g0 " + str(r0)
            # print "ratio of g2 " + str(r2)
            totalI = g0 + g2
        # g0 = g0t - g0
        # g2 = g2t - g2
        #4# Migrant pool formation (dispersal): Return the progeny of each group to the migrant pool.
        genotypes[0] = genotypes[0] + g0
        genotypes[2] = genotypes[2] + g2
    print "genotypes " + str(genotypes)

    # reproduce large groups
    for j in range(largeGroups.shape[0]):
        g1 = 0 # number of genotype 0
        g3 = 0 # number of genotype 2
        for i in range(40):
            if largeGroups[j][i] == 1:
                g1 += 1
            else:
                g3 += 1
        # print "g1 " + str(g1)
        # print "g3 " + str(g3)
        # r1 = float(g1)/(g1+g3)
        # r3 = float(g3)/(g1+g3)
        # print "ratio of g1 " + str(r1)
        # print "ratio of g3 " + str(r3)
        totalI = R_large
        for x in range(t):
            # compute resource consumed
            # print "size " + str(totalI)
            Rn = consume(totalI)
            # print "total resource " + str(Rn)
            r1 = genoResource(g1, g3, Gc, Cc, Gs, Cs, Rn)
            r3 = Rn - r1
            # print "resource consumed by g1 " + str(r1)
            # print "resource consumed by g3 " + str(r3)

            # reproduce genotypes
            g1 = (genoNumber(g1, r1, Cc, K)) # g1
            g3 = (genoNumber(g3, r3, Cs, K)) # g3
            # print "number of g1 " + str(g1)
            # print "number of g3 " + str(g3)
            # r1 = float(g1)/(g1+g3)
            # r3 = float(g3)/(g1+g3)
            # print "ratio of g1 " + str(r1)
            # print "ratio of g3 " + str(r3)
            totalI = g1 + g3
        # g1 = g1t - g1
        # g3 = g3t - g3
    #4# Migrant pool formation (dispersal)
        genotypes[1] = genotypes[1] + g1
        genotypes[3] = genotypes[3] + g3

    # print "genotypes " + str(genotypes)
    #5# Maintaining the global carrying capacity: Rescale the migrant pool back to size N, retaining the proportion of individuals with each genotype.
    groupSum = np.sum(genotypes)
    # print groupSum
    d = (groupSum)/N
    # d = math.ceil(d)
    # d = int(d)
    # print d
    genotypes = genotypes/d
    genoSum = np.sum(genotypes)
    # print "genoSum= " + str(genoSum)
    genoPlot[0] = genotypes[0]/genoSum
    genoPlot[1] = genotypes[1]/genoSum
    genoPlot[2] = genotypes[2]/genoSum
    genoPlot[3] = genotypes[3]/genoSum
    # print "genotypes= " + str(genotypes)
    genotypes = np.array(genotypes, dtype=int)
    genotypes = np.array(genotypes, dtype=float)

    # largeNum = genotypes[1] + genotypes[3]
    # lplot.append(largeNum)
    # selfNum = genotypes[2] + genotypes[3]
    # splot.append(selfNum)
    # print genoPlot
    geno0.append(genoPlot[0])
    geno1.append(genoPlot[1])
    geno2.append(genoPlot[2])
    geno3.append(genoPlot[3])

#6# Iteration: Repeat from step 2 onwards for a number of generations, T

plt.plot(geno0, label='cooperative + small')
plt.plot(geno1, label='cooperative + large')
plt.plot(geno2, label='selfish + small')
plt.plot(geno3, label='selfish + large')
plt.legend(loc='lower right', shadow=True, fontsize='medium')
plt.xlabel('Generation')
plt.ylabel('Global genotype frequency')
plt.savefig('plot.png', bbox_inches='tight')
plt.show()
