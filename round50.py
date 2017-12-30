import matplotlib.pyplot as plt
import numpy as np
import math
from random import *

N = 4000 #  global population
t = 4 # time spent in group
T = 1000 # number of generation
K = 0.1 # death rate
R_small = 4 # total resource, R for small group of 4
R_large = 50 # R for large group of 40
Gc = 0.018 # cooperative growth rate
Gs = 0.02 # selfish growth rate
Cc = 0.1 # cooperative consumption rate
Cs = 0.2 # selfish consumption rate

# lplot = []
geno0, geno1, geno2, geno3 = [], [], [], []

def genoResource(ni, nj, Gi, Ci, Gj, Cj, R): # group resource the genotype received
    ri = (ni*Gi*Ci*R) / ((nj*Gj*Cj) + (ni*Gi*Ci))
    return ri

def genoNumber(ni, ri, Ci, K): # number of individuals in a group with the genotype
    n = ni + (ri/Ci) - (K*ni)
    return n

#1# Initialise the migrant pool with N individuals.
# 4 possible genotypes
genoSize = N/4
genotypes = [genoSize, genoSize, genoSize, genoSize] # total = N
# 0 cooperative + small
# 1 cooperative + large
# 2 selfish + small
# 3 selfish + large

# for y in range(1):
for y in range(T):

#2# Group formation (aggregation): Assign individuals in the migrant pool to groups
# small size groups, size 4
    print "generation " + str(y+1)
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
        # print smallGroup
        smallGroups = np.append(smallGroups, smallGroup) # append into all small groups

    smallGroups = np.reshape(smallGroups, (-1, 4))
    if (j != 3):
        smallGroups = smallGroups[:-1]
    # print smallGroups

    # large size groups, size 40
    largeGroups = [] # all large grops
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
        largeGroups = np.append(largeGroups, largeGroup) # append into all small groups

    largeGroups = np.reshape(largeGroups, (-1, 40))
    if (j != 39):
        largeGroups = largeGroups[:-1]
    # print largeGroups

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

        for x in range(t):
            # compute resource consumed
            g0t = g0
            g2t = g2
            r0 = genoResource(g0t, g2t, Gc, Cc, Gs, Cs, R_small)
            r2 = R_small - r0
            # print "resource consumed by g0 " + str(r0)
            # print "resource consumed by g2 " + str(r2)

            # reproduce genotypes
            g0 = int(genoNumber(g0t, r0, Cc, K)) # g0
            # print "number of g0 " + str(g0)
            g2 = int(genoNumber(g2t, r2, Cs, K)) # g2
            # print "number of g2 " + str(g2)

        #4# Migrant pool formation (dispersal): Return the progeny of each group to the migrant pool.
        genotypes[0] = genotypes[0] + g0
        genotypes[2] = genotypes[2] + g2

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

        for x in range(t):
            g1t = g1
            g3t = g3
            # compute resource consumed
            r1 = genoResource(g1t, g3t, Gc, Cc, Gs, Cs, R_large)
            r3 = R_large - r1
            # print "resource consumed by g1 " + str(r1)
            # print "resource consumed by g3 " + str(r3)

            # reproduce genotypes
            g1 = int(genoNumber(g1t, r1, Cc, K)) # g1
            # print "number of g1 " + str(g1)
            g3 = int(genoNumber(g3t, r3, Cs, K)) # g3
            # print "number of g3 " + str(g3)

    #4# Migrant pool formation (dispersal)
        genotypes[1] = genotypes[1] + g1
        genotypes[3] = genotypes[3] + g3
    print genotypes

    #5# Maintaining the global carrying capacity: Rescale the migrant pool back to size N, retaining the proportion of individuals with each genotype.
    groupSum = np.sum(genotypes)
    # print groupSum
    d = float(groupSum)/N
    d = math.ceil(d)
    d = int(d)
    # print d
    genotypes = np.array(genotypes, dtype=float)
    genotypes = genotypes/d
    genotypes = np.array(genotypes, dtype=int)
    print genotypes
    # largeNum = genotypes[1] + genotypes[3]
    # lplot.append(largeNum)
    # selfNum = genotypes[2] + genotypes[3]
    # splot.append(selfNum)
    geno0.append(genotypes[0])
    geno1.append(genotypes[1])
    geno2.append(genotypes[2])
    geno3.append(genotypes[3])
    print np.sum(genotypes)

#6# Iteration: Repeat from step 2 onwards for a number of generations, T

# plt.plot(lplot)
plt.plot(geno0)
plt.plot(geno1)
plt.plot(geno2)
plt.plot(geno3)
plt.xlabel('Generation')
plt.ylabel('global genotype')
plt.savefig('plot.png', bbox_inches='tight')
