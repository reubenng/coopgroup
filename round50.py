import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.legend_handler import HandlerLine2D
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

def consume(x):
    y = (23.0 * x - 20)/18
    # y = ((23.0*x*x)-(184*x)+2960)/648
    return y

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

geno0, geno1, geno2, geno3 = [], [], [], []
genoPlot = [0, 0, 0, 0] # for plotting
geno0.append(0.25)
geno1.append(0.25)
geno2.append(0.25)
geno3.append(0.25)
type0, type1 = [], []

# for y in range(T):
for y in range(120):

#2# Group formation (aggregation): Assign individuals in the migrant pool to groups
# small size groups, size 4
    print "generation " + str(y+1)
    print "genotypes= " + str(genotypes)

    pool = genotypes
    pool = np.array(pool, dtype=float)

    smallGroups = [] # previous all small groups
    smallGroup = []
    for _ in range(4):
        smallGroup.append(0)

    cont = 1
    totalG = genotypes[0] + genotypes[2]
    if (totalG != 0):
        g0p = float(genotypes[0])/totalG
        g2p = float(genotypes[2])/totalG

        while (cont == 1) and (genotypes[0] > 0) and (genotypes[2] > 0):
            for j in range(4):
                totalG = genotypes[0] + genotypes[2]
                # print "totalG " + str(totalG)
                if (totalG != 0):
                    g0p = float(genotypes[0])/totalG
                    g2p = float(genotypes[2])/totalG
                    # print "genotypes[0] " + str(genotypes[0])
                    # print "genotypes[2] " + str(genotypes[2])
                    # print "g0p " + str(g0p)
                    # print "g2p " + str(g2p)
                    flip = random()
                    if (flip <= g0p):
                        if (genotypes[0] > 0):
                            smallGroup[j] = 0 # put coop into group
                            genotypes[0] = genotypes[0] - 1
                        else:
                            cont = 0
                            break
                    else:
                        if (genotypes[2] > 0):
                            smallGroup[j] = 2 # put selfish into group
                            genotypes[2] = genotypes[2] - 1
                        else:
                            cont = 0
                            break
            smallGroups = np.append(smallGroups, smallGroup) # append into all small groups
            skip = 0
        else:
            skip = 1

    smallGroups = np.reshape(smallGroups, (-1, 4))
    if (j != 3):
        smallGroups = smallGroups[:-1]
    # print smallGroups

    # large size groups, size 40
    largeGroups = [] # all large groups
    largeGroup = []
    for _ in range(40):
        largeGroup.append(0)

    cont = 1
    totalG = genotypes[1] + genotypes[3]
    if (totalG != 0):
        g1p = float(genotypes[1])/totalG
        g3p = float(genotypes[3])/totalG
        while (cont == 1) and (genotypes[1] > 0) and (genotypes[3] > 0):
            for j in range(40):
                totalG = genotypes[1] + genotypes[3]
                if (totalG != 0):
                    g1p = float(genotypes[1])/totalG
                    flip = random()
                    if (flip <= g1p):
                        if (genotypes[1] > 0):
                            largeGroup[j] = 1 # put coop into group
                            genotypes[1] = genotypes[1] - 1
                        else:
                            cont = 0
                            break
                    else:
                        if (genotypes[3] > 0):
                            largeGroup[j] = 3 # put selfish into group
                            genotypes[3] = genotypes[3] - 1
                        else:
                            cont = 0
                            break
            largeGroups = np.append(largeGroups, largeGroup) # append into all large groups
            skip = 0
        else:
            skip = 1

    largeGroups = np.reshape(largeGroups, (-1, 40))
    if (j != 39):
        largeGroups = largeGroups[:-1]
    # print largeGroups

    #3# Reproduction: Perform reproduction within groups for t time-steps
    if (skip == 0):
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
        totalI = g0 + g2
        for x in range(t):
            # compute resource consumed
            # Rn = consume(totalI)
            print "selfish " + str(pool[2]/(pool[0]+pool[2]))
            if (pool[2]/(pool[0]+pool[2]) >= 0.5): # if cheaters more than 80%
                nGs = Gc # growth and consumption rates become that of coop
                # nCs = Cc
            else:
                nGs = Gs
                # nCs = Cs
            r0 = genoResource(g0, g2, Gc, Cc, nGs, Cs, R_small)
            r2 = R_small - r0

            # reproduce genotypes
            g0 = int(genoNumber(g0, r0, Cc, K)) # g0
            g2 = int(genoNumber(g2, r2, Cs, K)) # g2
            # r0 = float(g0)/(g0+g2)
            # r2 = float(g2)/(g0+g2)
            # totalI = g0 + g2
            # pScale = totalI/4.0
            # g0 = g0/pScale
            # g2 = g2/pScale
        # g0 = g0t - g0
        # g2 = g2t - g2
        #4# Migrant pool formation (dispersal): Return the progeny of each group to the migrant pool.
        genotypes[0] = genotypes[0] + g0
        genotypes[2] = genotypes[2] + g2
        # print "genotypes " + str(genotypes)

    # reproduce large groups
    for j in range(largeGroups.shape[0]):
        g1 = 0 # number of genotype 0
        g3 = 0 # number of genotype 2
        for i in range(40):
            if largeGroups[j][i] == 1:
                g1 += 1
            else:
                g3 += 1
        totalI = g1 + g3
        for x in range(t):
            # compute resource consumed
            # Rn = consume(totalI)

            print "selfish " + str(pool[3]/(pool[1]+pool[3]))
            if (pool[3]/(pool[1]+pool[3]) >= 0.5): # if cheaters more than 80%
                nGs = Gc # growth and consumption rates become that of coop
                # nCs = Cc
                print "Cc "
            else:
                nGs = Gs
                # nCs = Cs
                print "Cs "
            r1 = genoResource(g1, g3, Gc, Cc, nGs, Cs, R_large)
            r3 = R_large - r1
            # rScale = Rn/5.0
            # p1 = r1/rScale # ratio in portion
            # p3 = r3/rScale
            # scale = (g1+g3)/4.0 # scale down to total 4
            # g1 = g1/scale
            # g3 = g3/scale

            # reproduce genotypes
            g1 = int((genoNumber(g1, r1, Cc, K))) # g1
            g3 = int((genoNumber(g3, r3, Cs, K))) # g3
            # g1 = g1 * scale#scale back up
            # g3 = g3 * scale
            # r1 = float(g1)/(g1+g3)
            # r3 = float(g3)/(g1+g3)
            # totalI = g1 + g3
            # pscale = totalI/40.0
            # g1 = g1/pScale
            # g3 = g3/pScale
        # g1 = g1t - g1
        # g3 = g3t - g3
    #4# Migrant pool formation (dispersal)
        genotypes[1] = genotypes[1] + g1
        genotypes[3] = genotypes[3] + g3
        # print "genotypes " + str(genotypes)

    #5# Maintaining the global carrying capacity: Rescale the migrant pool back to size N, retaining the proportion of individuals with each genotype.
    groupSum = np.sum(genotypes)
    # print groupSum
    d = float(groupSum)/N
    d = math.ceil(d)
    d = int(d)
    genotypes = np.array(genotypes, dtype=float)
    if (d != 0):
        genotypes = genotypes/d
    # print "genotypes= " + str(genotypes)
    genoSum = np.sum(genotypes)
    genoPlot[0] = genotypes[0]/genoSum
    genoPlot[1] = genotypes[1]/genoSum
    genoPlot[2] = genotypes[2]/genoSum
    genoPlot[3] = genotypes[3]/genoSum
    genotypes = np.array(genotypes, dtype=int)
    # print "genotypes= " + str(genotypes)

    largeNum = genoPlot[1] + genoPlot[3]
    type0.append(largeNum)
    selfNum = genoPlot[2] + genoPlot[3]
    type1.append(selfNum)

    geno0.append(genoPlot[0])
    geno1.append(genoPlot[1])
    geno2.append(genoPlot[2])
    geno3.append(genoPlot[3])

#6# Iteration: Repeat from step 2 onwards for a number of generations, T

# plt.plot(type0, label='Large group size')
# plt.plot(type1, label='Selfish resource usage')
# plt.legend(loc='center right', shadow=True, fontsize='xx-large')
# plt.xlabel('Generation')
# plt.ylabel('Global frequency')
# plt.savefig('plot.png', bbox_inches='tight')

plt.plot(geno0, label='cooperative + small')
plt.plot(geno1, label='cooperative + large')
plt.plot(geno2, label='selfish + small')
plt.plot(geno3, label='selfish + large')
# plt.legend(loc='lower right', shadow=True, fontsize='large')
plt.xlabel('Generation')
plt.ylabel('Global genotype frequency')
plt.savefig('plot.png', bbox_inches='tight')
plt.show()
