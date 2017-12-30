import numpy as np
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

def genoResource(): # group resource the genotype received

    return r

def genoNumber(): # number of individuals in a group with the genotype

    return n

## Initialise the migrant pool with N individuals.
# 4 possible genotypes
genoSize = N/4
genotypes = [genoSize, genoSize, genoSize, genoSize] # total = N
# 1 cooperative + small
# 2 cooperative + large
# 3 selfish + small
# 4 selfish + large

## Group formation (aggregation): Assign individuals in the migrant pool to groups
# small size groups, size 4
smallGroups = [] # previous all small groups
smallGroup = []
for _ in range(4):
    smallGroup.append(0)
print smallGroup
cont = 1
i = 0
while (cont == 1) and (genotypes[0] > 0) and (genotypes[2] > 0):
    print "i= " + str(i)
    for j in range(4):
        flip = random()
        print "random= " + str(flip)
        if (flip < 0.5):
            if (genotypes[0] > 0):
                print "coop"
                smallGroup[j] = 0 # put coop into group
                genotypes[0] = genotypes[0] - 1
                print "genotypes[0]= " + str(genotypes[0])
            else:
                cont = 0
                break
        else:
            if (genotypes[2] > 0):
                print "selfish"
                smallGroup[j] = 2 # put selfish into group
                genotypes[2] = genotypes[2] - 1
                print "genotypes[2]= " + str(genotypes[2])
            else:
                cont = 0
                break
    print smallGroup
    smallGroups = np.append(smallGroups, smallGroup) # append into all small groups
    # if (i == 10):
    #     break
    # print "smallGroups" + str(smallGroups)
    i += 1

smallGroups = np.reshape(smallGroups, (-1, 4))
print smallGroups
print genotypes

# large size groups, size 40
largeGroups = [] # all large grops
largeGroup = [None] * 40

# Reproduction: Perform reproduction within groups for t time-steps
# Migrant pool formation (dispersal): Return the progeny of each group to the migrant pool.
# Maintaining the global carrying capacity: Rescale the migrant pool back to size N, retaining the proportion of individuals with each genotype.
# Iteration: Repeat from step 2 onwards for a number of generations, T
