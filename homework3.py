import random
from math import sqrt, pow, floor
import numpy as np
from sys import maxsize
from time import time

def rouletteWheelRank(population, fitness):
    sumFitness = []
    start = 0
    for normal in fitness:
        start+=normal
        sumFitness.append(start)

    populationSize = len(population)
    parents = []
    for i in range(populationSize):
        randNum = random.uniform(0, 1)
        selected = 0
        for selectedStart in sumFitness:
            if(randNum<=selectedStart):
                parents.append(selected)
                break
            selected+=1
    selected = 0
    for p in parents:
        if fitness[p] > fitness[selected]:
            selected = p
    return population[selected]

def tournamentSelection(population, fitness, k=100):
    best = 0
    for i in range(1,k+1):
        ind = random.randint(0, len(population)-1)
        if fitness[ind] > fitness[best]:
            best = ind
    return population[best] 

def createInitialPopulation(size, sequence):
    initialPopulation = []
    for i in range(0,size):
        initialPopulation.append(random.sample(sequence, len(sequence)))
    return initialPopulation

def distanceFunction(n, sequence, locations):
    distance = 0.0
    for i in range(1,n):
        distance+= sqrt(pow(locations[sequence[i]][0]-locations[sequence[i-1]][0], 2)+ 
                        pow(locations[sequence[i]][1]-locations[sequence[i-1]][1], 2)+ 
                        pow(locations[sequence[i]][2]-locations[sequence[i-1]][2], 2)
        )
                            
    distance+=sqrt(pow(locations[sequence[0]][0]-locations[sequence[n-1]][0], 2)+
                   pow(locations[sequence[0]][1]-locations[sequence[n-1]][1], 2)+
                   pow(locations[sequence[0]][2]-locations[sequence[n-1]][2], 2)
    )
    return distance

def findPopulationFitness(n, population, locations):
    fitness =[]
    minDist = maxsize
    minSeq = population[0]
    for individual in population:
        individualDist = round(distanceFunction(n, individual, locations), 3)
        if(minDist>individualDist):
            minDist = individualDist
            minSeq = individual
        fitness.append(1/(individualDist+1))
    #normalize the fitness values
    sum =0.0
    for i in range(0,len(fitness)):
        sum +=fitness[i]
    fitness = np.array(fitness)/sum #fitness will be a probability between 0 and 1
    return fitness, minSeq[:], minDist
    
def oxCrossover(parent1, parent2, startIndex, endIndex):
    subsequence = parent1[startIndex: endIndex]
    subsequence2 = parent2[startIndex: endIndex]
    child1 = [0]*len(parent1)
    child2 = [0]*len(parent1)
    child1[startIndex: endIndex] = subsequence[:]
    child2[startIndex: endIndex] = subsequence2[:]
    i=0
    k=0
    for j in range(0, len(parent2)):
        while i>=startIndex and i<endIndex:
            i=i+1
        if(i>=len(parent2)):
            break
        if parent2[j] not in subsequence:
            child1[i] = parent2[j]
            i=i+1
    for j in range(0, len(parent2)):
        while k>=startIndex and k<endIndex:
            k=k+1
        if(k>=len(parent2)):
            break
        if parent1[j] not in subsequence2:
            child2[k] = parent1[j]
            k=k+1
    return child1[:], child2[:]

def mutateIndividual(individual, mutationRate=0.5):
    n = len(individual)
    for i in range(0, n):
        if(random.random() < mutationRate):
            firstIndex = random.randint(0, n-1)
            secondIndex = random.randint(0, n-1)

            temp = individual[firstIndex]
            individual[firstIndex] = individual[secondIndex]
            individual[secondIndex] = temp

    if(random.random() < mutationRate):
        firstIndex = random.randint(0,floor(n/2)-1)
        secondIndex = random.randint(firstIndex+1, floor(n/2))
        thirdIndex = random.randint(secondIndex+1, n-2)
        fourthIndex = random.randint(thirdIndex+1, n-1)

        individual = individual[0:firstIndex] + individual[thirdIndex:fourthIndex] + individual[secondIndex: thirdIndex] + individual[firstIndex: secondIndex] + individual[fourthIndex: n]

    else:
        firstIndex = random.randint(0, n-1)
        secondIndex = random.randint(0, n-1)

        individual[firstIndex:secondIndex] = random.sample(individual[firstIndex: secondIndex], len(individual[firstIndex: secondIndex]))

    return individual[:]

def checkSanity(individual):
    checkSet = set(individual)
    if(len(individual)!=len(checkSet)):
        raise Exception("Error: Duplicates")
    
def nextGenerationOfPopulation(iter, noOfCities, population, fitness):
    newPopulation = []
    for i in range(0, floor(len(population)/2)):
        # selectedParent1 = rouletteWheel(population, fitness)
        # selectedParent2 = rouletteWheel(population, fitness)
        selectedParent1 = tournamentSelection(population, fitness)
        selectedParent2 = tournamentSelection(population, fitness)
        # selectedParent1 = rouletteWheelRank(population, fitness)
        # selectedParent2 = rouletteWheelRank(population, fitness)
        startIndex = random.randint(0, noOfCities-1)
        endIndex = random.randint(startIndex+1, noOfCities)
        selectedParent1, selectedParent2 = oxCrossover(selectedParent1, selectedParent2, startIndex, endIndex)
        selectedParent1 = mutateIndividual(selectedParent1, 0.02)
        selectedParent2 = mutateIndividual(selectedParent2, 0.02)
        checkSanity(selectedParent1)
        checkSanity(selectedParent2)
        newPopulation.append(selectedParent1)
        newPopulation.append(selectedParent2)
    return newPopulation

def writeIntoFile(bestDist, sequence, locations):
    file = open('output.txt', 'w')
    if bestDist==-1:
        file.write(str(0))
        return
    file.write(str(bestDist)+'\n')
    for elem in sequence:
        file.write(str(locations[elem][0])+" "+str(locations[elem][1])+" "+str(locations[elem][2])+'\n')
    file.write(str(locations[sequence[0]][0])+" "+str(locations[sequence[0]][1])+" "+str(locations[sequence[0]][2])+'\n')

def readFrom(path):
    f = open(path, "r")
    f.readline()
    line = f.readline()
    locations = []
    while line:
        line = line[:-1]
        vector = line.split(" ")
        vector = [int(i) for i in vector]
        locations.append(vector)
        line = f.readline()
    f.close()
    
    return locations

def main():
    iterations = 20000
    locations = readFrom("input.txt")
    noOfCities = len(locations)
    populationSize = floor(noOfCities*1.3)
    if noOfCities==0:
        writeIntoFile(-1, [], locations)
        return
    if noOfCities==1:
        writeIntoFile(0, [0], locations)
        return
    timeLim=0
    if(noOfCities<=50):
        timeLim = 59.0
    elif(noOfCities<=100):
        timeLim = 74.0
    else:
        timeLim = 119.0

    sequence = [i for i in range(0,noOfCities)]
    population= createInitialPopulation(populationSize, sequence)
    fitness, bestEver, bestDist = findPopulationFitness(noOfCities, population, locations)

    numOfSameResult = 0
    sameResult = 0
    initialVal = 0
    startTime = time()
    for i in range(iterations):
        #print("Running iteration ", i)
        population = nextGenerationOfPopulation(i, noOfCities, population, fitness)
        fitness, currentBest, currentDist = findPopulationFitness(noOfCities, population, locations)
        if(currentDist == sameResult):
            numOfSameResult=numOfSameResult+1
        if(numOfSameResult==0):
            sameResult = currentDist
        if(numOfSameResult>5):
            break
        if(currentDist<bestDist):
            bestDist=currentDist
            bestEver=currentBest[:]
        if i==1:
            initialVal=bestDist
        
        checkDist  = round(distanceFunction(noOfCities, bestEver, locations),3)
        if(bestDist!=checkDist):
            raise Exception("Error: Wrong bestDist")
        
        #print("The length for the best Ever: ", bestEver, " is with the minimum distance of ", bestDist)
        timeElapsed = time()-startTime
        #print("Time taken: ", timeElapsed)
        if(timeElapsed > timeLim):
            break;
    #print("Initial -> Final: ", initialVal, " ", bestDist)
    writeIntoFile(bestDist, bestEver, locations)

if __name__ == "__main__":
    main()