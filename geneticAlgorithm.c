#define _GNU_SOURCE
#define NUM_OF_POPULATION 500
#define NUM_OF_FITNESS_INDICES 1
#define NUM_OF_GENS 2000
#define NUM_OF_PARENTS (2 * NUM_OF_POPULATION) //amt of parents each child has
#define AMT_OF_ERROR_PER_INDICE 100            //with our # range going from 0 to 100,000, a 100 error per indice means 1/1000 error
#define TOP_X 10

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

typedef struct
{
    double fitness[NUM_OF_FITNESS_INDICES];
    double fitnessIndex;

    double weight;

} individual;

void randomNumArr(double *arr);
void initializePopulation();
void printArr(double arr[]);
double fitnessComparison(double individual[], double goal[]);
void parentFunc();
void haveChildren();
double getError(int popNum);
void calculatePopulationWeights();
double algorithmInitialization();
void bestMutationChance(int runsToGetAvg, char filename[]);
void populationWeightTopIndividuals();
void pushArrayDownOneIndex(individual *topIndividual, int pushDownAmt);

int MUTATION_CHANCE = 0; // mutation_chance% of mutation, if 20, then 80

double bestFit[NUM_OF_FITNESS_INDICES];
individual *population[NUM_OF_POPULATION];
individual *initialPopulation[NUM_OF_POPULATION]; //duplicates population array to free the pointers

individual *parents[NUM_OF_POPULATION * 2];

int main()
{

    //dealloc
    for (int i = 0; i < NUM_OF_POPULATION; i++)
        free(initialPopulation[i]);
}

/**
 * Runs through a bunch of mutation chances to find best mutation chance with the current crossover function 
 * Takes in the amt of runs to get the average, and a string filename like text.txt to write to
*/
void bestMutationChance(int runsToGetAvg, char filename[])
{
    int bestMutationChance, initialMutationChance = MUTATION_CHANCE;

    FILE *fp;
    fp = fopen(filename, "w");

    if (fp == NULL)
        exit(1);
    fputs("# of indices,\t # of gens,\t # of pop,\t sample size,\t mutation chance,\t avgTime,\t netTime\n", fp);
    double bestRuntime = INFINITY;

    for (int i = 0; i < 100 - initialMutationChance; i++)
    {
        printf("run #%d with mutation chance %d\n", i, MUTATION_CHANCE);
        double totalTime = 0;

        for (int j = 0; j < runsToGetAvg; j++)
        {
            totalTime += algorithmInitialization();
        }

        double avgTime = totalTime / runsToGetAvg;
        printf("run #%d with avg time %f with a net time of %f\n", i, avgTime, totalTime);

        fprintf(fp, "%d,\t %d,\t %d,\t %d,\t %d, %f,\t %f\n", NUM_OF_FITNESS_INDICES, NUM_OF_GENS, NUM_OF_POPULATION, runsToGetAvg, MUTATION_CHANCE, avgTime, totalTime);

        if (avgTime < bestRuntime)
        {
            bestRuntime = avgTime;
            bestMutationChance = MUTATION_CHANCE;
        }

        MUTATION_CHANCE++;
    }

    printf("for #%d indices, best mutation chance is %d with avg runtime of %f\n", NUM_OF_FITNESS_INDICES, bestMutationChance, bestRuntime);
    fclose(fp);
}

//initializes the algorithm and calls it for each generation
double algorithmInitialization()
{
    struct timeval start, end;
    gettimeofday(&start, NULL); //start timer

    randomNumArr(bestFit); //generates what we will be regarding the best fit array

    initializePopulation();
    printf("pregen individual 1's err= %f\n", getError(0));

    // printf("gen #%d %f\n", 0, getError());

    for (int i = 0; i < NUM_OF_GENS; i++)
    {

        parentFunc();

        haveChildren();

        //find lowest error in population and break
        double lowestErr = INFINITY;
        int popWithLowestErr;
        int genCounter = i;

        for (int j = 0; j < NUM_OF_POPULATION; j++)
        {

            double err = getError(j);
            if (err < lowestErr)
            {
                lowestErr = err;
                popWithLowestErr = j;
            }

            if (err < (AMT_OF_ERROR_PER_INDICE * NUM_OF_FITNESS_INDICES))
            {
                printf("gen #%d had err %f and is converged!\n", i, err);
                i = NUM_OF_GENS;
                break;
            }
        }

        printf("gen #%d pop %d lowest err: %f\n", genCounter, popWithLowestErr, lowestErr);
    }

    gettimeofday(&end, NULL); //end timer

    // printf("Time passed %f seconds \n", (end.tv_sec - start.tv_sec) + ((end.tv_usec - start.tv_usec) * 1.0 / 1000000));

    return (end.tv_sec - start.tv_sec) + ((end.tv_usec - start.tv_usec) * 1.0 / 1000000);
}

//get total absolute val error on avg just for debugging
double getError(int popNum)
{
    double error = 0;

    for (int j = 0; j < NUM_OF_FITNESS_INDICES; j++)
    {
        error += fabs(population[popNum]->fitness[j] - bestFit[j]);
    }

    return error;
}

/**
 * crosses the arrays of each 2 parents into a new child, 
 * with each index having a 50% chance of being from parent A or from parent B
*/
void haveChildren()
{

    int popIndex = 0; //used to update population
    struct timeval time;
    int t;
    int crossOverPoint;

    double tempA[NUM_OF_FITNESS_INDICES];
    double tempB[NUM_OF_FITNESS_INDICES];

    for (int i = 0; i < NUM_OF_PARENTS; i += 2)
    {

        //copy data into temp var
        for (int copyIndex = 0; copyIndex < NUM_OF_FITNESS_INDICES; copyIndex++)
        {
            tempA[copyIndex] = parents[i]->fitness[copyIndex];
            tempB[copyIndex] = parents[i + 1]->fitness[copyIndex];
        }

        //for each indice, gets the average of the two parent's index and sets it to the index in population
        for (int j = 0; j < NUM_OF_FITNESS_INDICES; j++)
        {
            gettimeofday(&time, NULL);
            t = time.tv_usec; //random math to maybe make it more chaotic?
            crossOverPoint = (rand_r(&t) % 100);

            if (crossOverPoint < 50)
                population[popIndex]->fitness[j] = tempA[j];
            else
                population[popIndex]->fitness[j] = tempB[j];
        }

        population[popIndex]->fitnessIndex = fitnessComparison(population[popIndex]->fitness, bestFit); //update fitness

        popIndex++;
    }
}

void initializePopulation()
{

    for (int i = 0; i < NUM_OF_POPULATION; i++)
    {
        //allocate mem for each population, store original points in initialPopulation array to free them at end of run
        initialPopulation[i] = (individual *)malloc(sizeof(individual));
        population[i] = initialPopulation[i];

        // //initialize the indiivdual
        randomNumArr(population[i]->fitness);
        population[i]->fitnessIndex = fitnessComparison(population[i]->fitness, bestFit);
    }
}

//calculates a weight for each population, which is used to do weighted random reproduction
//for equal distribution among all individuals based on fitness
void calculatePopulationWeights()
{
    double totalSum = 0;

    for (int i = 0; i < NUM_OF_POPULATION; i++)
        totalSum += population[i]->fitnessIndex;

    for (int i = 0; i < NUM_OF_POPULATION; i++)
        population[i]->weight = population[i]->fitnessIndex / totalSum;
}

void parentFunc()
{
    calculatePopulationWeights();

    double num;
    struct timeval time;
    int t;

    for (int parent = 0; parent < NUM_OF_PARENTS; parent++)
    {

        num = 0; //sum of weights to find which parent to use

        gettimeofday(&time, NULL);
        t = time.tv_usec;

        double newParent = ((double)(rand_r(&t) % 100000)) / 100000; //0 to 999999 because our random doubles go to 0 to 100,000

        for (int j = 0; j < NUM_OF_POPULATION; j++)
        {
            num += population[j]->weight;

            if (newParent < num)
            {
                parents[parent] = population[j];

                break;
            }
        }
    }

    for (int parent = 0; parent < NUM_OF_PARENTS; parent++)
    {

        int mutationChance;

        for (int i = 0; i < NUM_OF_FITNESS_INDICES; i++)
        {
            gettimeofday(&time, NULL);
            t = time.tv_usec;

            mutationChance = rand_r(&t) % 100;
            if (mutationChance >= (100 - MUTATION_CHANCE)) //20% chance = 100 - 20
            {

                //50% chance to double, or 50% to half
                gettimeofday(&time, NULL);
                t = time.tv_usec;
                int binaryMutator = (rand_r(&t) % 2);

                if (binaryMutator == 0)
                    parents[parent]->fitness[i] *= 1.1; //increase it
                else
                    parents[parent]->fitness[i] *= .9; //decrease it
            }
        }
    }
}

// a return value of 0 means you hit the goal! higher than 0 means farther
double fitnessComparison(double individual[], double goal[])
{
    double fitnessFactor = 0;

    for (int i = 0; i < NUM_OF_FITNESS_INDICES; i++)
    {
        fitnessFactor += (fabs(individual[i] - goal[i]));
    }

    return (1 / fitnessFactor);
}

//initializes the array with pseudo random doubles
void randomNumArr(double *arr)
{

    //psuedo random seed with time
    struct timeval time;
    gettimeofday(&time, NULL);

    int t = time.tv_usec; //random math to maybe make it more chaotic?

    for (int i = 0; i < NUM_OF_FITNESS_INDICES; i++)
    {
        arr[i] = (rand_r(&t) % 100000) * 0.339 * 2.94988200619; //random double, more like random ints mapped to doubles. Max num is 100,000 exactly
    }
}

void printArr(double arr[])
{

    for (int i = 0; i < NUM_OF_FITNESS_INDICES; i++)
        printf("%f\n", arr[i]);
}

/**ALL CODE BENEATH HERE DOES NOT WORK AND IS COMPLETE GARBAGE
 * 
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * //////////////////////////////////////////////////////////
 * 
*/
void findParentsTopIndividuals(individual *topIndividuals[])
{
    double totalSum = 0;

    for (int i = 0; i < TOP_X; i++)
        totalSum += topIndividuals[i]->fitnessIndex;

    for (int i = 0; i < TOP_X; i++)
        topIndividuals[i]->weight = topIndividuals[i]->fitnessIndex / totalSum;
}

void populationWeightTopIndividuals()
{
    //top x individuals
    individual *topIndividuals[TOP_X]; //temporary array to hold all the top 10 individuals

    //initialize arr to NULL values
    for (int i = 0; i < TOP_X; i++)
    {
        topIndividuals[i] = NULL;
    }

    //for each individual, we get the fitnessIndex. We then go through each individual in our topIndividuals array
    //and check to see if our current individual has a larger fitnessIndex than any one in the array. If so,
    //we push down the array till that index, and put the new person in there.
    for (int i = 0; i < NUM_OF_POPULATION; i++)
    {
        double fitnessIndex = population[i]->fitnessIndex;
        for (int j = 0; j < TOP_X; j++)
        {
            if (topIndividuals[j] == NULL)
            {
                topIndividuals[j] = population[i];
            }
            else if (topIndividuals[j]->fitnessIndex > fitnessIndex)
            {
                for (int i = TOP_X - 1; i >= (TOP_X - j); j--)
                {
                    topIndividuals[i] = topIndividuals[i - 1];
                }
                topIndividuals[j] = population[i];
                break;
            }
        }
    }

    double totalSum = 0;

    for (int i = 0; i < TOP_X; i++)
        totalSum += topIndividuals[i]->fitnessIndex;

    for (int i = 0; i < TOP_X; i++)
        topIndividuals[i]->weight = topIndividuals[i]->fitnessIndex / totalSum;

    struct timeval time;
    double num;
    int t;

    for (int i = 0; i < NUM_OF_PARENTS; i++)
    {

        num = 0; //sum of weights to find which parent to use

        gettimeofday(&time, NULL);
        t = time.tv_usec;

        double newParent = ((double)(rand_r(&t) % 100000)) / 100000; //0 to 999999 because our random doubles go to 0 to 100,000

        for (int j = 0; j < TOP_X; j++)
        {
            num += topIndividuals[j]->weight;

            if (newParent < num)
            {
                parents[i] = topIndividuals[j];

                break;
            }
        }
    }

    for (int parent = 0; parent < NUM_OF_POPULATION; parent++)
    {
        int mutationChance;

        for (int i = 0; i < NUM_OF_FITNESS_INDICES; i++)
        {
            gettimeofday(&time, NULL);
            t = time.tv_usec;

            mutationChance = rand_r(&t) % 100;
            if (mutationChance >= (100 - MUTATION_CHANCE)) //20% chance = 100 - 20
            {

                //50% chance to double, or 50% to half
                gettimeofday(&time, NULL);
                t = time.tv_usec;
                int binaryMutator = (rand_r(&t) % 2);

                if (binaryMutator == 0)
                    parents[parent]->fitness[i] *= 1.1; //increase it
                else
                    parents[parent]->fitness[i] *= .9; //decrease it
            }
        }
    }

    int popIndex = 0; //used to update population

    int crossOverPoint;

    double tempA[NUM_OF_FITNESS_INDICES];
    double tempB[NUM_OF_FITNESS_INDICES];

    for (int i = 0; i < NUM_OF_PARENTS; i += 2)
    {

        //copy data into temp var
        for (int copyIndex = 0; copyIndex < NUM_OF_FITNESS_INDICES; copyIndex++)
        {
            tempA[copyIndex] = parents[i]->fitness[copyIndex];
            tempB[copyIndex] = parents[i + 1]->fitness[copyIndex];
        }

        //for each indice, gets the average of the two parent's index and sets it to the index in population
        for (int j = 0; j < NUM_OF_FITNESS_INDICES; j++)
        {
            gettimeofday(&time, NULL);
            t = time.tv_usec; //random math to maybe make it more chaotic?
            crossOverPoint = (rand_r(&t) % 100);

            if (crossOverPoint < 50)
                population[popIndex]->fitness[j] = tempA[j];
            else
                population[popIndex]->fitness[j] = tempB[j];
        }

        population[popIndex]->fitnessIndex = fitnessComparison(population[popIndex]->fitness, bestFit); //update fitness

        popIndex++;
    }
}

//pushes array down one index
void pushArrayDownOneIndex(individual *topIndividuals, int pushDownAmt)
{

    for (int i = TOP_X - 1; i >= (TOP_X - pushDownAmt); i--)
    {
        topIndividuals[i] = topIndividuals[i - 1];
    }
}
