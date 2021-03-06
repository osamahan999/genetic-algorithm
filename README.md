
graph of mutation chance vs average time unthreaded
https://docs.google.com/spreadsheets/d/18JoG38uJAAS2BCCCijN2On5PhXasAX27DEHTAxk79vY/edit#gid=0


This code has been refactored from my multithreaded genetic algorithm repo https://github.com/osamahan999/multithreaded-genetic-algorithm
I made this repo simply to see the different in speed between the multithreaded and a single threaded genetic algorithm, and unsurprisingly, 
the single threaded repo was way faster. I was making way too many threads, which I knew was not great, but I still wanted to do just to learn
how to use threading! 


Compile using gcc -o geneticAlgorithm geneticAlgorithm.c  

Run using ./geneticAlgorithm

What this program does is the following:

1. It initializes a best fit array which is our target. This is done with random numbers that are seeded using the current time, multiplied by some constants
which produce the max random value of 100000.

I chose this large number because I wanted a lot of variation in the numbers, and I wanted many chances of error so that I can develop the algorithm to be more efficient

2. It then initializes the population with random numbers, same way we initialized the best fit array
This is done without multithreading to compare to the multithreaded version

3. Then, utilizing a weighing function, I check to see which populations have the best value based on their errors. The weighing function then tells us a % that we should give to each population for a weighted random parent-choosing process. 

I also implemented a TOP_X parent choosing process. This only weighs the top X amout of individuals, and picks parents from those. 
I implemented this by using seleciton sort because i was lazy and didnt want to implement qsort, and found that its much slower. Much much
slower, such that even if I were using qsort which is a factor faster, would not make a difference. 

4. Utilize the weights to pseudo randomly decide the parents, and then have them have children using a crossover function
What this crossover function does is it chooses a gene from each parents and flips a coin for which one goes into the child.

5. We then mutated the child, with each index having some MUTATION_CHANCE percent of mutating by either 1.1 or .9 which are chosen with another coin flip

6. We do this for each generation

A quick notesheet runtime analysis produces
o(r) + o(nr^2) + 2 * o(n) + o(t * (6n + 2n^2 + 3nr 2nr^3))

where r is the amount of indices in our fitness array, n is the amount of indiviudals in the population, and t is the amount of generations we run.

