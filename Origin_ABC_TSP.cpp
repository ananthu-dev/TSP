/* ABC algorithm coded using C programming language */

/* Artificial Bee Colony (ABC) is one of the most recently defined algorithms by Dervis Karaboga in 2005,
motivated by the intelligent behavior of honey bees. */




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<time.h>
#include<unistd.h>


/* Control Parameters of ABC algorithm*/
#define NP 20 /* The number of colony size (employed bees+onlooker bees)*/
#define FoodNumber NP/2 /*The number of food sources equals the half of the colony size*/
#define limit 100  /*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
#define maxCycle 2500 /*The number of cycles for foraging {a stopping criteria}*/

/* Problem specific variables*/
#define D 20 /*The number of parameters of the problem to be optimized*/
#define lb -100 /*lower bound of the parameters. */
#define ub 100 /*upper bound of the parameters. lb and ub can be defined as arrays for the problems of which parameters have different bounds*/


#define runtime 30  /*Algorithm can be run many times in order to see its robustness*/

int cities[D][D];
double Foods[FoodNumber][D]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
double Foods1[FoodNumber][D]; // For TSP
double f[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
double f1[FoodNumber];
double fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
double fitness1[FoodNumber];
double trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
double prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
double solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
double solution1 [D];
double ObjValSol; /*Objective function value of new solution*/
double FitnessSol; /*Fitness value of new solution*/
int neighbour, param2change; /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
double GlobalMin; /*Optimum solution obtained by ABC algorithm*/
double GlobalParams[D]; /*Parameters of the optimum solution*/
double GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/
double r; /*a random number in the range [0,1)*/

/*a function pointer returning double and taking a D-dimensional array as argument */
/*If your function takes additional arguments then change function pointer definition and lines calling "...=function(solution);" in the code*/
typedef double (*FunctionCallback)(double sol[D]);

/*benchmark functions */
double sphere(double sol[D]);
double Rosenbrock(double sol[D]);
double Griewank(double sol[D]);
double Rastrigin(double sol[D]);

/*Write your own objective function name instead of sphere*/
FunctionCallback function = &sphere;

/*Fitness function*/
double CalculateFitness(double fun)
 {
	 double result=0;
	 if(fun>=0)
	 {
		 result=1/(fun+1);
	 }
	 else
	 {
		 result=1+fabs(fun);
	 }
	 return result;
 }

/*The best food source is memorized*/
void MemorizeBestSource()
{
   int i,j;

	for(i=0;i<FoodNumber;i++)
	{
	if (fitness1[i]<GlobalMin)
		{
        GlobalMin=fitness1[i];
        for(j=0;j<D;j++)
           GlobalParams[j]=Foods1[i][j];
        }
	}
}
void SPV1(double sol[D])
{
    int i,c=0,j,minindex=0;
    double min=0;
    for(j=0;j<D;j++)
    {
    minindex=0;
    min=sol[0];
    for(i=0;i<D;i++)
    {
        if(sol[i]<min)
        {
            min=sol[i];
            minindex=i;
        }
    }
    solution1[minindex]=c;
    sol[minindex]=10000000;
    c++;
    }
}

void SPV(int index)
{
    int i,c=0,j,minindex=0;
    double min=0;
    for(j=0;j<D;j++)
    {
    minindex=0;
    min=Foods[index][0];
    for(i=0;i<D;i++)
    {
        if(Foods[index][i]<min)
        {
            min=Foods[index][i];
            minindex=i;
        }
    }
    Foods1[index][minindex]=c;
    Foods[index][minindex]=10000000;
    c++;
    }
}

int CalculateFitness1(double sol[D])
{
    int i,sum=0;
    for(i=0;i<D-1;i++)
    {
        sum+=cities[(int)sol[i]][(int)sol[i+1]];
    }
    sum+=cities[(int)sol[D-1]][(int)sol[0]];
    return sum;
}
/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */
/* Counters of food sources are also initialized in this function*/
void init(int index)
{
   int j;
   for (j=0;j<D;j++)
		{
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        Foods[index][j]=r*(ub-lb)+lb;
		solution[j]=Foods[index][j];
		}
		SPV(index);
		for(j=0;j<D;j++)
		{
		    solution1[j]=Foods1[index][j];
		}
	f[index]=function(solution);
	fitness1[index]=CalculateFitness1(solution1);
	fitness[index]=CalculateFitness(f[index]);
	trial[index]=0;
}

/*All food sources are initialized */
void initial()
{
	int i;
	for(i=0;i<FoodNumber;i++)
	{
	init(i);
	}
	//GlobalMin=f[0];
	GlobalMin=fitness1[0];
    for(i=0;i<D;i++)
    {
    //GlobalParams[i]=Foods[0][i];
    GlobalParams[i]=Foods1[0][i];
    }

}

void SendEmployedBees()
{
  int i,j;
  /*Employed Bee Phase*/
   for (i=0;i<FoodNumber;i++)
        {
        /*The parameter to be changed is determined randomly*/
        r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        param2change=(int)(r*D);

        /*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);

        /*Randomly selected solution must be different from the solution i*/
        while(neighbour==i)
        {
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);
        }
        for(j=0;j<D;j++)
        solution[j]=Foods[i][j];

        /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        solution[param2change]=Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

        /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        if (solution[param2change]<lb)
           solution[param2change]=lb;
        if (solution[param2change]>ub)
           solution[param2change]=ub;
        ObjValSol=function(solution);
        SPV1(solution);
        FitnessSol=CalculateFitness1(solution1);

        /*a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>fitness1[i])
        {
        /*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
        trial[i]=0;
        for(j=0;j<D;j++)
        Foods1[i][j]=solution1[j];
        f[i]=ObjValSol;
        fitness1[i]=FitnessSol;
        }
        else
        {   /*if the solution i can not be improved, increase its trial counter*/
            trial[i]=trial[i]+1;
        }


        }

        /*end of employed bee phase*/

}

/* A food source is chosen with the probability which is proportioal to its quality*/
/*Different schemes can be used to calculate the probability values*/
/*For example prob(i)=fitness(i)/sum(fitness)*/
/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
void CalculateProbabilities()
{
     int i;
     double maxfit;
     maxfit=fitness1[0];
  for (i=1;i<FoodNumber;i++)
        {
           if (fitness1[i]<maxfit)
           maxfit=fitness1[i];
        }

 for (i=0;i<FoodNumber;i++)
        {
         prob[i]=(0.9*(fitness1[i]/maxfit))+0.1;
        }

}

void SendOnlookerBees()
{

  int i,j,t;
  i=0;
  t=0;
  /*onlooker Bee Phase*/
  while(t<FoodNumber)
        {

        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        if(r<prob[i]) /*choose a food source depending on its probability to be chosen*/
        {
        t++;

        /*The parameter to be changed is determined randomly*/
        r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        param2change=(int)(r*D);

        /*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);

        /*Randomly selected solution must be different from the solution i*/
        while(neighbour==i)
        {
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        neighbour=(int)(r*FoodNumber);
        }
        for(j=0;j<D;j++)
        solution[j]=Foods[i][j];

        /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
        r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
        solution[param2change]=Foods[i][param2change]+(Foods[i][param2change]-Foods[neighbour][param2change])*(r-0.5)*2;

        /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
        if (solution[param2change]<lb)
           solution[param2change]=lb;
        if (solution[param2change]>ub)
           solution[param2change]=ub;
        ObjValSol=function(solution);
        SPV1(solution);
        FitnessSol=CalculateFitness1(solution1);

        /*a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>fitness1[i])
        {
        /*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
        trial[i]=0;
        for(j=0;j<D;j++)
        Foods1[i][j]=solution[j];
        f[i]=ObjValSol;
        fitness1[i]=FitnessSol;
        }
        else
        {   /*if the solution i can not be improved, increase its trial counter*/
            trial[i]=trial[i]+1;
        }
        } /*if */
        i++;
        if (i==FoodNumber-1)
        i=0;
        }/*while*/

        /*end of onlooker bee phase     */
}

/*determine the food sources whose trial counter exceeds the "limit" value. In Basic ABC, only one scout is allowed to occur in each cycle*/
void SendScoutBees()
{
int maxtrialindex,i;
maxtrialindex=0;
for (i=1;i<FoodNumber;i++)
        {
         if (trial[i]>trial[maxtrialindex])
         maxtrialindex=i;
        }
if(trial[maxtrialindex]>=limit)
{
	init(maxtrialindex);
}
}


/*Main program of the ABC algorithm*/
int main()
{
int iter,run,i,j;
double mean;
mean=0;
srand(time(NULL));

//for(run=0;run<runtime;run++)
//{

//for(i=0;i<D;i++)
//{
//    for(j=0;j<D;j++)
//    {
//        scanf("%d",&cities[i][j]);
//    }
//}

for(i=0; i<D; ++i)
      for(j=0; j<D; ++j)
         cities[i][j] = 1; 

for(i=0; i<D; ++i){
    for(j=0; j<D; ++j) 
        if(cities[i][j] == 1)
        {
            cities[i][j] = 2;//1+rand()%100;
            cities[j][i] = 2;//cities[i][j];
        }
}
for(i=0;i<D;i++)
    cities[i][i]=0;

initial();
MemorizeBestSource();
clock_t t=clock();
for (iter=0;iter<maxCycle;iter++)
    {
    SendEmployedBees();
    CalculateProbabilities();
    SendOnlookerBees();
    MemorizeBestSource();
    SendScoutBees();
    }

    t=clock()-t;
//for(j=0;j<D;j++)
//		{
//			printf("GlobalParam[%d]: %f\n",j+1,GlobalParams[j]);
//		}

printf("Path is : \n\n");
for(j=0;j<D;j++)
		{
			printf("%d ",(int)GlobalParams[j]);
		}
		printf("%d",(int)GlobalParams[0]);
        printf("\n\n");
		printf("No. of Clicks is %d and time in sec is %f",t,(float)t/CLOCKS_PER_SEC);
        printf("\n");
//printf("%d. run: %e \n",run+1,GlobalMin);
//GlobalMins[run]=GlobalMin;
//mean=mean+GlobalMin;
//}
//mean=mean/runtime;
//printf("Means of %d runs: %e\n",runtime,mean);
//getchar();
//usleep(1000);
}


double sphere(double sol[D])
{
int j;
double top=0;
for(j=0;j<D;j++)
{
top=top+sol[j]*sol[j];
}
return top;
}

double Rosenbrock(double sol[D])
{
int j;
double top=0;
for(j=0;j<D-1;j++)
{
top=top+100*pow((sol[j+1]-pow((sol[j]),(double)2)),(double)2)+pow((sol[j]-1),(double)2);
}
return top;
}

 double Griewank(double sol[D])
 {
	 int j;
	 double top1,top2,top;
	 top=0;
	 top1=0;
	 top2=1;
	 for(j=0;j<D;j++)
	 {
		 top1=top1+pow((sol[j]),(double)2);
		 top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);

	 }
	 top=(1/(double)4000)*top1-top2+1;
	 return top;
 }

 double Rastrigin(double sol[D])
 {
	 int j;
	 double top=0;

	 for(j=0;j<D;j++)
	 {
		 top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
	 }
	 return top;
 }

