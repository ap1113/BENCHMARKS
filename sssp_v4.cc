#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include "carbon_user.h"
#include <time.h>
#include <sys/timeb.h>

#define BILLION 1E9
#define INT_MAX        100000000
//I quickly run out of stack space for any more than 2 threads

typedef struct		//define a structure (class) with attributes
{
  // int*      Q;
  // int**     W;
  // int**     W_index;
  int       tid;
  int       P;
  int       N;
  int       DEG;
  pthread_barrier_t* barrier;
} thread_arg_t;

// HERERERERERE //////
// int start = 64;
// int next_source = -1;
thread_arg_t thread_arg[1024];
pthread_t thread_handle[1024];
 // int* avg;
 // double* delta;

pthread_barrier_t  barrier;       //barrier to synchronize threads
pthread_mutex_t id;               //single lock for 1 variable
pthread_mutex_t* locks;   			//locks for each vertex in the graph
int tid_g = -1;                     //to increment global thread IDs later

int* D;              
// int* Q;
// int** W;                          //graph weights (not used here, but may be random for other graphs)
// int** W_index;  
long W[131072][16] = {0};   //graph (adjacency list)  
long W_index[131072][16] = {0};   //graph (adjacency list)   
//int Total_tid[1024] = {0};        //triangles for each thread
int P_max=2; //256;
int terminate = 0;
long range=1;


int N=131072;                          //global, otherwise softbound complains
int DEG=16;
int P=0;
//long long Total = 0;              //Total triangles (a large number so long long)


volatile int passed = 0; // Number of barriers, passed by all threads.
int bar=0;
int PMAX=0;
void* barrier_wait();		//Initialize wait function
void* threadMain(void* arg);                //function initialization

void* barrier_wait()
{
    int passed_old = passed; // Should be evaluated before incrementing *bar*!

    if(__sync_fetch_and_add(&bar,1) == (PMAX - 1)) // Not happening ...
    {   
        // The last thread, faced barrier.
        bar = 0;
        // *bar* should be reseted strictly before updating of barriers counter.
        __sync_synchronize(); 
        passed++; // Mark barrier as passed.
    }   
    else
    {   
        // Not the last thread. Wait others.
        while(passed == passed_old) {}; 
        // Need to synchronize cache with other threads, passed barrier.
        __sync_synchronize();
    } 
	return 0;
}

void* threadMain(void* arg)                 //main parallel function
{
  //each thread gets it IDs
  int tid=0;
  tid = __sync_fetch_and_add(&tid_g,1);
  tid++;
  // pthread_mutex_lock(&id);
	// tid_g++;
  // int tid = tid_g;
  // pthread_mutex_unlock(&id);
  
  
  int local_count          = N;
  int i, j, po;
  long uu = 0;

  int cntr = 0;
  // int i_start =  0;  //tid    * DEG / (arg->P);
  // int i_stop  = 0;   //(tid+1) * DEG / (arg->P);
  int start = 0;
  int stop = 1;
  int neighbor=0;
    printf("\n %d %d %d",tid,start,stop);
  // pthread_barrier_wait(&barrier);
  barrier_wait();

	while(terminate==0){
	  while(terminate==0)
	  {
		for(uu=start;uu<stop;uu++)
		{

		// if(test[uu]==0)
		  // continue;

		for(int i = 0; i < DEG; i++)
		{
		  if(uu<N)
			neighbor = W_index[uu][i];

		  if(neighbor>=N)
			break;

		  //pthread_mutex_lock(&locks[neighbor]);



				//relax
		  if((D[W_index[uu][i]] > (D[uu] + W[uu][i])))
			D[W_index[uu][i]] = D[uu] + W[uu][i];

		  //pthread_mutex_unlock(&locks[neighbor]);
		}
		}

		// pthread_barrier_wait(&barrier);
        // barrier_wait();
			
		if(tid==0)
		{
		  //range heuristic here
		  //old_range=range;
		  range = range*DEG; //change this for range heuristic e.g. range = range+DEG; 
		  if(range>=N)
			range=N;
		}

		// pthread_barrier_wait(&barrier);
        barrier_wait();


		start = tid * (range/P);
		stop = (tid+1) * (range/P);

		if(stop>range)
		  stop=range;	

		//{ pthread_mutex_lock(&lock);
			if(start==N || uu>N-1)
			  terminate=1;
		//} pthread_mutex_unlock(&lock);

		// pthread_barrier_wait(&barrier);
        // barrier_wait();

		printf("\n TID:%d   start:%d stop:%d terminate:%d",tid,start,stop,terminate);
	  }
	  // pthread_barrier_wait(&barrier);
      // barrier_wait();
	  if(tid==0)
	  {
		printf("\nAt 0");
		cntr++;
		if(cntr<P_max)
		{
		  terminate=0;
		  range=1;
		}
	  }
	  start=0;
	  stop=1;
	  // pthread_barrier_wait(&barrier);
      barrier_wait();
	}
  
  // barrier_wait();
   return 0;
}

int main(int argc, char *argv[])
{

   int N1 = 131072; //atoi(argv[2]);  //1024;         //input arguments
   int DEG1 = 16; //atoi(argv[3]); //16;
   int P1 = atoi(argv[1]); //8;

     // N=N1;
	 // DEG=DEG1;
	 P=P1;
     PMAX = P;
	 if (DEG > N)
	{
		fprintf(stderr, "DEG cannot be greater than N (0 P N DEG)\n");
		exit(EXIT_FAILURE);
	}
  
   // D = (int*) malloc((N)*sizeof(int));            //memory allocations
   int D1[N];
   if(!D) {D = D1;}
   // Q = (int*) malloc((N)*sizeof(int));
   // int Q1[N];
   // if(!Q) {Q = Q1;}
   // W = (int**)malloc((N)*sizeof(int*));
   // W_index = (int**)malloc((N)*sizeof(int*));
   locks = new pthread_mutex_t[N];

   // for (int i = 0; i < N; i++)
   // {
	  // W[i] = (int*) malloc((DEG)*sizeof(int));
	  // W_index[i] = (int*) malloc((DEG)*sizeof(int));
   // }

pthread_mutex_init(&id, NULL);                         //barrier and lock inits
pthread_barrier_init(&barrier, NULL, P);
    for(long i=0; i<N; i++)
		pthread_mutex_init(&locks[i], NULL);

//fill D and Q
for(long i = 0; i < N; i++)
{
    D[i] = 0;
    // Q[i] = 1;
	for(int j=0;j<DEG;j++)
    {
      // W[i][j] = 100000;
      // W_index[i][j] = 100000;
    }
}
	
//Now fill the graph with synthetic stuff (same as init_weights in orig)

          int range = DEG + (N - DEG)/16;
  
          // Initialize to -1
          for(long i = 0; i < N; i++)
                  for(int j = 0; j < DEG; j++)
                          W_index[i][j]= -1;
  
          // Populate Index Array
          for(long i = 0; i < N; i++)
          {
                  long last = 0;
                  int min = 0;
                  int max = DEG;
                  for(int j = 0; j < DEG; j++)
                  {
                          if(W_index[i][j] == -1)
                          {        
								long neighbor = i+j;																	
                                  if(neighbor > last)
                                  {
                                          W_index[i][j] = neighbor;
                                          last = W_index[i][j];
                                  }
                                  else
                                  {
                                          if(last < (N-1))
                                          {
                                                  W_index[i][j] = (last + 1);
                                                  last = W_index[i][j];
                                          }
                                  }
                          }
                          else
                          {
                                  last = W_index[i][j];
                          }
                          if(W_index[i][j]>=N)
                          {
                                  W_index[i][j] = N-1;
                          }
                  }
          }
  
          // Populate Cost Array
          for(long i = 0; i < N; i++)
          {
                  for(int j = 0; j < DEG; j++)
                  {
                          double v = drand48();
  
                          if(W_index[i][j] == i)
                                  W[i][j] = 0;
  
                          else
                                  W[i][j] = (int) (v*100) + 1;	
                  }
          }
		  
//Graph filling finished


   // Enable performance and energy models
   CarbonEnableModels();
   
   struct timespec requestStart, requestEnd;
   clock_gettime(CLOCK_REALTIME, &requestStart);  //to measure time

   pthread_t thread_handles[P];              //spawning threads
   for (int i = 1; i < P; i++)
   {
     pthread_create(thread_handle+i, NULL, threadMain, (void*)&thread_arg[i]); 
   }
	 threadMain((void*) &thread_arg[0]);

   /*for (int i = 1; i < P; i++)
   {
      pthread_join(thread_handles[i], NULL);
   }*/

   clock_gettime(CLOCK_REALTIME, &requestEnd);
	  double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;
	  printf( "\nTime:%lf\n", accum );   //prints time
   
   
   // Disable performance and energy models
   CarbonDisableModels();
   
   	


   printf("Single Source: Done.\n");
   
   // Free memory

   //CarbonStopSim();
   return 0;
}
