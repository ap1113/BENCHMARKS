#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include "carbon_user.h"
#include <time.h>
#include <sys/timeb.h>

#define BILLION 1E9

typedef struct		//define a structure (class) with attributes
{
  int*      Q;
  //int**     W;
  // int     W_index[16384][32];
  int       tid;
  int       P;
  long       N;
  int       DEG;
  pthread_barrier_t* barrier;
} thread_arg_t;

// HERERERERERE //////
int start = 64;
int next_source = -1;
thread_arg_t thread_arg[1024];
pthread_t thread_handle[1024];
 // int* avg;
 // double* delta;

pthread_barrier_t  barrier;       //barrier to synchronize threads
pthread_mutex_t id;               //single lock for 1 variable
pthread_mutex_t* locks;   			//locks for each vertex in the graph
int tid_g = -1;                     //to increment global thread IDs later

int* D;              
int* Q;
//SPECIFIC TO BFS
  // int change = 0;
  // int *test;
  // int *test1;
  int *temporary;
volatile int terminate = 0;
int** W;                          //graph weights (not used here, but may be random for other graphs)
//int**     sigma; 
// int Total_tid[1024] = {0};        //triangles for each thread


long N=262144;                          //global, otherwise softbound complains
int DEG=16;
int P=0;
// int Total = 0;   

long W_index[262144][16] = {0};   //graph (adjacency list)           


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
  // pthread_mutex_lock(&id);
	// tid_g++;
  // int tid = tid_g;
   int tid=0;
  tid = __sync_fetch_and_add(&tid_g,1);
  tid++;
  printf("\n TID: %d",tid);
  // pthread_mutex_unlock(&id);

	// pthread_barrier_wait(&barrier);

   //divide graph among threads
	long start = tid    *  (N) / (P);
	long stop = (tid+1) *  (N) / (P);
	// int node = 0;
	
    long iter = 0;
   
   // pthread_barrier_wait(&barrier);
   printf("\n %d b0",tid);
   barrier_wait();
	 printf("\n %d %ld %ld",tid,start,stop);

   //add edges to D array

	while(terminate==0)
      {   
        for(long uu=start; uu<stop; uu++)
        {
          // if(test1[uu]==0)
            // continue;
          //if(D[uu]==0 || D[uu]==2) //THIS
           // continue;				//THIS
           
          for(int i = 0; i < DEG; i++)
          {   
            long neighbor = W_index[uu][i];
            // pthread_mutex_lock(&locks[neighbor]);
            if(Q[neighbor]==1)
              Q[neighbor]=0;
            temporary[neighbor] = 1;
            // pthread_mutex_unlock(&locks[neighbor]);
          }
          
          // //termination condition
          if(Q[N-1]==0 || iter>=N) //largest vertex done, or all vertices done
            terminate=1;	
        }
        // // pthread_barrier_wait(&barrier);
        printf("\n %d T=%d b1",tid, terminate);
        // barrier_wait();
        for(long uu=start;uu<stop;uu++)
        {
          D[uu] = temporary[uu];
        }
        iter++;
        // // pthread_barrier_wait(&barrier);
        // // barrier_wait();
      }
	
	
    printf("\n %d bf",tid);
	// pthread_barrier_wait(&barrier);
    barrier_wait();
   return 0;
}

int main(int argc, char *argv[])
{

   long N1 = 8192; //atoi(argv[2]);  //1024;         //input arguments
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
   int Q1[N];
   if(!Q) {Q = Q1;}
   //delta = (double*) malloc((N)*sizeof(double));   
  
   // // test = (int*) malloc((N)*sizeof(int));
   // int test11[N];
   // if(!test) {test = test11;}

   // test1 = (int*) malloc((N)*sizeof(int));
   // int test21[N];
   // if(!test1) {test1 = test21;}
   

   // temporary = (int*) malloc((N)*sizeof(int));
   int temporary1[N];
   if(!temporary) {temporary = temporary1;}
   
   // //W = (int**)malloc((N)*sizeof(int*));
    // W_index = (int **)malloc(sizeof(int *) * N);
       

    // return 1;

   locks = new pthread_mutex_t[N];

   // for (int i = 0; i < N; i++)
   // {
	  // // //W[i] = (int*) malloc((DEG)*sizeof(int));
        // // W_index[i] = (int*) malloc((DEG)*sizeof(int));
        
   // }

pthread_mutex_init(&id, NULL);                         //barrier and lock inits
pthread_barrier_init(&barrier, NULL, P);
    for(long i=0; i<N; i++)
		pthread_mutex_init(&locks[i], NULL);

// //fill D and Q
for(long i = 0; i < N; i++)
{
    D[i] = 0;
    Q[i] = 1;
	// if(i < N-1){
		// // test[i]=DEG;
		// // test1[i]=1;
		// Total++;
	// }
	temporary[i]=0;	
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
          // for(int i = 0; i < N; i++)
          // {
                  // for(int j = 0; j < DEG; j++)
                  // {
                          // double v = drand48();
  
                          // if(W_index[i][j] == i)
                                  // W[i][j] = 0;
  
                          // else
                                  // W[i][j] = (int) (v*100) + 1;	
                  // }
          // }
		  
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

   // Disable performance and energy models
   
   
   	clock_gettime(CLOCK_REALTIME, &requestEnd);
	double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;
	printf( "\nTime:%lf\n", accum );   //prints time
	
	CarbonDisableModels();

   printf("BFS: Done.\n");
   
   // Free memory

   //CarbonStopSim();
   return 0;
}
