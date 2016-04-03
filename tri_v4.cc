#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include "carbon_user.h"
#include <time.h>
#include <sys/timeb.h>

#define BILLION 1E9

typedef struct
{
	int*      Q;  
	// int**     W_index;
	int       tid;
	int       P;  
	long       N;  
    int       DEG;
	pthread_barrier_t* barrier;
} thread_arg_t;

thread_arg_t thread_arg[1024];
pthread_t thread_handle[1024];

pthread_barrier_t  barrier;       //barrier to synchronize threads
pthread_mutex_t id;               //single lock for 1 variable
pthread_mutex_t* locks;   			//locks for each vertex in the graph
int tid_g = -1;                     //to increment global thread IDs later

int* D;              
int* Q;
// int** W;                          //graph weights (not used here, but may be random for other graphs)
// int** W_index;                    //graph connections (adjacency list)
long W_index[262144][16] = {0}; 
// long W[262144][16] = {0}; 
int Total_tid[1024] = {0};        //triangles for each thread


long N=262144;                          //global, otherwise softbound complains
int DEG=16;
int P=0;
long long Total = 0;              //Total triangles (a large number so long long)


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
  printf("\n %d",tid);
  // pthread_mutex_unlock(&id);

	// pthread_barrier_wait(&barrier);
	barrier_wait();

   //divide graph among threads
	int start =  tid    *  (N) / (P);
	int stop =  (tid+1) *  (N) / (P);
   
   // pthread_barrier_wait(&barrier);
	barrier_wait();
	 printf("\n %d %d %d",tid,start,stop);

   //add edges to D array
   for(int uu=start;uu<stop;uu++)
   {
     for(int i = 0; i < DEG; i++)
     {
      int neighbor = W_index[uu][i];
			
		  //pthread_mutex_lock(&locks[neighbor]);   //THESE LOCKS ARE CAUSING SEG FAULTS IN SOFTBOUND

			D[W_index[uu][i]]++;
		  Q[W_index[uu][i]]=0;
	  
		  //pthread_mutex_unlock(&locks[neighbor]);
     }
	}
	
	// pthread_barrier_wait(&barrier);
	barrier_wait();
	
	
	//each thread adds its triangles into an array
	for(int uu=start;uu<stop;uu++)
	{
    unsigned int ret = -1;
	  ret = D[uu]/3;
	  D[uu]=ret;
	  if(D[uu]>=1)
	  {
		  Total_tid[tid] = Total_tid[tid]+D[uu];
	  }
	}
	printf("\n %d %d",tid, Total_tid[tid]);
	
	// pthread_barrier_wait(&barrier);
	barrier_wait();
	
	//finally the master thread adds up all the triangles from all the threads
	if(tid==0)
	{
	  for(int i=0;i<P;i++)
	  {
        Total = Total + Total_tid[i];
	  }
	}
   
   // pthread_barrier_wait(&barrier);
	barrier_wait();
   return 0;
}

int main(int argc, char *argv[])
{

   int N1 = 16384;//atoi(argv[1]);           //input arguments
   int DEG1 = 16;//atoi(argv[2]);
   int P1 = atoi(argv[1]);    //8

     // N=N1;
	 // DEG=DEG1;
	 P=P1;
     PMAX = P;

   
   // D = (int*) malloc((N)*sizeof(int));            //memory allocations
    int D1[N];
   if(!D) {D = D1;}
   // Q = (int*) malloc((N)*sizeof(int));
   int Q1[N];
   if(!Q) {Q = Q1;}
   // W = (int**)malloc((N)*sizeof(int*));
   // W_index = (int**)malloc((N)*sizeof(int*));
   //locks = (pthread_mutex_t*) malloc(N*(sizeof(pthread_mutex_t)));
   locks = new pthread_mutex_t[N];

   // for (int i = 0; i < N; i++)
   // {
	  // W[i] = (int*) malloc((DEG)*sizeof(int));
	  // W_index[i] = (int*) malloc((DEG)*sizeof(int));
   // }

pthread_mutex_init(&id, NULL);                         //barrier and lock inits
pthread_barrier_init(&barrier, NULL, P);
    for(int i=0; i<N; i++)
		pthread_mutex_init(&locks[i], NULL);

//fill D and Q
for(long i = 0; i < N; i++)
{
    D[i] = 0;
    Q[i] = 1;
}
		
//Now fill the graph with synthetic stuff

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
  
          // // Populate Cost Array
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

   
   
   	clock_gettime(CLOCK_REALTIME, &requestEnd);
	  double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;
	  printf( "\nTime:%lf\n", accum );   //prints time

// Disable performance and energy models
   CarbonDisableModels();
   
   printf("Triangle Counting: Done.\n");
   // printf("\nTriangles:%lld",Total);
   
   // Free memory

   //CarbonStopSim();
   return 0;
}
