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
  // int**     W_index;
  int       tid;
  int       P;
  int       N;
  int       DEG;
  pthread_barrier_t* barrier;
} thread_arg_t;

// HERERERERERE //////
// int start = 64;
int next_source = -1;
thread_arg_t thread_arg[1024];
pthread_t thread_handle[1024];

pthread_barrier_t  barrier;       //barrier to synchronize threads
pthread_mutex_t id;               //single lock for 1 variable
pthread_mutex_t* locks;   			//locks for each vertex in the graph
int tid_g = -1;                     //to increment global thread IDs later

// int* D;              
int* Q;
// int** W_index;                    //graph connections (adjacency list)
long W_index[262144][32] = {0};   //graph (adjacency list)  
// int Total_tid[1024] = {0};        //triangles for each thread

long long Total = 0;
// int *edges;                       //deg of a given vertex
int *exist;                      //whether vertex in graph
int N=262144;                          //global, otherwise softbound complains
int DEG=32;
int P=0;
// long long Total = 0;              //Total triangles (a large number so long long)


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
  long v = 0;
  
  
  // pthread_mutex_lock(&id);
	// tid_g++;
  // int tid = tid_g;
  
  // pthread_mutex_unlock(&id);
  int tid=0;
  tid = __sync_fetch_and_add(&tid_g,1);
  tid++;
    
    
	// pthread_barrier_wait(&barrier);
    
    // long ptr = 0;
    long ptr = 1;
   //divide graph among threads
	long start = tid    *  (N) / (P);
	long stop = (tid+1) *  (N) / (P);
	int node = 0;
    
    long stack[N/P];
    // stack = (int*) malloc((N/P)*sizeof(int));
    stack[0] = start;
    // printf("\n %d %ld == %ld",tid,stack[ptr],start);
    
	
   // pthread_barrier_wait(&barrier);
	 printf("\n %d %ld %ld",tid,start,stop);
     barrier_wait();
     
	   for(long vv=start;vv<stop;vv++)
      {
          
				if (ptr > 0)
                    ptr--;
                if (ptr >= (N/P))
                    break;
                
				v = stack[ptr];
                // printf(" %ld \n",stack[ptr]);
                // printf("\n %d v: %ld",tid, v);

         if(exist[v]==0)
             // printf("\n exist check %d v: %ld",tid, v);
             continue;
         // }                 
		int a;
         do{
             a=1;
             int *p;
              p = &Q[v];
             a = __sync_val_compare_and_swap(p,0,1);
            }while(!a);

         for(int i = 0; i < DEG; i++)
         {   
            long neighbor = W_index[v][i];
            stack[ptr] = neighbor;
            // printf("\n %d Using ptr %ld st %ld",tid,ptr, stack[ptr]);
            if(ptr < (N/P)-1)
             ptr++;
         }
      }
   // //pthread_barrier_wait(arg->barrier_total);
   printf("\n BF %d",tid);
   barrier_wait();
	return 0;
}

int main(int argc, char *argv[])
{

   int N1 = 16384; //atoi(argv[2]);  //1024;         //input arguments
   int DEG1 = 32; //atoi(argv[3]); //16;
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
  
   // // D = (int*) malloc((N)*sizeof(int));            //memory allocations
   // Q = (int*) malloc((N)*sizeof(int));
   int Q1[N];
   if(!Q) {Q = Q1;}
   // exist = (int*) malloc((N)*sizeof(int));
   int exist1[N];
   if(!exist) {exist = exist1;}
   // edges = (int*) malloc((N)*sizeof(int));
   // int edges1[N];
   // if(!edges) {edges = edges1;}
   // W_index = (int**)malloc((N)*sizeof(int*));
   locks = new pthread_mutex_t[N];
   


pthread_mutex_init(&id, NULL);                         //barrier and lock inits
pthread_barrier_init(&barrier, NULL, P);
    for(long i=0; i<N; i++)
		pthread_mutex_init(&locks[i], NULL);

// //fill D and Q
for(long i = 0; i < N; i++)
{
    // edges[i] = DEG;
    exist[i] = 1;
    Q[i] = 1;
}
// // Q[0] = 0;
	
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
								// long neighbor = i+j;
                                 W_index[i][j] = rand()%(N);			
                                
                                 /*									
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
                                   */
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

   // /*for (int i = 1; i < P; i++)
   // {
      // pthread_join(thread_handles[i], NULL);
   // }*/

   	clock_gettime(CLOCK_REALTIME, &requestEnd);
	  double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec ) / BILLION;
	  printf( "\nTime:%lf\n", accum );   //prints time

	   // Disable performance and energy models
   CarbonDisableModels();

   printf("Depth First Search: Done.\n");
   
   // Free memory

   //CarbonStopSim();
   return 0;
}
