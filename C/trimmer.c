#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


/* Maximum allowed number of points*/
#define MAXN 100
/* we use a minimum number of allowed points 
   in order to have enough points for our demo*/
#define DEFAULTMINPO 30
/* Randomly generated data involve at most 10 segments, each with Gaussian error whose
 standard deviance which is at most 0.01*/
#define SEGMENTS 5
#define MAXMEAN 10
#define STD 0.01
#define PI 3.14159265358979323846
/* trimmer.c 
   This program implements the O(n^2) dynamic programming algorithm of Tsourakakis et al. (ACM JEA '11)
   and some demos. 
   author: Charalampos E. Tsourakakis  

   Paper: Approximation Algorithms for Speeding up Dynamic Programming and Denoising aCGH data
   C. Tsourakakis, R. Peng, M. Tsiarli, G. Miller, R. Schwartz
   Howto compile: gcc -o trimmer trimmer.c and if it doesn't work (because of the math.h library) 
                  compile using gcc -lm -Wall  -o trimmer trimmer.c   
	Tested in Ubuntu  10.x and Windows 7
*/



/* Number of points*/
int n=-1; 

/* Normalization Constant*/
double C=1; 

/* Input signal */
double aCGH[MAXN]; 

/*Means of randomly generated segments*/ 
int rMeans[SEGMENTS];

/*first and second moments*/ 
double firstmoment[MAXN+1]; /*think of first index as -1*/
double secondmoment[MAXN+1];

/*breakpoints*/
int breakpoints[MAXN];

/* OPT_i values, see paper for their meaning */ 
double OPT[1+MAXN];


/* Fitted Signal, i.e., FIT[i] <-> aCGH[i], but F tends to have piecewise constant segments in it*/ 
double FIT[MAXN]; 

/* absolute*/ 
int abs(int);
/* This function reads from the standard input the number of points and the normalization 
constant C and subsequently reads the values of the n points from input. */
void input(); 
/* This function prints the aCGH signal.*/
void print_acgh();
/* This function generates a random piecewise constant signal. */
void generate_random_acgh(int nopoints); 
/* This function generates a random piecewise constant signal with small noise along each segment. */
void generate_random_noisy_acgh(int nopoints,double std);
/* These functions print a vector*/
void vprint(double v[],int size);
void v2print(int v[],int size);
/* Dynamic programming solution, see Vanilla DP algorithm in Tsourakakis et al. ACM JEA '11*/ 
void vanillaDP();
/* square function */
double square(double x);
/* print fitted Segments from index 0 until index ind*/
void printsegments(int ind);
/* generate fitted signal */
void generateF();
/* menu */
void menu();
/* get choice */ 
char get_choice();
/* generate N(0,1) using Box Muller */
double generatenoise(double mu,double sigma);

int main()
{
   int choice, nopoints;
   menu();  
   while( (choice = get_choice())!='q')
   { 
    switch(choice){
      case 'a':  input();
                 vanillaDP(); 
                 generateF();
				 if( n <=100){
				    printsegments(n-1);                  
					printf("Input Signal\n");
					vprint(aCGH,n); 
					printf("Output Signal\n");
					vprint(FIT,n);
                }
                else{
				    printf("Signal is too large. Printing the first 20 points \n");  
					printf("Input Signal\n");
					vprint(aCGH,20); 
					printf("Output Signal\n");
					vprint(FIT,20);
				}
                break;
      case 'b': printf("Enter #points: ");
                scanf("%d", &nopoints);
				printf("\n");
	            generate_random_acgh(nopoints);
                vanillaDP();   
                generateF();
				if( n <=100){
				        printsegments(n-1);
				  	printf("Input Signal\n");
					vprint(aCGH,n); 
					printf("Output Signal\n");
					vprint(FIT,n);
                }
                else{
				    printf("Signal is too large. Printing the first 20 points \n");  
					printf("Input Signal\n");
					vprint(aCGH,20); 
					printf("Output Signal\n");
					vprint(FIT,20);
				}				
                break;
      case 'c': printf("Enter #points: ");
                scanf("%d", &nopoints);
				printf("\n");
	            generate_random_noisy_acgh(nopoints,STD);
                vanillaDP();   
                generateF(); 
				if( n <=100){
				    printsegments(n-1);
					printf("Input Signal\n");
					vprint(aCGH,n); 
					printf("Output Signal\n");
					vprint(FIT,n);
                }
                else{
				    printf("Signal is too large. Printing the first 20 points \n");  
					printf("Input Signal\n");
					vprint(aCGH,20); 
					printf("Output Signal\n");
					vprint(FIT,20);
				}					
                break;
      case 'q': printf("Bye bye \n"); 
                break;
      }
   }
   return 1; 
} 




void input() {
  int i;
  printf("Enter #points n and the normalization constant C\n");
  scanf("%d%lf", &n, &C);
  printf("C=%lf n=%d.\n Enter now the values of the n points\n",C,n);
  for(i = 0; i < n; ++i)
    scanf("%lf", &aCGH[i]);
}



void print_acgh()
{ 
   int i;
   if(n==-1) 
   {
     printf("aCGH signal not initialized yet \n"); 
     return; 
   }
   printf("[");
   for(i=0;i<n;i++)
   {
       printf("%lf ",aCGH[i]); 
   } 
   printf("]\n"); 
  
}


void generate_random_acgh(int nopoints)
{    
  int j,i; 
  if( nopoints >MAXN)
  { 
     printf("Number of points asked %d is greater than maximum allowed %d\n",nopoints,MAXN);
     n=MAXN; 
  }
  else if( nopoints < 30)
  {
       printf("Number of points too small, i.e., %d less than 30. \n #POINTS automatically set to %d.\n", nopoints, DEFAULTMINPO);
       n=DEFAULTMINPO; 
  }
  else
	   n = nopoints;
  srand( time(0) );
  /* number of segments */
  int seg = rand()%SEGMENTS+1; 
  printf("%d Segments shall be generated\n",seg);
  /* Generate #seg different means*/   
  for(i=0;i<seg;i++)
     rMeans[i] = rand()%MAXMEAN+1; 
  int points = (int)MAXN/seg; 
  int counter=0; 
  for(i=0;i<(seg-1);i++)
  {
     for(j=0;j<points;j++)
     {   
        aCGH[counter]= rMeans[i]; 
        counter+=1;
     }
  }
  /* the last shall become first, i.e., the last segments gets rest */
  for(i=counter;i<MAXN;i++)
        aCGH[i] = rMeans[seg-1]; 
} 

void generate_random_noisy_acgh(int nopoints,double std)
{    
  int j,i; 
  if( nopoints >MAXN)
  { 
     printf("%Number of points asked %d is greater than MAX allows %d\n",nopoints,MAXN);
     n=MAXN; 
  }
  else if( nopoints < 30)
  {
      printf("Number of points too small, i.e., %d less than 30. \n #POINTS automatically set to %d.\n", nopoints, DEFAULTMINPO);
       n=DEFAULTMINPO; 
  }
  else
	n = nopoints;
  srand( time(0) );
  /* number of segments */
  int seg = rand()%SEGMENTS+1; 
  printf("%d Segments shall be generated\n",seg);
  /* Generate #seg different means*/   
  for(i=0;i<seg;i++)
     rMeans[i] = rand()%MAXMEAN+1; 
  int points = (int)MAXN/seg; 
  int counter=0; 
  for(i=0;i<(seg-1);i++)
  {
     for(j=0;j<points;j++)
     {   
        aCGH[counter]= rMeans[i]+generatenoise(0,std); 
        counter+=1;
     }
  }
  /* the last segments gets rest */
  for(i=counter;i<MAXN;i++)
        aCGH[i] = rMeans[seg-1]+generatenoise(0,std); ; 
} 


void vprint(double v[],int size)
{ 
   int i;
   printf("[");
   for(i=0;i<size;i++)
   {
      printf("%lf ",v[i]);
   }
   printf("]\n"); 

}


void v2print(int v[],int size)
{ 
   int i;
   printf("[");
   for(i=0;i<size;i++)
   {
      printf("%d ",v[i]);
   }
   printf("]\n"); 

}

void calc_moments()
{
  
   int i;
   firstmoment[0]=0;
   secondmoment[0]=0;
   for(i=1;i<n+1;i++) 
   {
        firstmoment[i]=firstmoment[i-1]+aCGH[i-1];
        secondmoment[i]=secondmoment[i-1]+square(aCGH[i-1]);
   }  
   
  /* printf("Moments calculated\n");
   vprint(firstmoment,n+1);
   vprint(secondmoment,n+1);
   */
}




void vanillaDP()
{   
    int i,j; 
    double tmp;
    calc_moments();
    for(i =1; i<n+1;i++)
         OPT[i]=RAND_MAX;
    OPT[0]=0; 
    breakpoints[0]=0; 
    for(i=1;i<=n;i++) 
    {  
         for(j=0;j<i;j++)
         {  
               
              tmp = OPT[j]+(secondmoment[i]-secondmoment[j]) - square( (firstmoment[i]-firstmoment[j]) )/(i-j)+C;
              if(OPT[i]>tmp)
              { 
                   breakpoints[i] = j; 
                   OPT[i]=tmp;
				   
              }
		}   
         
    }     
}


double square(double x)
{ 
   return x*x; 
}


/* low has to be >=1 when this is being called, index 0 served as a dummy index,
   Function uses recursion to print them in correct order. */

void printsegments(int ind)
{
   
   if(!breakpoints[ind+1]) 
       printf("Interval from %d to %d\n", 0,ind); 
   else
   {
       printsegments(breakpoints[ind+1]-1);
       printf("Interval from %d to %d\n", breakpoints[ind+1],ind);    
   } 
}  


void generateF()
{
      int up=n; 
      int b=breakpoints[up]; 
      int i;
      while(b) 
      {
          for(i=b;i<up;i++) 
          { 
               FIT[i]=(firstmoment[up]-firstmoment[b])/(up-b); 
//               printf("firstmoment[%d]=%lf,firstmoment[%d]=%lf, FIT[i]=%lf\n",up,firstmoment[up],b,firstmoment[b],FIT[i]);
          }
          up = b; 
          b = breakpoints[up];
      }
     for(i=0;i<up;i++) 
         FIT[i]=(firstmoment[up]-firstmoment[b])/(up-b);   
     

}


void menu()
{ 
   printf("Welcome! \n");
   printf("[a]Enter input signal from keyboard\n[b]Generate Piecewise constant signal\n[c]Generate Piecewise constant signal with noise\n[q] quit \n");
}


char get_choice()
{
    int ch;
    printf("Enter the letter of your choice:\n");
    ch=getchar();
    while( (ch<'a' || ch>'c') && ch!='q')
    {  
        printf("[a]Enter input signal from keyboard\n[b]Generate Piecewise constant signal\n[c]Generate Piecewise constant signal with noise\n[q] quit \n");
        ch = getchar();
    }
    return ch;
}


int abs(int a)
{
   return (a<0?-a:a);
}


double generatenoise(double mu,double sigma)
{
    double dist = sqrt( -2.0 * log( (double)rand() / (double)RAND_MAX ) ) ;
    double angle = 2.0 * PI * (double)rand() / (double)RAND_MAX;
    return dist * sin(angle) * sigma + mu;
    return 1;
}      