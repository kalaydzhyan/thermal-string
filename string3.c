#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#define dim         3	   		     	 // Dimensionality of the space
#define MAX_LENGTH  100000		      	 // Maximal length of the string,
						 // otherwise use malloc()
#define N_UPDATES   300				 // Number of update cicles over the string length
#define NTHERM      60				 // Number of thermalization steps
#define nmeas	    10 				 // Interval between outputs
#define TENSION     0.42*0.42			 // String tension (in GeV^2)
#define a           2.73 			 // Lattice spacing (in GeV^{-1})
#define Elink       TENSION*a			 // Energy per link

float g = 1.0;					 // Gravitational constant
int   length = 0;				 // String length
float temperature;
float old_energy;				 // Current energy of the string
float new_energy;				 // ... and the new one after an update
float A, sigma;					 // Max value of the temperature and the size of the distribution
float pt0[dim] = {5.0, 5.0, 0.0};		 // Center of the temperature distribution

//////// Shift of the tail of the array ////////
void shift_array(int *array, char lr, int begin_position, int offset)
{
  int k;					 // Loop variable
  
  switch(lr)
  {
  case 'r':					 // Shift to the right
    if ((length + offset)*dim > MAX_LENGTH)
        { 
          printf("Can't shift, the string is too long!\n");
        };
        
    for (k = (length)*dim - 1; k >= begin_position*dim; k--) 
    {   
        *(array + k + offset*dim) = *(array + k);
    };
    length += offset;
    break;
  
  case 'l':					// Shift to the left
    if ((length-offset) < 0) printf("Too long shift to the left!\n");
    for (k = begin_position*dim; k < length*dim; k++)
    {
        *(array + k  - offset*dim) = *(array + k);
    };
    length -= offset;
    break;
  
  default:
    printf("Wrong parameter of array shifting!\n");
    break;
  }
} 

//////// Check if the point "stpoint" is occupied by string ////////
bool occupied(int array[][dim], int stpoint[dim])
{
  int k, d;					// Loop variables

  for(k = 0; k < length; k++)
  {
    if (memcmp(*(array + k), stpoint, dim*sizeof(int)) == 0) return true;
  }

  return false;
}

//////// Distance between two points ////////
double distance(int array[][dim], int p1, int p2)
{
  int d;					
  double dist = 0.0;				 
    
  for(d=0; d < dim; d++)
  {
    dist += pow(*(*(array + p1) + d) - *(*(array + p2) + d), 2);
  };

  return sqrt(dist);
}

//////// Gravitational energy of the string ////////
double gravit_energy(int array[][dim])
{
  int p, i;					 // Loop variables
  double energy = 0.0;

  for(p = 0; p < length; p++)
  {
    for(i = 1; i < length; i++)
    energy += -g*pow(TENSION*a, 2)/pow(a*distance(array, p, (p+i) % length), dim-2);
  }
  
  return energy / 2.0;                           // 2.0 because we count each interaction twice
}

// Space-dependent temperature T = A*exp(-(x-x0)^2 / sigma)
float temp(int pt[dim])
{
  float dist = 0;
  int d;
  
  for(d=0; d < dim; d++)
  {
    dist += pow((*(pt + d) - pt0[d]), 2);
  }
  
  return A*exp(-dist /(2*sigma*sigma));
}

//////// Gravitational energy of one node ////////
double gravit_energy_one_node(int array[][dim], int p1)
{
  int p, i;					 // Loop variables
  double energy = 0.0;

  for(i = 1; i < length; i++)
  energy += -g*pow(TENSION*a, 2)/pow(a*distance(array, p1, (p1+i) % length), dim-2);
  
  return energy;

}

//////// Gravitational energy between two nodes ////////
double gravit_energy_two_nodes(int array[][dim], int p1, int p2)
{
  return -g*pow(TENSION*a, 2)/pow(a*distance(array, p1, p2), dim-2);
}

//////// Average distance between points = radius of the string ball
double average_distance(int array[][dim])
{
  int p, i;					 // Loop variables
  double dist = 0.0;

  for(p = 0; p < length; p++)
  {
    for(i = 1; i < length; i++)
    dist += pow(distance(array, p, (p+i) % length),2);
  }
  
  return a*sqrt(dist / 2.0)/length;                // 2.0 because we count each distance twice
}

//////// Metropolis algorithm ////////

bool metropolis(float T, float old_E, float new_E)
{

  float p; 						// probability of the update
  float energy_difference = new_E - old_E;
  
  if(energy_difference < 0) p = 1.0;
  else p = exp(-energy_difference / T);
  
  return (drand48() <= p);
}

//////// Update of type 1 ///////
// This function also returns a bool = true if the update was performed 
bool update1(int array[][dim], int curpoint)
{
  int  d, d1, k;				// Loop variables 
  bool flags[2*dim];				// Allowed staples (x+, x-, y+, y- ...)
  int  staples[2*dim][2][dim];			// Staples themselves (pairs of points)
  int  orientation, increment;			// Orientation of the link
  int  n_staples = 2*(dim - 1);			// Number of allowed staples
  
  for(k = 0; k < 2*dim; k++) flags[k] = true;
  
  // Here we understand the orientation of the link
  for(d = 0; d < dim; d++)
    {
      increment =(*(*(array + curpoint) + d) - *(*(array + curpoint + 1) + d));
      if(abs(increment) == 1)
      {
        orientation = d;
        flags[2*d] = flags[2*d + 1] = false;
      };
    };
    
  // Generating staples
  for(d = 0; d < dim; d++)
  {
    for(d1 = 0; d1 < dim; d1++)
    {
      if(d1 == d)
      {
        staples[2*d][0][d1] = *(*(array + curpoint) + d1) + 1;
        staples[2*d][1][d1] = *(*(array + curpoint + 1) + d1) + 1;
        staples[2*d + 1][0][d1] = *(*(array + curpoint) + d1) - 1;
        staples[2*d + 1][1][d1] = *(*(array + curpoint + 1) + d1) - 1;
      }
      else
      {
        staples[2*d][0][d1] = *(*(array + curpoint) + d1);
        staples[2*d][1][d1] = *(*(array + curpoint + 1) + d1);
        staples[2*d + 1][0][d1] = *(*(array + curpoint) + d1);
        staples[2*d + 1][1][d1] = *(*(array + curpoint + 1) + d1);
      };
    }    
  };
  
  // Identify allowed staples
  for(k = 0; k <= 1; k++)
  for(d = 0; d < dim; d++)
  {
    if((occupied(array, &staples[2*d + k][0][0]) || occupied(array, &staples[2*d + k][1][0])) && (d != orientation))
      {
        flags[2*d + k] = false;
        n_staples--;
      }
  };
  
  if(n_staples == 0) return false;
  
  // Update. Adding a random staple.
  
  int staple_number = rand() % n_staples;
  int staple = -1;

  while(staple_number >=0)
  {
    staple++;
    if(flags[staple]) staple_number--;
  }

  shift_array(*array, 'r', curpoint + 1, 2);
  for(k = 0; k <= 1; k++)
  for(d = 0; d < dim; d++)
      *(*(array + curpoint + k + 1) + d) = staples[staple][k][d];

  new_energy = old_energy + gravit_energy_one_node(array, curpoint + 1) + gravit_energy_one_node(array, curpoint + 2)
               - gravit_energy_two_nodes(array, curpoint + 1, curpoint + 2) + 2*Elink;
    
  return metropolis(temperature, old_energy, new_energy);
}

//////// Update of type 2 ///////
// This function also returns a bool = true if the update was performed 
bool update2(int array[][dim], int curpoint)
{
  int d;				 	// Loop variable 
  bool flag = true;
  float median[dim] = { 0.0 };			
  int moved_point[dim];				// New position of the second point in the corner

  for(d = 0; d < dim; d++) 
    median[d] = (*(*(array + curpoint) + d) + *(*(array + curpoint + 2) + d))/2.0;
  
  // Pattern recognition: if 3 points do not lie on a straight line, then they form a corner 
  // flag = true if they are on a straight line
  for(d = 0; d < dim; d++)
    flag &= (median[d] == (float)*(*(array + curpoint + 1) + d));

  if(!flag) 
    {
      for(d = 0; d < dim; d++)
        moved_point[d] =(int)(2*median[d] - *(*(array + curpoint + 1) + d));
        
      if(!occupied(array, moved_point))
      {
        new_energy = old_energy - gravit_energy_one_node(array, curpoint + 1);
	 for(d = 0; d < dim; d++)  *(*(array + curpoint + 1) + d) = moved_point[d];
	 new_energy += gravit_energy_one_node(array, curpoint + 1);
        return metropolis(temperature, old_energy, new_energy);
      }
    };

  return false;
}

//////// Update of type 3 ///////
// This function also returns a bool = true if the update was performed 
bool update3(int array[][dim], int curpoint)
{

  // Pattern recognition: if 4th point lies in the neighborhood of the 1st, 
  // then they form a staple 
  if(distance(array, curpoint, curpoint + 3) == 1.0) 
    {
      new_energy = old_energy - gravit_energy_one_node(array, curpoint + 1) - gravit_energy_one_node(array, curpoint + 2) 
  		   + gravit_energy_two_nodes(array, curpoint + 1, curpoint + 2) - 2*Elink;
      shift_array(*array, 'l', curpoint + 3, 2); // Delete 2 points in between
      return metropolis(temperature, old_energy, new_energy);
    };

  return false;
}

//////////// Initial string shape from a file. TODO: arbitrary dimensions

void string_init(int array[][dim])
{
  FILE *init = fopen("./string.ini", "r");
  
  int test[dim];
  int status = 1;
  int d;
  
  if (init == NULL) 
  {
        fprintf(stderr, "Can't open input file!\n");
  }

  while(status==1)
  {
	for(d = 0; d < dim; d++)
	{
	  status = fscanf(init, "%i", &array[length][d]);
	}

    length++;
  }
  
  length--;
  
  fclose(init);
}

////////                            Main routine              /////////
int main(int argc, char *argv[])
{
  int i, j, i_update;				 // Loop variables
  int string[MAX_LENGTH][dim];			 // String itself
  int string_old[MAX_LENGTH][dim]; 		 // Here we save an old string configuration
  int length_old;  

  string_init(string);				 // Initialization of the string shape from a file

  old_energy = TENSION*a*length + gravit_energy(string);
  
  srand(time(NULL)); 			   	 // Initializing the randomizer
  srand48(time(NULL)); 			   	 // Initializing the randomizer

  // TEST: 
  // shift_array(&string[0][0], 'r' , 2 ,2);
  // int pt[3] = {4, 3, 0};
  // if(occupied(string, pt)) printf ("occupied!\n");
  // update2(string, 6);
 
  // int pnt[1][1][3]={{{4, 3, 0}}};
  // int pnt_copy[3];
  // memcpy(&pnt_copy, &pnt[0], dim*sizeof(int));
  // printf("%i\n", pnt_copy[0]);
  // if(occupied(string, &pnt[0][0][0])) printf("ok\n");

  // Test of Metropolis
  // for(i=0; i < 100; i++) printf(metropolis(temperature, 1.0, 1.1) ? "+" : "-");
  // printf("\n");
    

  //  printf("Enter the temperature (e.g. 1.0): ");
  //  scanf("%f", &temperature);

  if(argc != 4)
  {
    printf("This version of the program has T=A*exp(-(x-x0)^2)/(2*sigma^2) is physical units.\n");
    printf("Enter A (e.g. 1.0): ");
    scanf("%f", &A);
    printf("Enter sigma (e.g. 1.0): ");
    scanf("%f", &sigma);
    printf("Enter the gravitational constant (e.g. 1.0): ");
    scanf("%f", &g);
  }
  else
  {
    sscanf(argv[1], "%f",&A);
    sscanf(argv[2], "%f",&sigma);
    sscanf(argv[3], "%f",&g);
  };
  
  //////// Output initialization ////////

  char name[30];

  snprintf(name,sizeof(name),"./string_A%1.4f_g%1.2f.dat",A,g);
  FILE *fp = fopen(name, "w");

  if (fp == NULL) 
  {
        fprintf(stderr, "Can't open output file!\n");
  }
    
  int   curpoint[dim];				  // coordinates of the current point
  float update_en = old_energy;			  // energy before an update cycle

// the main updates start

  for(i_update = 0; i_update < N_UPDATES; i_update++)
  {  
    for(i = 0; i < length - 1; i++)
    { 
      memcpy(&curpoint[0], &string[i][0], dim*sizeof(int));
      temperature = temp(curpoint);
      memcpy(&string_old[0][0], &string[0][0], length*dim*sizeof(int));
      length_old = length;
      //if update is ok, then jump a bit forward, otherwise restore the old configuration
      if(update1(string, i) && length >= 4) 
      {
	 i+=2;
	 old_energy = new_energy;
      }
      else 
      {
        memcpy(&string[0][0], &string_old[0][0], length_old*dim*sizeof(int));
        length = length_old;
      }
     };
 
    for(i = 0; i < length - 2; i++)
    {
      memcpy(&curpoint[0], &string[i][0], dim*sizeof(int));
      temperature = temp(curpoint);
      memcpy(&string_old[0][0], &string[0][0], length*dim*sizeof(int));
      length_old = length;
      //if update is ok, then jump a bit forward, otherwise restore the old configuration
      if(update2(string, i) && length >= 4) 
      {
	 i++;
	 old_energy = new_energy;
      }
      else 
      {
        memcpy(&string[0][0], &string_old[0][0], length_old*dim*sizeof(int));
        length = length_old;
      }

    };

    for(i = 0; i < length - 3; i++)
    {
      memcpy(&curpoint[0], &string[i][0], dim*sizeof(int));
      temperature = temp(curpoint);
      memcpy(&string_old[0][0], &string[0][0], length*dim*sizeof(int));
      length_old = length;
      //if update is ok, then update energy, otherwise restore the old configuration
      if(update3(string, i) && length >= 4)
      {
	 old_energy = new_energy;
      }
      else 
      {
        memcpy(&string[0][0], &string_old[0][0], length_old*dim*sizeof(int));
        length = length_old;
      }
    };

      ////// Output begins /////
    if(i_update % nmeas==0 && i_update >= NTHERM)
      {
          fprintf(fp, "%i\n", length);			 // Write the string length

          for(i=0; i < length; i++)			
          {
              for(j=0; j < dim; j++)
              {
                  fprintf(fp, "%i ", string[i % length][j]);
              }
              fprintf(fp, "\n");
          };
      }
      
    printf("Energy increment at the step %2.i : %6.2f, string length: %i\n", i_update + 1, old_energy - update_en, length);
    update_en = old_energy;
      
  }
 // end of updates
  
  fclose(fp);
    
    //  system("./string.plt");				// Export the string to an image file

  printf("\nRadius of the ball: %f \n", average_distance(string));
  printf("Total energy (gravit. energy): %f  (%f) \n", TENSION*a*length + gravit_energy(string), gravit_energy(string));

  return 0;
}
