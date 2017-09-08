#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#define dim	    3	   		     	 // Dimensionality of the space
#define MAX_LENGTH  100000		      	 // Maximal length of the string,
						 // otherwise use malloc()
#define N_UPDATES   15				 // Number of update cicles over the string length

float g = 1.0;					 // Gravitational constant
int length = 10;				 // String length
float temperature = 1;
float old_energy;				 // Current energy of the string
float new_energy;				 // ... and the new one after an update

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
    dist += pow(*(*(array + p1) + d) - *(*(array + p2) + d),2);
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
    energy += -g/pow(distance(array, p, (p+i) % length), dim-2);
  }
  
  return energy / 2.0;                           // 2.0 because we count each interaction twice
}

//////// Gravitational energy of one node ////////
double gravit_energy_one_node(int array[][dim], int p1)
{
  int p, i;					 // Loop variables
  double energy = 0.0;

    for(i = 1; i < length; i++)
    energy += -g/pow(distance(array, p1, (p1+i) % length), dim-2);
  
  return energy;

}

//////// Gravitational energy between two nodes ////////
double gravit_energy_two_nodes(int array[][dim], int p1, int p2)
{
  return -g/pow(distance(array, p1, p2), dim-2);
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
  
  return sqrt(dist / 2.0)/length;                // 2.0 because we count each distance twice
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
  int d, d1, k;				 	// Loop variables 
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
               - gravit_energy_two_nodes(array, curpoint + 1, curpoint + 2) + 2;
    
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
        new_energy = old_energy - gravit_energy_one_node(array, curpoint);
	 for(d = 0; d < dim; d++)  *(*(array + curpoint + 1) + d) = moved_point[d];
	 new_energy += gravit_energy_one_node(array, curpoint);
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
  		   + gravit_energy_two_nodes(array, curpoint + 1, curpoint + 2) - 2;
      shift_array(*array, 'l', curpoint + 3, 2); // Delete 2 points in between
      return metropolis(temperature, old_energy, new_energy);
    };
  // TODO: Maybe also the last (closing) part of the sting can be also updated.
  return false;
}

//////// Main routine /////////
int main (void)
{
  int i, j, i_update;				 // Loop variables
  int string[MAX_LENGTH][dim] = {		 // Initial shape of the string, cold start
			        {4, 3, 0},
				{4, 4, 0},
		  		{4, 5, 0},
		  		{4, 6, 0},
		  		{5, 6, 0},
		  		{6, 6, 0},
		  		{6, 5, 0},
		  		{6, 4, 0},
		  		{6, 3, 0}, 
				{5, 3, 0}
				}; 
  
  int string_old[MAX_LENGTH][dim]; 		 // Here we save an old string configuration
  int length_old;  

  old_energy = length + gravit_energy(string);
  
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
    

  printf("Enter the temperature (e.g. 1.0): ");
  scanf("%f", &temperature);
  printf("Enter the gravitational constant (e.g. 1.0): ");
  scanf("%f", &g);

  // TEST: energy increment after an update loop
  int update_en = old_energy;

  for(i_update = 0; i_update < N_UPDATES; i_update++)
  {  
    for(i = 0; i < length - 1; i++)
    {
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

    printf("Energy increment at the step %2.i : %6.2f, string length: %i\n", i_update + 1, new_energy - update_en, length);
  }
  
  //////// BEGIN Output
  FILE *fp = fopen("./string.dat", "w");
  
  if (fp == NULL) {
		    fprintf(stderr, "Can't open output file!\n");
	          }
  
  for(i=0; i <= length; i++)			 // The first points is written twice
  {
     for(j=0; j < dim; j++)
     {
       fprintf(fp, "%i ", string[i % length][j]);
     }
     
     fprintf(fp, "\n");
  };
  
  fclose(fp);
  
//  system("./string.plt");			// Export the string to an image file

  //////// END Output

  printf("\nRadius of the ball: %f \n", average_distance(string));
  printf("Total energy (gravit. energy): %f  (%f) \n", length + gravit_energy(string), gravit_energy(string));

  return 0;
}
