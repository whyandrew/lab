/*
  Particle filters implementation for a simple robot.

  Your goal here is to implement a simple particle filter
  algorithm for robot localization.

  A map file in .ppm format is read from disk, the map
  contains empty spaces and walls.

  A simple robot is randomly placed on this map, and
  your task is to write a particle filter algorithm to
  allow the robot to determine its location with
  high certainty.

  You must complete all sections marked

  TO DO:

  NOTE: 2 keyboard controls are provided:

  q -> quit (exit program)
  r -> reset particle set during simulation

  Written by F.J.E. for CSC C85, May 2012. Updated Aug. 15, 2014

  This robot model inspired by Sebastian Thrun's
  model in CS373.
*/

#include "ParticleFilters.h"
/**********************************************************
 GLOBAL DATA
**********************************************************/
unsigned char *map;		// Input map
unsigned char *map_b;		// Temporary frame
struct particle *robot;		// Robot
struct particle *list;		// Particle list
int sx,sy;			// Size of the map image
char name[1024];		// Name of the map
int n_particles;		// Number of particles
int windowID;               	// Glut window ID (for display)
int Win[2];                 	// window (x,y) size
int RESETflag;			// RESET particles
double sigma = 20.0; // Sigma value for Gaussian PDR evaluation
/**********************************************************
 PROGRAM CODE
**********************************************************/
int main(int argc, char *argv[])
{
 /*
   Main function. Usage for this program:

   ParticleFilters map_name n_particles

   Where:
    map_name is the name of a .ppm file containing the map. The map
             should be BLACK on empty (free) space, and coloured
             wherever there are obstacles or walls. Anythin not
             black is an obstacle.

    n_particles is the number of particles to simulate in [100, 50000]

   Main loads the map image, initializes a robot at a random location
    in the map, and sets up the OpenGL stuff before entering the
    filtering loop.
 */

 if (argc!=3)
 {
  fprintf(stderr,"Wrong number of parameters. Usage: ParticleFilters map_name n_particles.\n");
  exit(0);
 }

 strcpy(&name[0],argv[1]);
 n_particles=atoi(argv[2]);

 if (n_particles<100||n_particles>50000)
 {
  fprintf(stderr,"Number of particles must be in [100, 50000]\n");
  exit(0);
 }

 fprintf(stderr,"Reading input map\n");
 map=readPPMimage(name,&sx, &sy);
 if (map==NULL)
 {
  fprintf(stderr,"Unable to open input map, or not a .ppm file\n");
  exit(0);
 }

 // Allocate memory for the temporary frame
 fprintf(stderr,"Allocating temp. frame\n");
 map_b=(unsigned char *)calloc(sx*sy*3,sizeof(unsigned char));
 if (map_b==NULL)
 {
  fprintf(stderr,"Out of memory allocating image data\n");
  free(map);
  exit(0);
 }

 srand48((long)time(NULL));		// Initialize random generator from timer
 fprintf(stderr,"Init robot...\n");
 robot=initRobot(map,sx,sy);
 if (robot==NULL)
 {
  fprintf(stderr,"Unable to initialize robot.\n");
  free(map);
  free(map_b);
  exit(0);
 }
 sonar_measurement(robot,map,sx,sy);	// Initial measurements...

 // Initialize particles at random locations
 fprintf(stderr,"Init particles...\n");
 list=NULL;
 initParticles();

 // Done, set up OpenGL and call particle filter loop
 fprintf(stderr,"Entering main loop...\n");
 Win[0]=800;
 Win[1]=800;
 glutInit(&argc, argv);
 initGlut(argv[0]);
 glutMainLoop();

 // This point is NEVER reached... memory leaks must be resolved by OpenGL main loop
 exit(0);

}

void initParticles(void)
{
 /*
   This function creates and returns a linked list of particles
   initialized with random locations (not over obstacles or walls)
   and random orientations.

   There is a utility function to help you find whether a particle
   is on top of a wall.

   Use the global pointer 'list' to keep trak of the *HEAD* of the
   linked list.

   Probabilities should be uniform for the initial set.
 */
 
 particle *prev, *curr; // Pointers to the previous particle in the list (if any),
                        // and the most recently created one

 double initial_prob = 1.0 / n_particles; // Initial probablitity for each particle
                                           // (distributed evenly for a total of 1.0)
  
 prev = list = NULL; // Initially the list has no head, and there is no previous particle

 for (int i=0; i<n_particles; i++)
 {
  curr = initRobot(map, sx, sy); // Try to create new particle to add to the list
  
  if (curr == NULL) // initRobot() returned a null pointer so we can't proceed
  {
   fprintf(stderr, "Init particles failed, exiting. (Out of memory?)\n");
   free(map);
   free(map_b);
   free(robot);
   deleteList(list);
   exit(0);
  }
  
  curr->prob = initial_prob;     // Set particle probablility
   
  if (prev==NULL) // First iteration, so the new particle will be head of the list
  {
   list = curr;
  }
  else            // Not the first particle in the list, so link to the new particle 
  {               // from the previous one
   prev->next = curr;
  }
  prev = curr; // The particle created this time will be the "previous" one in
 }             // the next loop iteration
}

void computeLikelihood(struct particle *p, struct particle *rob, double noise_sigma)
{
 /*
   This function computes the likelihood of a particle p given the sensor
   readings from robot 'robot'.

   Both particles and robot have a measurements array 'measureD' with 16
   entries corresponding to 16 sonar slices.

   The likelihood of a particle depends on how closely its measureD
   values match the values the robot is 'observing' around it. Of
   course, the measureD values for a particle come from the input
   map directly, and are not produced by a simulated sonar.

   Assume each sonar-slice measurement is independent, and if

   error_i = (p->measureD[i])-(sonar->measureD[i])

   is the error for slice i, the probability of observing such an
   error is given by a Gaussian distribution with sigma=20 (the
   noise standard deviation hardcoded into the robot's sonar).

   You may want to check your numbers are not all going to zero...

   This function updates the likelihood for the particle in
   p->prob
 */

 /****************************************************************
 // TO DO: Complete this function to calculate the particle's
 //        likelihood given the robot's measurements
 ****************************************************************/
 double error_i; // Error value for sensor direction i
 for (int i=0; i<16; i++)
 {
  error_i = (p->measureD[i])-(rob->measureD[i]);
  p->prob *= GaussEval(error_i, sigma);
 }

}

void move_bounce(struct particle *p, double dist)
{
 /* Move the particle forward by distance dist using the provided function
 move(), and if this movement places the particle in contact with a solid object, 
 move the particle in some random direction that takes the particle into an open
 space again.
 */
 int hit_wall = 0;       // 1 if a wall has been hit, 0 otherwise
 int i = 0;              // Count of tries for debugging purposes
 struct particle *copy;  // Copy of the particle to move
 
 copy = (particle*)calloc(1, sizeof(particle));
 do // Runs once if particle doesn't hit a wall, otherwise however many
 {  // times it takes to randomly find a direction with no hit
  *copy = *p;
  if (hit_wall) // Wall has been hit, only executes after first loop iteration
  {
   printf("Hit obstacle. Try#%d, dist %f\n", ++i, dist);
   copy->theta = 12.0 * drand48(); // Assign a random direction
  }
  move(copy, dist);
  hit_wall = hit(copy, map, sx, sy);
 } while (hit_wall);
 *p = *copy;
 free(copy);
// move(p, dist);
// if (hit(p, map, sx, sy))
// {
//  copy = (particle*)calloc(1, sizeof(particle));
//  
//  int i = 0;
// *copy = *p;
//  while (hit(copy, map, sx, sy))  // Loop until it isn't hitting anything
//  {
//   *copy = *p;                     // Make a fresh copy
//   //copy->theta = 12.0 * drand48(); // Assign a random direction
//  
//   // Reverse direction w/ slight random adjustment, and correct if we go over 12
//   copy->theta += 6.0 + 6.0 * (drand48() - 0.5);
//   copy->theta = (copy->theta <= 12) ? copy->theta : copy->theta - 12.0;
//   
//   printf("Hit obstacle. Try#%d: direction %f\n", i++, copy->theta);
//   move(copy, dist * 5.0);               // Move forward in the new direction
                                        // (farther to increase chance of clearing it)
//  }

//  *p = *copy; // After completion of the loop, the copy is free of obstacles,
              // so we can copy it to p.
//  free(copy); // Release memory used for the copy
// }
}

void ParticleFilterLoop(void)
{
 /*
    Main loop of the particle filter
 */

  // OpenGL variables. Do not remove
  unsigned char *tmp;
  GLuint texture;
  static int first_frame=1;
  double max;
  struct particle *p,*pmax;
  char line[1024];
  // Add any local variables you need right below.
  
  int i; // Loop counter 
  double total_prob = 0.0;
  double move_distance = 10.0;
  particle *curr = list; // Current particle, initially the one at the head of the list
  
  double normalized_total = 0.0; // Stores running total of normalized particle weights
  double weights[n_particles];   // Stores successive running totals of 
                                  // normalized particle weights between 0 and 1. Used for
                                  // binary search, to find indices that correspond to
                                  // randomly chosen particles.
  particle *particles[n_particles]; // Array of pointers to particles, ordered as in the
                                    // linked list                  

  
  particle *prev;    // Previous particle while traversing original list (during resampling)
  particle *new_list = NULL; // Pointer to head of new resampled list
  double selection; // Holds random value in 0,1 used in selecting a particle to
                     // copy to the new list during resampling
  particle *curr_p_newlist; // Pointer to most newly created particle during resampling
  particle *prev_p_newlist; // Pointer to the particle preceding the most newly created
                            // one during resampling

  if (!first_frame)
  {
   // Step 1 - Move all particles a given distance forward (this will be in
   //          whatever direction the particle is currently looking).
   //          To avoid at this point the problem of 'navigating' the
   //          environment, any particle whose motion would result
   //          in hitting a wall should be bounced off into a random
   //          direction.
   //          Once the particle has been moved, we call ground_truth(p)
   //          for the particle. This tells us what we would
   //          expect to sense if the robot were actually at the particle's
   //          location.
   //
   //          Don't forget to move the robot the same distance!
   //
   move_bounce(robot, move_distance);
   //if (hit(robot, map, sx, sy)) {
   // robot->theta = (robot->theta + 6.0) % 12.0
   // move(robot, move_distance);
   //}
   while (curr != NULL)
   {
    move_bounce(curr, move_distance);
    ground_truth(curr, map, sx, sy);
    //if (hit(curr, map, sx, sy)) {
    // curr->theta = (curr->theta + 6.0) % 12.0
    // move(robot, move_distance);
    //}
    curr = curr->next;
   }
   /******************************************************************
   // TO DO: Complete Step 1 and test it
   //        You should see a moving robot and sonar figure with
   //        a set of moving particles.
   ******************************************************************/

   // Step 2 - The robot makes a measurement - use the sonar
   sonar_measurement(robot,map,sx,sy);

   // Step 3 - Compute the likelihood for particles based on the sensor
   //          measurement. See 'computeLikelihood()' and call it for
   //          each particle. Once you have a likelihood for every
   //          particle, turn it into a probability by ensuring that
   //          the sum of the likelihoods for all particles is 1.

   /*******************************************************************
   // TO DO: Complete Step 3 and test it
   //        You should see the brightness of particles change
   //        depending on how well they agree with the robot's
   //        sonar measurements. If all goes well, particles
   //        that agree with the robot's position/direction
   //        should be brightest.
   *******************************************************************/
   
   curr = list;
   while (curr != NULL)
   {
    computeLikelihood(curr, robot, sigma);
    total_prob += curr->prob;
    curr = curr->next;
   }
   
   i = 0;
   curr = list;
   while (curr != NULL)
   {
    particles[i] = curr;
    curr->prob /= total_prob;
    normalized_total += curr->prob;
    weights[i++] = normalized_total;
    curr = curr->next;
   }

   // Step 4 - Resample particle set based on the probabilities. The goal
   //          of this is to obtain a particle set that better reflect our
   //          current belief on the location and direction of motion
   //          for the robot. Particles that have higher probability will
   //          be selected more often than those with lower probability.
   //
   //          To do this: Create a separate (new) list of particles,
   //                      for each of 'n_particles' new particles,
   //                      randomly choose a particle from  the current
   //                      set with probability given by the particle
   //                      probabilities computed in Step 3.
   //                      Initialize the new particles (x,y,theta)
   //                      from the selected particle.
   //                      Note that particles in the current set that
   //                      have high probability may end up being
   //                      copied multiple times.
   //
   //                      Once you have a new list of particles, replace
   //                      the current list with the new one. Be sure
   //                      to release the memory for the current list
   //                      before you lose that pointer!
   //
   
   /*******************************************************************
   // TO DO: Complete and test Step 4
   //        You should see most particles disappear except for
   //        those that agree with the robot's measurements.
   //        Hopefully the largest cluster will be on and around
   //        the robot's actual location/direction.
   *******************************************************************/
   
   prev_p_newlist = NULL; // Pointer to current particle in list being created
   
   for (i=0; i<n_particles; i++)
   {
    total_prob = 0.0;      // Running total of particle probabilities
    selection = drand48(); // Random double in 0,1
    curr = list;           // Current particle
    
    // Do a binary search for the particle corresponding to selection,
    // based on a precalculated table of running total of particle probabilities.
    // The target particle is the one for which total_prob has the least value
    // not greater than selection. When the while loop terminates, upper_bound
    // is the index of the selected particle.
     
    // To do: replace with something more efficient? Maybe put precalculated
    // running totals into an array and do binary search?
    
    //int j = 0;
    //while(total_prob<selection && curr != NULL)
    //{
    // total_prob += curr->prob;
    // prev = curr;
    // curr = curr->next;
    // j++;
    //}
    int diff, guess;
    int lower_bound = 0;
    int upper_bound;
    upper_bound = n_particles-1;
    while ((diff = upper_bound - lower_bound) > 1)
    {
     guess = lower_bound + diff / 2;
     if (weights[guess] <= selection)
     {
      lower_bound = guess;
     }
     else
     {
      upper_bound = guess;
     }
     //printf("Lower: %d upper: %d selection: %f guess: %f guess-1: %f\n", lower_bound, upper_bound, selection, weights[guess], weights[guess-1]);
    }
    //printf("%f *%f* %f\n", particles[guess-1]->prob, particles[guess]->prob, particles[guess+1]->prob);
    //printf("lower %d upper %d guess %d\n", lower_bound, upper_bound, guess);
    
    
    
    // Allocate memory for new particle and check for failure
    if ((curr_p_newlist = (particle*)calloc(1, sizeof(particle))) == NULL)
    {
     fprintf(stderr, "Resample particles failed, exiting. (Out of memory?)\n");
     free(map);
     free(map_b);
     free(robot);
     deleteList(list);
     deleteList(new_list);
     exit(0);
    }
    // Copy selected particle to the new one
    *curr_p_newlist = *(particles[upper_bound]);//*prev;
    
    if (prev_p_newlist == NULL) // First particle in new list so we link it as head
    {
     new_list = curr_p_newlist;
    }
    else // Otherwise link the previous particle to the new one
    {
     prev_p_newlist->next = curr_p_newlist;
    }
    prev_p_newlist = curr_p_newlist; // Advance to next particle for next iteration
   }
   curr_p_newlist->next = NULL; // Make sure the final particle links to NULL,
                                // instead of whatever the particle it was copied
                                // from pointed to

   deleteList(list); // Free memory from original list
   list = new_list;  // Link "list" pointer to the new list
   
  }  // End if (!first_frame)

  /***************************************************
   OpenGL stuff
   You DO NOT need to read code below here. It only
   takes care of updating the screen.
  ***************************************************/
  if (RESETflag)	// If user pressed r, reset particles
  {
   deleteList(list);
   list=NULL;
   initParticles();
   RESETflag=0;
  }
  renderFrame(map,map_b,sx,sy,robot,list);

  // Clear the screen and depth buffers
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glEnable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);

  glGenTextures( 1, &texture);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glBindTexture( GL_TEXTURE_2D, texture);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

  glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, sx, sy, 0, GL_RGB, GL_UNSIGNED_BYTE, map_b);

  // Draw box bounding the viewing area
  glBegin (GL_QUADS);
  glTexCoord2f (0.0, 0.0);
  glVertex3f (0.0, 100.0, 0.0);
  glTexCoord2f (1.0, 0.0);
  glVertex3f (800.0, 100.0, 0.0);
  glTexCoord2f (1.0, 1.0);
  glVertex3f (800.0, 700.0, 0.0);
  glTexCoord2f (0.0, 1.0);
  glVertex3f (0.0, 700.0, 0.0);
  glEnd ();

  p=list;
  max=0;
  while (p!=NULL)
  {
   if (p->prob>max)
   {
    max=p->prob;
    pmax=p;
   }
   p=p->next;
  }

  if (!first_frame)
  {
   sprintf(&line[0],"X=%3.2f, Y=%3.2f, th=%3.2f, EstX=%3.2f, EstY=%3.2f, Est_th=%3.2f, Error=%f",robot->x,robot->y,robot->theta,\
           pmax->x,pmax->y,pmax->theta,sqrt(((robot->x-pmax->x)*(robot->x-pmax->x))+((robot->y-pmax->y)*(robot->y-pmax->y))));
   glColor3f(1.0,1.0,1.0);
   glRasterPos2i(5,22);
   for (int i=0; i<strlen(&line[0]); i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,(int)line[i]);
  }

  // Make sure all OpenGL commands are executed
  glFlush();
  // Swap buffers to enable smooth animation
  glutSwapBuffers();

  glDeleteTextures( 1, &texture );

  // Tell glut window to update ls itself
  glutSetWindow(windowID);
  glutPostRedisplay();

  if (first_frame)
  {
   fprintf(stderr,"All set! press enter to start\n");
   gets(&line[0]);
   first_frame=0;
  }
}

/*********************************************************************
 OpenGL and display stuff follows, you do not need to read code
 below this line.
*********************************************************************/
// Initialize glut and create a window with the specified caption
void initGlut(char* winName)
{
    // Set video mode: double-buffered, color, depth-buffered
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Create window
    glutInitWindowPosition (0, 0);
    glutInitWindowSize(Win[0],Win[1]);
    windowID = glutCreateWindow(winName);

    // Setup callback functions to handle window-related events.
    // In particular, OpenGL has to be informed of which functions
    // to call when the image needs to be refreshed, and when the
    // image window is being resized.
    glutReshapeFunc(WindowReshape);   // Call WindowReshape whenever window resized
    glutDisplayFunc(ParticleFilterLoop);   // Main display function is also the main loop
    glutKeyboardFunc(kbHandler);
}

void kbHandler(unsigned char key, int x, int y)
{
 if (key=='r') {RESETflag=1;}
 if (key=='q') {deleteList(list); free(map); free(map_b); exit(0);}
}

void WindowReshape(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();			// Initialize with identity matrix
    gluOrtho2D(0, 800, 800, 0);
    glViewport(0,0,w,h);
    Win[0] = w;
    Win[1] = h;
}
