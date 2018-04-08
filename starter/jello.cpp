/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Your name:
  <write your name here>

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"
#include <string>


// camera parameters
double Theta = piJello / 6;
double Phi = piJello / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=0, viewingMode=0, saveScreenToFile=0;

struct point applyForceDeltaMouse;

struct world jello;

int windowWidth, windowHeight;

void myinit()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90.0,1.0,0.01,1000.0);

  // set background color to grey
  glClearColor(0.5, 0.5, 0.5, 0.0);

  //glCullFace(GL_BACK); - turn off face culling for plane
  //glEnable(GL_CULL_FACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  
  double maxLength = sqrt(3 * pow((4.0/jello.resolution), 2));
  
  
  jello.fallOffEpsilon = log(10E-2)/(-1 * maxLength); //makes fall off such that it goes to 0.01 for largest possible distance between points (in fact slightly larger than largest possible distance)
  
  
  //setup values for inclined plane
  jello.a = 1;
  jello.b = 0.0;
  jello.c = 3.0;
  jello.d = 6.0;
  
  return; 
}

void reshape(int w, int h) 
{
  // Prevent a divide by zero, when h is zero.
  // You can't make a window of zero height.
  if(h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Set the perspective
  double aspectRatio = 1.0 * w / h;
  gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);
  
  
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 

  windowWidth = w;
  windowHeight = h;

  glutPostRedisplay();
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  // camera parameters are Phi, Theta, R
  gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
	        0.0,0.0,0.0, 0.0,0.0,1.0);


  /* Lighting */
  /* You are encouraged to change lighting parameters or make improvements/modifications
     to the lighting model . 
     This way, you will personalize your assignment and your assignment will stick out. 
  */

  // global ambient light
  GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };
  
  // light 's ambient, diffuse, specular
  GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd0[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd1[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs1[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd2[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs2[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd3[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs3[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd4[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs4[] = { 0.0, 0.0, 0.0, 1.0 };

  GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd5[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs5[] = { 1.0, 1.0, 1.0, 1.0};

  GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd6[] = { 0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs6[] = { 1.0, 1.0, 1.0, 1.0};

  GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd7[] = {  0.578125, 0.0, 0.078125, 1.0 };
  GLfloat lKs7[] = {1.0, 1.0, 1.0, 1.0};

  // light positions and directions
  GLfloat lP0[] = { -1.999, -1.999, -1.999, 1.0 };
  GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
  GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
  GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
  GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
  GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
  GLfloat lP6[] = { 1.999, 1.999, 1.999, 1.0 };
  GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
  
  // jelly material color

  GLfloat mKa[] = { 0.0, 0.0, 0.0, 0.5 };
  GLfloat mKd[] = { 0.6, 0.6, 0.6, 0.5 };
  GLfloat mKs[] = { 0.3, 0.3, 0.3, 0.3 };
  GLfloat mKe[] = { 0.0, 0.0, 0.0, 1.0 };

  /* set up lighting */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  // set up cube color
  glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
  glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
  glMaterialf(GL_FRONT, GL_SHININESS, 120);
    
  // macro to set up light i
  #define LIGHTSETUP(i)\
  glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
  glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
  glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
  glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
  glEnable(GL_LIGHT##i)
  
  LIGHTSETUP (0);
  LIGHTSETUP (1);
  LIGHTSETUP (2);
  LIGHTSETUP (3);
  LIGHTSETUP (4);
  LIGHTSETUP (5);
  LIGHTSETUP (6);
  LIGHTSETUP (7);

  // enable lighting
  glEnable(GL_LIGHTING);    
  glEnable(GL_DEPTH_TEST);

  // show the cube
  showCube(&jello);

  glDisable(GL_LIGHTING);

  // show the bounding box
  showBoundingBox();
  
  showPlane();
 
  glutSwapBuffers();
}

void showPlane()
{
  
  
  glColor4f(0,0.6,0.6,0);
  glBegin(GL_TRIANGLE_STRIP);
  double a = jello.a;
  double b = jello.b;
  double c = jello.c;
  double d = jello.d;
  

  
  //z
  //calculateIntersect(a, b, d, c, -2, 2);
  //calculateIntersect(a, b, d, c, 2, 2);
  
  //y
  calculateIntersect(a, c, d, b, -2, 1);
  calculateIntersect(a, c, d, b, 2, 1);
  
  //x
  
  //calculateIntersect(b, c, d, a, 2, 0);
  //calculateIntersect(b, c, d, a, -2, 0);
  
  
   glEnd();
//  glColor4f(0.6,0,0,0);
//  glBegin(GL_LINES);
//  glVertex3f(0, 0, -2.0);
//  glVertex3f(1, 0, 1.0);
//  glEnd();
  
  return;
  
  
  
  
}


void calculateIntersect(double dim1,double dim2, double d, double dim3, double planeValue, int whichPlane) {
  double xmaxIntersect = (-1.0 * ((dim1 * 2) + dim3 * planeValue  + d))/dim2;
  double xminIntersect = (-1.0 * ((dim1 * -2) + dim3 * planeValue + d))/dim2;
  double x1;
  double y1;
  
  double x2;
  double y2;
  bool set1 = false;
  bool set2 = false;
  if(xminIntersect<=2 && xminIntersect>= -2)
  {
    x1 = -2;
    y1 = xminIntersect;
    set1= true;
    if (xmaxIntersect <= 2 && xmaxIntersect >= -2) {
      x2 = 2;
      y2 = xmaxIntersect;
      set2 = true;
    }
  } else if (xmaxIntersect <= 2 && xmaxIntersect >= -2) {
    x1 = 2;
    y1 = xmaxIntersect;
    set1 = true;
  }
  
  
  
  
  int yminIntersect = (-1.0 * ((dim2 * -2) + dim3 * planeValue + d))/dim1;
  int ymaxIntersect = (-1.0 * ((dim2 * 2) + dim3 * planeValue + d))/dim1;
 
  if(!set2) {
    if(yminIntersect<= 2 && yminIntersect>= -2)
    {
      
      if(set1) {
        y2 = -2;
        x2 = yminIntersect;
      } else {
        y1 = -2;
        x1 = yminIntersect;
        y2 = 2;
        x2 = ymaxIntersect;
      }
      
      
      
    } else {
      y2 = 2;
      x2 = ymaxIntersect;
    }
  }
  
  if(whichPlane == 0) {
    glVertex3f(planeValue,x1,y1);
  } else if (whichPlane == 1) {
    glVertex3f(x1,planeValue,y1);
  } else {
    glVertex3f(x1,y1, planeValue);
  }
  
  if(whichPlane == 0) {
    glVertex3f(planeValue,x2,y2);
  } else if (whichPlane == 1) {
    glVertex3f(x2,planeValue,y2);
  } else {
    glVertex3f(x2,y2, planeValue);
  }
}

void doIdle()
{
  char s[20]="picxxxx.ppm";
  //int i;
  
  // save screen to file
  s[3] = 48 + (sprite / 1000);
  s[4] = 48 + (sprite % 1000) / 100;
  s[5] = 48 + (sprite % 100 ) / 10;
  s[6] = 48 + sprite % 10;

  if (saveScreenToFile==1)
  {
    saveScreenshot(windowWidth, windowHeight, s);
    saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
  }

  if (sprite >= 300) // allow only 300 snapshots
  {
    exit(0);	
  }

  if (pause == 0)
  {
    
    if(std::string(jello.integrator).compare("RK4") == 0) {
      
      RK4(&jello);
      
    } else {
      Euler(&jello);
    }
    
    // insert code which appropriately performs one step of the cube simulation:
  }

  

  
  glutPostRedisplay();
}

int main (int argc, char ** argv)
{
  if (argc<2)
  {  
    printf ("Oops! You didn't say the jello world file!\n");
    printf ("Usage: %s [worldfile]\n", argv[0]);
    exit(0);
  }

  readWorld(argv[1],&jello);

  glutInit(&argc,argv);
  
  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Jello cube");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  return(0);
}

