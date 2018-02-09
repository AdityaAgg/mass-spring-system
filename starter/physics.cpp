/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <iostream>


const double GLOBAL_REST_LENGTH = 1.0/7.0;
const double GLOBAL_REST_LENGTH_SHEAR = sqrt(2.0/49);
const double GLOBAL_REST_LENGTH_BEND = 2.0/7.0;
double prev_length_diff = -1000;
bool isIncreasing = false;
/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  
  
  calculateStructuralForces(jello, a);
  calculateShearForces(jello, a);
  calculateBendForces(jello, a);
  
  
  
}
//calculate hooke's law between two points
//note force on point a is -1 * calculateHookesLaw
struct point calculateHookesLaw(struct point a, struct point b, double restLength, double elasticConstant) {
  
  //difference between two points
  struct point diff;
  
  pDIFFERENCE(a, b , diff);
  
  
  //normalize for direction
  double lengthDouble = sqrt((diff).x * (diff).x + (diff).y * (diff).y + (diff).z * (diff).z);
  pNORMALIZE(diff);
  
  
  
  struct point finalOutput;
  
  pMULTIPLY(diff, -1 * (lengthDouble - restLength) * elasticConstant, finalOutput);
  
  return finalOutput;
}


void calculateStructuralForces(struct world * jello, struct point a[8][8][8]) { //array is passed by reference to function since it is a pointer
  
  double massReciprocal = 1/(jello->mass);
  
  
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        
        if(k<7)
        {

          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j][k+1], GLOBAL_REST_LENGTH, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j][k+1], accelOnB, a[i][j][k+1]); //add acceleration on point B
          
        }
        if(j<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j+1][k], GLOBAL_REST_LENGTH,  jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j+1][k], accelOnB, a[i][j+1][k]); //add acceleration on point B
        }
        
        if(i<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j][k], GLOBAL_REST_LENGTH, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j][k], accelOnB, a[i+1][j][k]); //add acceleration on point B
        }
      }
    }
  }
}

void calculateBendForces(struct world * jello, struct point a[8][8][8]) {
  double massReciprocal = 1/(jello->mass);



  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {

        //top right to bottom left diagonal
        if(k<6)
        {

          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j][k+2], GLOBAL_REST_LENGTH_BEND, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j][k+2], accelOnB, a[i][j][k+2]); //add acceleration on point B

        }
        if(j<6)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j+2][k], GLOBAL_REST_LENGTH_BEND,  jello->kElastic);
          pMULTIPLY(f , massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j+2][k], accelOnB, a[i][j+2][k]); //add acceleration on point B
        }
        if(i<6)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+2][j][k], GLOBAL_REST_LENGTH_BEND,  jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+2][j][k], accelOnB, a[i+2][j][k]); //add acceleration on point B
        }
      }
    }
  }
}

void calculateShearForces(struct world * jello, struct point a[8][8][8]) { //array is passed by reference to function since it is a pointer
  
  double massReciprocal = 1/(jello->mass);
  
  
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        
        //top right to bottom left diagonal
        if(k<7 && j<7)
        {
          
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j+1][k+1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j+1][k+1], accelOnB, a[i][j+1][k+1]); //add acceleration on point B
          
        }
        if(i<7 && j<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j+1][k], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j+1][k], accelOnB, a[i+1][j+1][k]); //add acceleration on point B
        }
        
        if(i<7 && k<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j][k+1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic);
          pMULTIPLY(f , massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j][k+1], accelOnB, a[i+1][j][k+1]); //add acceleration on point B
        }
        
        
        
        //top left to bottom right diagonal
        if(k<7 && j>0)
        {

          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j-1][k+1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j-1][k+1], accelOnB, a[i][j-1][k+1]); //add acceleration on point B

        }
        if(i<7 && j>0)
        {

          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j-1][k], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j-1][k], accelOnB, a[i+1][j-1][k]); //add acceleration on point B

        }

        if(i<7 && k>0)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j][k-1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j][k-1], accelOnB, a[i+1][j][k-1]); //add acceleration on point B
        }
        
       
        
      }
    }
  }
}



  




/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        a[i][j][k].x = 0;
        a[i][j][k].y = 0;
        a[i][j][k].z =0;
      }
    }
  }

  
  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        double changeX = jello->dt * jello->v[i][j][k].x;
        double changeY =  jello->dt * jello->v[i][j][k].y;
        double changeZ = jello->dt * jello->v[i][j][k].z;
        
        
        jello->p[i][j][k].x += changeX;
        jello->p[i][j][k].y += changeY;
        jello->p[i][j][k].z += changeZ;
        
        
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        a[i][j][k].x = 0;
        a[i][j][k].y = 0;
        a[i][j][k].z =0;
      }
    }
  }

  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }
  
  return;  
}
