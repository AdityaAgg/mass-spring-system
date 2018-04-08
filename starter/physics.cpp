/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <iostream>


const double GLOBAL_REST_LENGTH = 1.0/7.0;
const double GLOBAL_REST_LENGTH_SHEAR = sqrt(2.0/49);
const double GLOBAL_REST_LENGTH_SHEAR_DIAGONAL = sqrt(3.0/49);
const double GLOBAL_REST_LENGTH_BEND = 2.0/7.0;
double prev_length_diff = -1000;
bool isIncreasing = false;
/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        a[i][j][k].x = 0;
        a[i][j][k].y = 0;
        a[i][j][k].z =0;
      }
    }
  }
  calculateStructuralForces(jello, a);
  calculateShearForces(jello, a);
  calculateBendForces(jello, a);
  calculateCollisionForces(jello, a);
  addForceFieldForces(jello, a);
  
  
  
}


void addForceFieldForces(struct world * jello, struct point a[8][8][8]) {
  
  
  for(int i = 0; i<8; ++i)
  {
    for(int j = 0; j<8; ++j)
    {
      for(int k =0; k<8; ++k)
      {
        calculateForceforPoint(a[i][j][k], jello, i, j, k);
        pSUM(a[i][j][k], applyForceDeltaMouse, a[i][j][k]);
        
        
      }
    }
  }
  
  applyForceDeltaMouse.x  = (applyForceDeltaMouse.x > 10E-2)?(applyForceDeltaMouse.x * 0.99):0.0;
  applyForceDeltaMouse.y = (applyForceDeltaMouse.y > 10E-2)?(applyForceDeltaMouse.y * 0.99):0.0;
  applyForceDeltaMouse.z = (applyForceDeltaMouse.z > 10E-2)?(applyForceDeltaMouse.z * 0.99):0.0;
  
  
  
}

void calculateForceforPoint(struct point & vec, struct world * jello, int l, int m, int n) {
  int resolution = jello->resolution;
  
  double halfStep = 2.0/resolution;
  struct point totals = {0, 0, 0};
  double fullStep = halfStep * 2;
  int voxelX = (jello->p[l][m][n].x + 2)/(fullStep) - 1;
  int voxelY = (jello->p[l][m][n].y + 2)/(fullStep) - 1;
  int voxelZ = (jello->p[l][m][n].z + 2)/(fullStep) - 1;
  
  
  
  

  
  //printf("function x:%f y:%f z:%f\n", jello->p[l][m][n].x, jello->p[l][m][n].y, jello->p[l][m][n].z);
  
  
  
  
  
  
  //printf("function called voxelX: %d voxelY: %d voxelZ: %d\n", voxelX, voxelY, voxelZ);
  
  
  
  
  for(int i = voxelX; i<voxelX + 2; ++i)
  {
    for(int j = voxelY; j<voxelY + 2; ++j)
    {
      for(int k =voxelZ; k<voxelZ + 2; ++k)
      {
        
        
        if( !(((i < 0 || j < 0) || (k<0)) || ((i>=resolution) || (j>=resolution || k>= resolution))))
        {

          
          struct point voxelCenter = {(i * (fullStep) + halfStep) - 2, (j * (fullStep) + halfStep) - 2, (k * (fullStep) + halfStep) - 2};
          struct point diff;
          pDIFFERENCE(voxelCenter, jello->p[l][m][n], diff);
          
          double distance = sqrt((diff).x * (diff).x + (diff).y * (diff).y + (diff).z * (diff).z);

          double weight = exp(-1 * jello->fallOffEpsilon * distance);
          struct point dest;
          
          pMULTIPLY(jello->forceField[i*resolution + j*resolution + k],weight, dest);
          pSUM(dest, totals, totals);
          
        }
        
      }
    }
  }
  
  


    pMULTIPLY(totals, (1/(jello->mass)), totals);
    //printf("before x:%f y:%f z:%f \n", vec.x, vec.y, vec.z);
    pSUM(vec, totals, vec);
    

}



//calculate hooke's law between two points
//note force on point a is -1 * calculateHookesLaw
struct point calculateHookesLaw(struct point a, struct point b, struct point vA, struct point vB, double restLength, double elasticConstant, double dampingConstant) {
  
  
  

  //difference between two points
  struct point diff;
  pDIFFERENCE(a, b , diff);
  
  
  //normalize for direction
  double lengthDouble = sqrt((diff).x * (diff).x + (diff).y * (diff).y + (diff).z * (diff).z);
  struct point regDiff; // to deep copy original vector/non-normalized
  pCPY(diff, regDiff);
  pNORMALIZE(diff);
  
  
  //find hooke's force
  struct point finalOutput;
  pMULTIPLY(diff, -1 * (lengthDouble - restLength) * elasticConstant, finalOutput);
  
  
  
  
  //damping calculations
  
  struct point damping;
  double result;
  struct point diffV;
  pDIFFERENCE(vA, vB, diffV);
  DOTPRODUCT(diffV, regDiff, result);
  
  result /= lengthDouble;
  pMULTIPLY(diff, -1 * result * dampingConstant, damping);
  
  
  
  
  //sum damping and hooke's force vectors together
  pSUM(damping, finalOutput, finalOutput);

  return finalOutput;
}




void calculateCollisionForces(struct world * jello, struct point a[8][8][8]) {
  double massReciprocal = (1/jello->mass);
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        struct point position = jello->p[i][j][k];
        if(position.x < -2)  {
          double difference = -2 - position.x;
          double magnitude = sqrt(difference * difference);
          a[i][j][k].x += massReciprocal * (jello->kCollision * magnitude);
          a[i][j][k].x -= massReciprocal * jello->dCollision * jello->v[i][j][k].x;
        }
        
        if(position.y < -2)  {
          double difference = -2 - position.y;
          double magnitude = sqrt(difference * difference);
          a[i][j][k].y += massReciprocal * jello->kCollision * magnitude;
          a[i][j][k].y -= massReciprocal * jello->dCollision * jello->v[i][j][k].y;
        }
        
        if(position.z < -2)  {
          double difference = -2 - position.z;
          double magnitude = sqrt(difference * difference);
          a[i][j][k].z += massReciprocal * jello->kCollision * magnitude;
          a[i][j][k].z -= massReciprocal * jello->dCollision * jello->v[i][j][k].z;
        }
        
        
        
        
        if(position.x > 2)  {
          double difference = 2 - position.x;
          double magnitude = sqrt(difference * difference);
          a[i][j][k].x -= massReciprocal * jello->kCollision * magnitude;
          a[i][j][k].x -= massReciprocal * jello->dCollision * jello-> v[i][j][k].x;
        }
        
        if(position.y > 2)  {
          double difference = 2 - position.y;
          double magnitude = sqrt(difference * difference);
          a[i][j][k].y -= massReciprocal * jello->kCollision * magnitude;
          a[i][j][k].y -= massReciprocal * jello->dCollision * jello-> v[i][j][k].y;
        }
        
        
        if(position.z > 2)  {
          double difference = 2 - position.z;
          double magnitude = sqrt(difference * difference);
          a[i][j][k].z -= massReciprocal * jello->kCollision * magnitude;
          a[i][j][k].z -= massReciprocal * jello->dCollision * jello-> v[i][j][k].z;
        }
        
        
        //inclined plane collision calculations
        struct point normal = {jello->a, jello->b, jello->c};
        struct point otherValue = {jello->p[i][j][k].x, jello->p[i][j][k].y, jello->p[i][j][k].z + 2.0};
        pNORMALIZE(otherValue);
        pNORMALIZE(normal);
        double dotProduct = normal.x * (otherValue.x) + normal.y * (otherValue.y) + normal.z * (otherValue.z);
        
        if(dotProduct < 0) {
          
          double t = -1 * (position.x + 3.0 * position.z + 6.0)/(10.0);
          struct point normal= {jello->a, jello->b, jello->c};
          struct point projectionOnPlane;
          pMULTIPLY(normal, t, projectionOnPlane);
          pSUM(projectionOnPlane, position, projectionOnPlane);
          
          struct point planeVelocity = {0.0, 0.0, 0.0};
          
          struct point f = calculateHookesLaw(projectionOnPlane, position, planeVelocity, jello->v[i][j][k], 0, jello->kCollision, jello->dCollision);
          struct point ans;
          
          pMULTIPLY(f,   100 * -1 * massReciprocal, ans);
          pSUM(a[i][j][k], ans, a[i][j][k]);
          
        }
        
        
      }
    }
  }
}

void calculateStructuralForces(struct world * jello, struct point a[8][8][8]) { //array is passed by reference to function since it is a pointer
  
  double massReciprocal = 1/(jello->mass);
  
  
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        
        if(k<7)
        {

          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j][k+1], jello->v[i][j][k], jello->v[i][j][k+1], GLOBAL_REST_LENGTH, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j][k+1], accelOnB, a[i][j][k+1]); //add acceleration on point B
          
        }
        if(j<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j+1][k], jello->v[i][j][k], jello->v[i][j+1][k], GLOBAL_REST_LENGTH,  jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j+1][k], accelOnB, a[i][j+1][k]); //add acceleration on point B
        }
        
        if(i<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j][k], jello->v[i][j][k], jello->v[i+1][j][k], GLOBAL_REST_LENGTH, jello->kElastic, jello->dElastic);
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
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j][k+2], jello->v[i][j][k], jello->v[i][j][k+2], GLOBAL_REST_LENGTH_BEND, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j][k+2], accelOnB, a[i][j][k+2]); //add acceleration on point B

        }
        if(j<6)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j+2][k], jello->v[i][j][k], jello->v[i][j+2][k], GLOBAL_REST_LENGTH_BEND,  jello->kElastic, jello->dElastic);
          pMULTIPLY(f , massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j+2][k], accelOnB, a[i][j+2][k]); //add acceleration on point B
        }
        if(i<6)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+2][j][k], jello->v[i][j][k], jello->v[i+2][j][k], GLOBAL_REST_LENGTH_BEND,  jello->kElastic, jello->dElastic);
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
  
  double massReciprocal = 1.0/(jello->mass);
  
  
  
  for(int i = 0; i < 8; ++i) {
    for(int j = 0; j < 8; ++j) {
      for(int k = 0; k < 8; ++k) {
        
        //top right to bottom left diagonal
        if(k<7 && j<7)
        {
          
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j+1][k+1], jello->v[i][j][k], jello->v[i][j+1][k+1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j+1][k+1], accelOnB, a[i][j+1][k+1]); //add acceleration on point B
          
        }
        if(i<7 && j<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j+1][k], jello->v[i][j][k], jello->v[i+1][j+1][k], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A
          
          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j+1][k], accelOnB, a[i+1][j+1][k]); //add acceleration on point B
        }
        
        if(i<7 && k<7)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j][k+1], jello->v[i][j][k], jello->v[i+1][j][k+1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic, jello->dElastic);
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
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i][j-1][k+1],jello->v[i][j][k], jello->v[i][j-1][k+1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i][j-1][k+1], accelOnB, a[i][j-1][k+1]); //add acceleration on point B

        }
        if(i<7 && j>0)
        {

          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j-1][k], jello->v[i][j][k], jello->v[i+1][j-1][k], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j-1][k], accelOnB, a[i+1][j-1][k]); //add acceleration on point B

        }

        if(i<7 && k>0)
        {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j][k-1], jello->v[i][j][k], jello->v[i+1][j][k-1], GLOBAL_REST_LENGTH_SHEAR, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j][k-1], accelOnB, a[i+1][j][k-1]); //add acceleration on point B
        }
        
        
       //3D Diagonal
        if(i<7 && j>0 && k>0) {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j-1][k-1], jello->v[i][j][k], jello->v[i+1][j-1][k-1], GLOBAL_REST_LENGTH_SHEAR_DIAGONAL, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j-1][k-1], accelOnB, a[i+1][j-1][k-1]); //add acceleration on point B
        }

        if(i<7 && j>0 && k<7) {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j-1][k+1], jello->v[i][j][k], jello->v[i+1][j-1][k+1], GLOBAL_REST_LENGTH_SHEAR_DIAGONAL, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j-1][k+1], accelOnB, a[i+1][j-1][k+1]); //add acceleration on point B
        }

        if(i<7 && j<7 && k>0) {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j+1][k-1], jello->v[i][j][k], jello->v[i+1][j+1][k-1], GLOBAL_REST_LENGTH_SHEAR_DIAGONAL, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j+1][k-1], accelOnB, a[i+1][j+1][k-1]); //add acceleration on point B
        }

        if(i<7 && k<7 && j<7) {
          struct point accelOnA;
          struct point f = calculateHookesLaw(jello->p[i][j][k], jello->p[i+1][j+1][k+1], jello->v[i][j][k], jello->v[i+1][j+1][k+1], GLOBAL_REST_LENGTH_SHEAR_DIAGONAL, jello->kElastic, jello->dElastic);
          pMULTIPLY(f ,  massReciprocal, accelOnA);
          pSUM(a[i][j][k], accelOnA, a[i][j][k]); //add acceleration on point A

          struct point accelOnB;
          pMULTIPLY(f,  -1 * massReciprocal, accelOnB);
          pSUM(a[i+1][j+1][k+1], accelOnB, a[i+1][j+1][k+1]); //add acceleration on point B
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
         //F2p = dt * buffer.v;
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
