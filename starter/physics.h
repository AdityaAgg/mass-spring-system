/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);
struct point calculateHookesLaw(struct point a, struct point b, double restLength, double elasticConstant);
void calculateStructuralForces(struct world * jello, struct point  a[8][8][8]);
void calculateShearForces(struct world * jello, struct point a[8][8][8]);
void calculateBendForces(struct world * jello, struct point a[8][8][8]);
#endif

