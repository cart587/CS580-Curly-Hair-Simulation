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

// ****ADDED*****
point computeVectorAddition(struct point vector1, struct point vector2);
point computeVectorSubtraction(struct point vector1, struct point vector2);
point computeVectorScale(struct point vector, double scale);
void springSmoothingFunction(struct world* jello, double beta);
// **************

#endif

