/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[numMass]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
//void RK4(struct world * jello);

// ****ADDED*****
point computeVectorAddition(struct point vector1, struct point vector2);
point computeVectorSubtraction(struct point vector1, struct point vector2);
point computeVectorScale(struct point vector, double scale);
void springSmoothingFunction(struct world* jello, double beta);
void multiplyFrameByEdgeInit(double matrix[][3], double edge[][1], point &dst);
void ComputeInitialReferenceVectors(struct world* jello);
void transposeMatrix(double src[][3], double dst[][3]);
point computeVectorNormalized(point vector);
void ComputeReferenceVectors(struct world* jello);
void ComputeFrames(struct world* jello);
// **************

#endif

