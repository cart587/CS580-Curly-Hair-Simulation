/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

enum CollisionPlane
{
	positiveX,
	positiveY,
	positiveZ,
	negativeX,
	negativeY,
	negativeZ
};

/*Checks if a mass located at currI, currJ, currK can have a neighbor at the location given by
offsets in the i,j, and k direction*/
bool validNeighbor(int currI, int currJ, int currK, int iOffset, int jOffset, int kOffset)
{
	int neighborI = currI + iOffset;
	int neighborJ = currJ + jOffset;
	int neighborK = currK + kOffset;

	return !((neighborI > 7) || (neighborI < 0) || (neighborJ > 7) || (neighborJ < 0) || (neighborK > 7) || (neighborK < 0));
}

/*Returns the magnitude of a given vector*/
double computeVectorMagnitude(struct point vector)
{
	return sqrt(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2));
}

/*ADDED: Returns vector1-vector2*/
point computeVectorAddition(struct point vector1, struct point vector2)
{
	point res;
	res.x = vector1.x + vector2.x;
	res.y = vector1.y + vector2.y;
	res.z = vector1.z + vector2.z;
	return res;
}

/*ADDED: Returns vector1-vector2*/
point computeVectorSubtraction(struct point vector1, struct point vector2)
{
	point res;
	res.x = vector1.x - vector2.x;
	res.y = vector1.y - vector2.y;
	res.z = vector1.z - vector2.z;
	return res;
}

/*ADDED: Returns vector1*scale*/
point computeVectorScale(struct point vector, double scale)
{
	point res;
	res.x = vector.x * scale;
	res.y = vector.y * scale;
	res.z = vector.z * scale;
	return res;
}

/*ADDED: Apply smoothing function to spring points*/
void springSmoothingFunction(struct world* jello, double alpha)
{
	// compute avg rest length
	double totalLength = 0.0;

	for (int i = 1; i < 8; i++)
	{
		point currVec = computeVectorSubtraction(jello->p_init[i][0][0], jello->p_init[i - 1][0][0]);

		totalLength += computeVectorMagnitude(currVec);
	}

	double avgRestLength = totalLength / 7.0;
	double beta = min(1, 1 - exp(-avgRestLength / alpha));

	point d[8][8][8];

	// compute d
	for (int i = 0; i < 7; i++)
	{
		point dminus1 = (i - 1 < 0) ? computeVectorSubtraction(jello->p[1][0][0], jello->p[0][0][0]) : d[i - 1][0][0];
		point dminus2 = (i - 2 < 0) ? computeVectorSubtraction(jello->p[1][0][0], jello->p[0][0][0]) : d[i - 2][0][0];
		point nextPointMinusThis = computeVectorSubtraction(jello->p[i + 1][0][0], jello->p[i][0][0]);

		point firstTerm = computeVectorScale(dminus1, (2 * (1 - beta)));
		point secondTerm = computeVectorScale(dminus2, (1 - beta)*(1 - beta));
		point thirdTerm = computeVectorScale(nextPointMinusThis, beta*beta);

		d[i][0][0] = computeVectorSubtraction(firstTerm, secondTerm);
		d[i][0][0] = computeVectorAddition(d[i][0][0], thirdTerm);
	}

	// recompute p after the root
	for (int i = 1; i < 8; i++)
	{
		jello->p[i][0][0] = computeVectorAddition(jello->p[i - 1][0][0], d[i - 1][0][0]);
	}
}

/*Computes the force that the origin mass has on the neighbor mass and adds the result to the total force for the neighbor*/
void applySpringForce(struct point origin, struct point neighbor, struct point * force, struct world * jello, double restLength)
{
	struct point vector;
	double kCoefficient = jello->kElastic;

	//Length vector from origin to neighbor
	vector.x = neighbor.x - origin.x;
	vector.y = neighbor.y - origin.y;
	vector.z = neighbor.z - origin.z;

	//calculation of the scalars in Hook's law equation
	double currLength = computeVectorMagnitude(vector);
	double scalarResult = -1.0 * kCoefficient * (currLength - restLength) * (1.0 / currLength);

	force->x += (scalarResult * vector.x);
	force->y += (scalarResult * vector.y);
	force->z += (scalarResult * vector.z);

}

/*Calculates the damping force that should be applied to the force that the origin mass has on the neighbor mass.
Result is accumulated into the neighbor's corresponding location in the force array*/
void applyDampingForce(struct point origin, struct point neighbor, struct point * force, \
						struct world * jello, struct point oVelocity, struct point nVelocity, double dCoefficient)
{
	struct point lengthVector, velocityVector;
	double scalarResult, currLength;

	//subtract velocity vectors
	velocityVector.x = nVelocity.x - oVelocity.x;
	velocityVector.y = nVelocity.y - oVelocity.y;
	velocityVector.z = nVelocity.z - oVelocity.z;

	//calculate length vector
	lengthVector.x = neighbor.x - origin.x;
	lengthVector.y = neighbor.y - origin.y;
	lengthVector.z = neighbor.z - origin.z;
	currLength = computeVectorMagnitude(lengthVector);

	//dot product
	scalarResult = (velocityVector.x * lengthVector.x) + (velocityVector.y * lengthVector.y) + (velocityVector.z * lengthVector.z);
	scalarResult *= -1.0 * dCoefficient * (1.0 / pow(currLength, 2));
	
	force->x += (scalarResult * lengthVector.x);
	force->y += (scalarResult * lengthVector.y);
	force->z += (scalarResult * lengthVector.z);
}

/*Takes the mass represented by the index [bi][bj][bk] and applies spring forces to the 6 structural neighbors of the mass.
Results are stored in the forces array*/
void applyStructuralForceToNeighbors(int bi, int bj, int bk, struct world *jello, struct point forces[8][8][8])
{
	struct point neighborPoint, nVelocity;
	struct point originPoint = jello->p[bi][bj][bk];
	struct point originVelocity = jello->v[bi][bj][bk];
	double kCoefficient = jello->kElastic;
	double restLength = 1.0 / 7.0;

	//Neighbor in positive x direction
	if (validNeighbor(bi, bj, bk, 1, 0, 0)) 
	{
		neighborPoint = jello->p[bi + 1][bj][bk];
		nVelocity = jello->v[bi + 1][bj][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}
}

/*Takes the mass represented by the index [bi][bj][bk] and applies spring forces to the 20 shear neighbors of the mass.
Results are stored in the forces array*/
void applyShearForceToNeighbors(int bi, int bj, int bk, struct world *jello, struct point forces[8][8][8])
{
	struct point neighborPoint, nVelocity;
	struct point originPoint = jello->p[bi][bj][bk];
	struct point originVelocity = jello->v[bi][bj][bk];
	double kCoefficient = jello->kElastic;
	double restLength = sqrt(2.0) / 7.0;
	double longRestLength = sqrt(3.0) / 7.0;

	//There are 20 Shear neighbors

	if (validNeighbor(bi, bj, bk, 1, -1, -1))
	{
		neighborPoint = jello->p[bi + 1][bj - 1][bk - 1];
		nVelocity = jello->v[bi + 1][bj - 1][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj - 1][bk - 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj - 1][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, -1, -1))
	{
		neighborPoint = jello->p[bi - 1][bj - 1][bk - 1];
		nVelocity = jello->v[bi - 1][bj - 1][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj - 1][bk - 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj - 1][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, 1, -1))
	{
		neighborPoint = jello->p[bi - 1][bj + 1][bk - 1];
		nVelocity = jello->v[bi - 1][bj + 1][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj + 1][bk - 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj + 1][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, 1, -1))
	{
		neighborPoint = jello->p[bi + 1][bj + 1][bk - 1];
		nVelocity = jello->v[bi + 1][bj + 1][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj + 1][bk - 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj + 1][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, -1, 1))
	{
		neighborPoint = jello->p[bi + 1][bj - 1][bk + 1];
		nVelocity = jello->v[bi + 1][bj - 1][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj - 1][bk + 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj - 1][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, -1, 1))
	{
		neighborPoint = jello->p[bi - 1][bj - 1][bk + 1];
		nVelocity = jello->v[bi - 1][bj - 1][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj - 1][bk + 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj - 1][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, 1, 1))
	{
		neighborPoint = jello->p[bi - 1][bj + 1][bk + 1];
		nVelocity = jello->v[bi - 1][bj + 1][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj + 1][bk + 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj + 1][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, 1, 1))
	{
		neighborPoint = jello->p[bi + 1][bj + 1][bk + 1];
		nVelocity = jello->v[bi + 1][bj + 1][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj + 1][bk + 1], jello, longRestLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj + 1][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, 1, 0))
	{
		neighborPoint = jello->p[bi + 1][bj + 1][bk];
		nVelocity = jello->v[bi + 1][bj + 1][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj + 1][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj + 1][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, 1, 0))
	{
		neighborPoint = jello->p[bi - 1][bj + 1][bk];
		nVelocity = jello->v[bi - 1][bj + 1][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj + 1][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj + 1][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, -1, 0))
	{
		neighborPoint = jello->p[bi - 1][bj - 1][bk];
		nVelocity = jello->v[bi - 1][bj - 1][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj - 1][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj - 1][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, -1, 0))
	{
		neighborPoint = jello->p[bi + 1][bj - 1][bk];
		nVelocity = jello->v[bi + 1][bj - 1][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj - 1][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj - 1][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, 1, 1))
	{
		neighborPoint = jello->p[bi][bj + 1][bk + 1];
		nVelocity = jello->v[bi][bj + 1][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj + 1][bk + 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj + 1][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, -1, 1))
	{
		neighborPoint = jello->p[bi][bj - 1][bk + 1];
		nVelocity = jello->v[bi][bj - 1][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj - 1][bk + 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj - 1][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, -1, -1))
	{
		neighborPoint = jello->p[bi][bj - 1][bk - 1];
		nVelocity = jello->v[bi][bj - 1][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj - 1][bk - 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj - 1][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, 1, -1))
	{
		neighborPoint = jello->p[bi][bj + 1][bk - 1];
		nVelocity = jello->v[bi][bj + 1][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj + 1][bk - 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj + 1][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, 0, 1))
	{
		neighborPoint = jello->p[bi + 1][bj][bk + 1];
		nVelocity = jello->v[bi + 1][bj][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj][bk + 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, 0, 1))
	{
		neighborPoint = jello->p[bi - 1][bj][bk + 1];
		nVelocity = jello->v[bi - 1][bj][bk + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj][bk + 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj][bk + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -1, 0, -1))
	{
		neighborPoint = jello->p[bi - 1][bj][bk - 1];
		nVelocity = jello->v[bi - 1][bj][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 1][bj][bk - 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 1][bj][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 1, 0, -1))
	{
		neighborPoint = jello->p[bi + 1][bj][bk - 1];
		nVelocity = jello->v[bi + 1][bj][bk - 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1][bj][bk - 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1][bj][bk - 1], jello, originVelocity, nVelocity, jello->dElastic);
	}
}

/*Takes the mass represented by the index [bi][bj][bk] and applies spring forces to the 6 bend neighbors of the mass.
Results are stored in the forces array*/
void applyBendForceToNeighbors(int bi, int bj, int bk, struct world *jello, struct point forces[8][8][8])
{
	struct point neighborPoint, nVelocity;
	struct point originPoint = jello->p[bi][bj][bk];
	struct point originVelocity = jello->v[bi][bj][bk];
	double kCoefficient = jello->kElastic;
	double restLength = 2.0 / 7.0;

	//There are 6 Bend neighbors

	if (validNeighbor(bi, bj, bk, 0, 0, -2))
	{
		neighborPoint = jello->p[bi][bj][bk - 2];
		nVelocity = jello->v[bi][bj][bk - 2];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj][bk - 2], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj][bk - 2], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, -2, 0))
	{
		neighborPoint = jello->p[bi][bj - 2][bk];
		nVelocity = jello->v[bi][bj - 2][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj - 2][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj - 2][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, -2, 0, 0))
	{
		neighborPoint = jello->p[bi - 2][bj][bk];
		nVelocity = jello->v[bi - 2][bj][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi - 2][bj][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi - 2][bj][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, 0, 2))
	{
		neighborPoint = jello->p[bi][bj][bk + 2];
		nVelocity = jello->v[bi][bj][bk + 2];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj][bk + 2], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj][bk + 2], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 0, 2, 0))
	{
		neighborPoint = jello->p[bi][bj + 2][bk];
		nVelocity = jello->v[bi][bj + 2][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi][bj + 2][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi][bj + 2][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}

	if (validNeighbor(bi, bj, bk, 2, 0, 0))
	{
		neighborPoint = jello->p[bi + 2][bj][bk];
		nVelocity = jello->v[bi + 2][bj][bk];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 2][bj][bk], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 2][bj][bk], jello, originVelocity, nVelocity, jello->dElastic);
	}
}

/*Applies a collision spring force to the mass represented by index [bi][bj][bk] and accumulates the result into the
the mass's corresponding force index. CollisionPlane parameter specifies what plane the mass collided with. This 
functionn should only be called if an actual collision has been detected. This function does not check if a mass
point has actually collided. The damping force for the collision is computed at the end*/
void applyCollisionForce(int bi, int bj, int bk, struct world * jello, struct point forces[8][8][8], CollisionPlane cp)
{
	double scalarResult = -1.0 * jello->kCollision;
	struct point lengthVector, position, projection, projVelocity;

	position = jello->p[bi][bj][bk];
	lengthVector.x = 0.0;
	lengthVector.y = 0.0;
	lengthVector.z = 0.0;

	projection.x = position.x;
	projection.y = position.y;
	projection.z = position.z;
	projVelocity.x = 0.0;
	projVelocity.y = 0.0;
	projVelocity.z = 0.0;

	//Sets the projected point and length vector of the mass that crossed box boundaries
	//based on the plane the mass collided with
	switch (cp)
	{
	case positiveX:
		lengthVector.x = position.x - 2;
		projection.x = 2;
		break;
	case positiveY:
		lengthVector.y = position.y - 2;
		projection.y = 2;
		break;
	case positiveZ:
		lengthVector.z = position.z - 2;
		projection.z = 2;
		break;
	case negativeX:
		lengthVector.x = position.x + 2;
		projection.x = -2;
		break;
	case negativeY:
		lengthVector.y = position.y + 2;
		projection.y = -2;
		break;
	case negativeZ:
		lengthVector.z = position.z + 2;
		projection.z = -2;
		break;
	default:
		break;
	}

	//accumulate spring force into array
	forces[bi][bj][bk].x += (scalarResult * lengthVector.x);
	forces[bi][bj][bk].y += (scalarResult * lengthVector.y);
	forces[bi][bj][bk].z += (scalarResult * lengthVector.z);

	//Applies the damping force for collisions
	applyDampingForce(projection, position, &forces[bi][bj][bk], jello, projVelocity, jello->v[bi][bj][bk], jello->dCollision);
}

/*This function determines if the mass represented by index [bi][bj][bk] has passed the [-2,2] boundary in the
X,Y,Z directions. If so, a collision force will be applied to the mass*/
void checkForCollision(int bi, int bj, int bk, struct world * jello, struct point forces[8][8][8])
{
	struct point position = jello->p[bi][bj][bk];

	//Collision in x plane boundaries
	if (position.x > 2.0)
	{
		applyCollisionForce(bi, bj, bk, jello, forces, positiveX);
	} 
	else if (position.x < -2.0)
	{
		applyCollisionForce(bi, bj, bk, jello, forces, negativeX);
	}

	//Collision in y plane boundaries
	if (position.y > 2.0)
	{
		applyCollisionForce(bi, bj, bk, jello, forces, positiveY);
	}
	else if (position.y < -2.0)
	{
		applyCollisionForce(bi, bj, bk, jello, forces, negativeY);
	}

	//Collision in z plane boundaries
	if (position.z > 2.0)
	{
		applyCollisionForce(bi, bj, bk, jello, forces, positiveZ);
	}
	else if (position.z < -2.0)
	{
		applyCollisionForce(bi, bj, bk, jello, forces, negativeZ);
	}
}

/*Uses trilinear interpolation to apply the force field, if any, acting on the mass represented by index
[bi][bj][bk].*/
void checkForForceField(int mi, int mj, int mk, struct world *jello, struct point forces[8][8][8])
{
	struct point position = jello->p[mi][mj][mk];
	const int resolution = jello->resolution;
	const double boundingBoxMin = -2.0;
	const double boundingBoxMax = 2.0;
	const double boundingBoxLength = 4.0;
	double dist = boundingBoxLength / resolution;

	//No force field if resolution is less than 1
	if (resolution <= 1)
		return;

	//Mass point is out of bounds of force field
	if (position.x > 2.0 || position.x < -2.0 || position.y > 2.0 || position.y < -2.0 || position.z > 2 || position.z < -2.0)
		return;
	
	struct point totalForceField = { 0.0, 0.0, 0.0 };
	struct point f000, f001, f010, f011, f100, f101, f110, f111;
	double xLow, xHigh, yLow, yHigh, zLow, zHigh;
	double distanceFromX, distanceFromY, distanceFromZ;
	double x, y, z;
	int tempI, tempJ, tempK;

	//Calculcate the boundaries of the voxel
	distanceFromX = position.x - boundingBoxMin;
	distanceFromY = position.y - boundingBoxMin;
	distanceFromZ = position.z - boundingBoxMin;

	xLow = boundingBoxMin + floor(distanceFromX / dist) * dist;
	yLow = boundingBoxMin + floor(distanceFromY / dist) * dist;
	zLow = boundingBoxMin + floor(distanceFromZ / dist) * dist;
	xHigh = boundingBoxMin + ceil(distanceFromX / dist) * dist;
	yHigh = boundingBoxMin + ceil(distanceFromY / dist) * dist;
	zHigh = boundingBoxMin + ceil(distanceFromZ / dist) * dist;

	//If mass is perfectly in voxel plane, xLow and xHigh will come out to same value,
	//so we update the values to select proper boundaries and prevent divide by 0 errors
	double maxVoxel = boundingBoxMin + ceil(boundingBoxLength / dist) * dist;
	if (xLow == xHigh) {
		if (xHigh + dist <= maxVoxel)
			xHigh += dist;
		else
			xLow -= dist;
	}

	if (yLow == yHigh) {
		if (yHigh + dist <= maxVoxel)
			yHigh += dist;
		else
			yLow -= dist;
	}

	if (zLow == zHigh) {
		if (zHigh + dist <= maxVoxel)
			zHigh += dist;
		else
			zLow -= dist;
	}

	//Calculcate how far along the x,y,z voxel the mass is
	x = (position.x - xLow) / (xHigh - xLow);
	y = (position.y - yLow) / (yHigh - yLow);
	z = (position.z - zLow) / (zHigh - zLow);

	//Get the forces at each corner of the voxel and interpolate the force
	//that it applies to the mass point
	// temp = (distance from -2 to xLow/xHigh) * (distance between voxel points)
	//f000 coordinates
	tempI = (int)((xLow - boundingBoxMin)/boundingBoxLength * (resolution - 1));
	tempJ = (int)((yLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f000 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (1 - x)*(1 - y)*(1 - z)*f000.x;
	totalForceField.y += (1 - x)*(1 - y)*(1 - z)*f000.y;
	totalForceField.z += (1 - x)*(1 - y)*(1 - z)*f000.z;

	//f001 coordinates
	tempI = (int)((xLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f001 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (1 - x)*(1 - y)*(z)*f001.x;
	totalForceField.y += (1 - x)*(1 - y)*(z)*f001.y;
	totalForceField.z += (1 - x)*(1 - y)*(z)*f001.z;

	//f010 coordinates
	tempI = (int)((xLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f010 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (1 - x)*(y)*(1 - z)*f010.x;
	totalForceField.y += (1 - x)*(y)*(1 - z)*f010.y;
	totalForceField.z += (1 - x)*(y)*(1 - z)*f010.z;

	//f011 coordinates
	tempI = (int)((xLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f011 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (1 - x)*(y)*(z)*f011.x;
	totalForceField.y += (1 - x)*(y)*(z)*f011.y;
	totalForceField.z += (1 - x)*(y)*(z)*f011.z;

	//f100
	tempI = (int)((xHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f100 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (x)*(1 - y)*(1 - z)*f100.x;
	totalForceField.y += (x)*(1 - y)*(1 - z)*f100.y;
	totalForceField.z += (x)*(1 - y)*(1 - z)*f100.z;

	//f101
	tempI = (int)((xHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f101 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (x)*(1 - y)*(z)*f101.x;
	totalForceField.y += (x)*(1 - y)*(z)*f101.y;
	totalForceField.z += (x)*(1 - y)*(z)*f101.z;

	//f110
	tempI = (int)((xHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zLow - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f110 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (x)*(y)*(1 - z)*f110.x;
	totalForceField.y += (x)*(y)*(1 - z)*f110.y;
	totalForceField.z += (x)*(y)*(1 - z)*f110.z;

	//f111
	tempI = (int)((xHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempJ = (int)((yHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	tempK = (int)((zHigh - boundingBoxMin) / boundingBoxLength * (resolution - 1));
	f111 = jello->forceField[tempI * resolution * resolution + tempJ * resolution + tempK];
	totalForceField.x += (x)*(y)*(z)*f111.x;
	totalForceField.y += (x)*(y)*(z)*f111.y;
	totalForceField.z += (x)*(y)*(z)*f111.z;

	//Add the total force from the voxel to the mass's accumulated forces
	forces[mi][mj][mk].x += totalForceField.x;
	forces[mi][mj][mk].y += totalForceField.y;
	forces[mi][mj][mk].z += totalForceField.z;
}

/*Computes the force field applied to each mass*/
void computeExternalForces(struct world *jello, struct point forces[8][8][8])
{
	int i, j = 0, k = 0;
	for (i = 0; i < 8; i++) {
		checkForForceField(i, j, k, jello, forces);
	}
}

/*Checks and applies spring force if any of the masses have collided with the boundary planes*/
void computeCollisionSpringForces(struct world * jello, struct point forces[8][8][8])
{
	int i, j = 0, k = 0;
	for (i = 0; i < 8; i++) {
		checkForCollision(i, j, k, jello, forces);
	}
}

/*	Calculates the force each mass-point has on its bend neighbors, and
	sums the total force at each point. Returns forces on each mass
	in result array 'forces'*/
void computeBendSpringForces(struct world * jello, struct point forces[8][8][8])
{
	int i, j, k;
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			for (k = 0; k < 8; k++) {
				applyBendForceToNeighbors(i, j, k, jello, forces);
			}
		}
	}
}

/*	Calculates the force each mass-point has on its shear neighbors, and
	sums the total force at each point. Returns forces on each mass
	in result array 'forces'*/
void computeShearSpringForces(struct world * jello, struct point forces[8][8][8])
{
	int i, j, k;
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			for (k = 0; k < 8; k++) {
				applyShearForceToNeighbors(i, j, k, jello, forces);
			}
		}
	}
}

/*	Calculates the force each mass-point has on its structural neighbors, and
	sums the total force at each point. Returns forces on each mass
	in result array 'forces'*/
void computeStructuralSpringForces(struct world * jello, struct point forces[8][8][8])
{
	int i, j = 0, k = 0;
	for (i = 0; i < 8; i++) {
		applyStructuralForceToNeighbors(i, j, k, jello, forces);
	}
}

/*Takes the total force applied to each mass and computes the acceleration by dividing each component
by the mass. Function should be called after all shear, bend, structural, collision, and external forces
have been calculated and summed*/
void computeEachMassAcceleration(struct world * jello, struct point forces[8][8][8], struct point a[8][8][8])
{
	int i, j = 0, k = 0;

	for (i = 0; i < 8; i++) {
		a[i][j][k].x = forces[i][j][k].x / jello->mass;
		a[i][j][k].y = forces[i][j][k].y / jello->mass;
		a[i][j][k].z = forces[i][j][k].z / jello->mass;
	}
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
	struct point forces[8][8][8];
	int i, j, k;
	double alpha = 6.0;

	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			for (k = 0; k < 8; k++) {
				forces[i][j][k].x = 0.0;
				forces[i][j][k].y = 0.0;
				forces[i][j][k].z = 0.0;
			}
		}
	}

	computeStructuralSpringForces(jello, forces);
	springSmoothingFunction(jello, alpha);
	// computeShearSpringForces(jello, forces);
	// computeBendSpringForces(jello, forces);
	computeCollisionSpringForces(jello, forces);
	computeExternalForces(jello, forces);

	computeEachMassAcceleration(jello, forces, a);

}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j = 0,k = 0;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  //keep first mass frozen
  jello->p[0][j][k].x = 0;
  jello->p[0][j][k].y = 0;
  jello->p[0][j][k].z = 0;
  jello->v[0][j][k].x = 0;
  jello->v[0][j][k].y = 0;
  jello->v[0][j][k].z = 0;

  for (i=1; i<=7; i++)
  {
	jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
	jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
	jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
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
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
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
