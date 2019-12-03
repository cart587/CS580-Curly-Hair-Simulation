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
bool validNeighbor(int currI, int iOffset)
{
	int neighborI = currI + iOffset;

	return !((neighborI >= numMass) || (neighborI < 0));
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
	alpha = 10.0;
	// compute avg rest length
	double totalLength = 0.0;

	for (int i = 1; i < numMass; i++)
	{
		point currVec = computeVectorSubtraction(jello->p_init[i], jello->p_init[i - 1]);

		totalLength += computeVectorMagnitude(currVec);
	}

	double avgRestLength = totalLength / 7.0;
	double beta = (alpha <= 1e-5) ? 1 : min(1, 1 - exp(-avgRestLength / alpha));
	point d[numMass];

	// compute d
	for (int i = 0; i < numMass - 1; i++)
	{
		point dminus1 = (i - 1 < 0) ? computeVectorSubtraction(jello->p[1], jello->p[0]) : d[i - 1];
		point dminus2 = (i - 2 < 0) ? computeVectorSubtraction(jello->p[1], jello->p[0]) : d[i - 2];
		point nextPointMinusThis = computeVectorSubtraction(jello->p[i + 1], jello->p[i]);

		point firstTerm = computeVectorScale(dminus1, (2 * (1 - beta)));
		point secondTerm = computeVectorScale(dminus2, (1 - beta)*(1 - beta));
		point thirdTerm = computeVectorScale(nextPointMinusThis, beta*beta);

		d[i] = computeVectorSubtraction(firstTerm, secondTerm);
		d[i] = computeVectorAddition(d[i], thirdTerm);
	}

	jello->p_smoothed[0] = jello->p_init[0];

	// recompute p after the root
	for (int i = 1; i < numMass; i++)
	{
		jello->p_smoothed[i] = computeVectorAddition(jello->p_smoothed[i - 1], d[i - 1]);
	}
}

/*Computes the force that the origin mass has on the neighbor mass and adds the result to the total force for the neighbor*/
void applySpringForce(struct point origin, struct point neighbor, struct point * force, struct world * jello, double restLength)
{
	struct point vector;
	double kCoefficient = 500000.0;//orig:5000.0; demo: 0.0;

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

void applyBendSpringForce(struct point origin, struct point neighbor, struct point * force, struct world * jello)
{
	struct point vector;
	double kCoefficient = 10;//orig: 10.0; demo: 0.0,; //jello->kElastic;

	//Length vector from origin to neighbor
	vector.x = neighbor.x - origin.x;
	vector.y = neighbor.y - origin.y;
	vector.z = neighbor.z - origin.z;

	struct point result;
	pMULTIPLY(vector, kCoefficient, result);

	force->x += result.x;
	force->y += result.y;
	force->z += result.z;
}

/*Calculates the damping force that should be applied to the force that the origin mass has on the neighbor mass.
Result is accumulated into the neighbor's corresponding location in the force array*/
void applyDampingForce(struct point origin, struct point neighbor, struct point * force, \
						struct world * jello, struct point oVelocity, struct point nVelocity, double dCoefficient)
{
	dCoefficient = 10.25;//orig: 10.25; demo:

	struct point lengthVector; 
	struct point velocityVector;
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

void apply520DampingForce(struct point edge, struct point * force, \
	struct world * jello, struct point oVelocity, struct point nVelocity, double dCoefficient)
{
	dCoefficient = 1.25;//orig:1.25; demo:0.0;

	struct point velocityVector;
	double velocityEHatDot;

	//subtract velocity vectors
	velocityVector.x = nVelocity.x - oVelocity.x;
	velocityVector.y = nVelocity.y - oVelocity.y;
	velocityVector.z = nVelocity.z - oVelocity.z;

	//calculate length vector
	struct point eHat = computeVectorNormalized(edge);

	//dot product
	velocityEHatDot = (velocityVector.x * eHat.x) + (velocityVector.y * eHat.y) + (velocityVector.z * eHat.z);

	struct point projectedEHat, dampingDirection, result;
	pMULTIPLY(eHat, velocityEHatDot, projectedEHat);
	pDIFFERENCE(velocityVector, projectedEHat, dampingDirection);
	pMULTIPLY(dampingDirection, dCoefficient, result);

	force->x += result.x;
	force->y += result.y;
	force->z += result.z;
}

/*Takes the mass represented by the index [bi][bj][bk] and applies spring forces to the 6 structural neighbors of the mass.
Results are stored in the forces array*/
void applyStructuralForceToNeighbors(int bi, struct world *jello, struct point forces[numMass])
{
	struct point neighborPoint, nVelocity;
	struct point originPoint = jello->p[bi];
	struct point originVelocity = jello->v[bi];
	double kCoefficient = jello->kElastic;
	double restLength = computeVectorMagnitude(computeVectorSubtraction(jello->p_init[1], jello->p_init[0]));

	//Neighbor in positive x direction
	if (validNeighbor(bi, 1)) 
	{
		neighborPoint = jello->p[bi + 1];
		nVelocity = jello->v[bi + 1];
		applySpringForce(originPoint, neighborPoint, &forces[bi + 1], jello, restLength);
		applyDampingForce(originPoint, neighborPoint, &forces[bi + 1], jello, originVelocity, nVelocity, jello->dElastic);
	}
}

/*Takes the mass represented by the index [bi] and applies spring forces to the 6 bend neighbors of the mass.
Results are stored in the forces array*/
void applyBendForceToNeighbors(int bi, struct world *jello, struct point forces[numMass])
{
	struct point nVelocity = jello->v[bi + 1];
	struct point originVelocity = jello->v[bi];
	struct point edge;

	pDIFFERENCE(jello->p[bi], jello->p[bi+1], edge);
	applyBendSpringForce(jello->t[bi], edge, &forces[bi], jello);
	apply520DampingForce(edge, &forces[bi], jello, originVelocity, nVelocity, jello->dElastic);
}

/*Applies a collision spring force to the mass represented by index [bi][bj][bk] and accumulates the result into the
the mass's corresponding force index. CollisionPlane parameter specifies what plane the mass collided with. This 
functionn should only be called if an actual collision has been detected. This function does not check if a mass
point has actually collided. The damping force for the collision is computed at the end*/
void applyCollisionForce(int bi, struct world * jello, struct point forces[numMass], CollisionPlane cp)
{
	double scalarResult = -1.0 * jello->kCollision;
	struct point lengthVector, position, projection, projVelocity;

	position = jello->p[bi];
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
	forces[bi].x += (scalarResult * lengthVector.x);
	forces[bi].y += (scalarResult * lengthVector.y);
	forces[bi].z += (scalarResult * lengthVector.z);

	//Applies the damping force for collisions
	applyDampingForce(projection, position, &forces[bi], jello, projVelocity, jello->v[bi], jello->dCollision);
}

/*This function determines if the mass represented by index [bi][bj][bk] has passed the [-2,2] boundary in the
X,Y,Z directions. If so, a collision force will be applied to the mass*/
void checkForCollision(int bi, struct world * jello, struct point forces[numMass])
{
	struct point position = jello->p[bi];

	//Collision in x plane boundaries
	if (position.x > 2.0)
	{
		applyCollisionForce(bi, jello, forces, positiveX);
	} 
	else if (position.x < -2.0)
	{
		applyCollisionForce(bi, jello, forces, negativeX);
	}

	//Collision in y plane boundaries
	if (position.y > 2.0)
	{
		applyCollisionForce(bi, jello, forces, positiveY);
	}
	else if (position.y < -2.0)
	{
		applyCollisionForce(bi, jello, forces, negativeY);
	}

	//Collision in z plane boundaries
	if (position.z > 2.0)
	{
		applyCollisionForce(bi, jello, forces, positiveZ);
	}
	else if (position.z < -2.0)
	{
		applyCollisionForce(bi, jello, forces, negativeZ);
	}
}

/*Uses trilinear interpolation to apply the force field, if any, acting on the mass represented by index
[bi][bj][bk].*/
void checkForForceField(int mi, struct world *jello, struct point forces[numMass])
{
	struct point position = jello->p[mi];
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
	forces[mi].x += totalForceField.x;
	forces[mi].y += totalForceField.y;
	forces[mi].z += totalForceField.z;
}

/*Computes the force field applied to each mass*/
void computeExternalForces(struct world *jello, struct point forces[numMass])
{
	int i;
	for (i = 0; i < numMass; i++) {
		forces[i].z += -9.86;
		//checkForForceField(i, j, k, jello, forces);
	}
}

/*Checks and applies spring force if any of the masses have collided with the boundary planes*/
void computeCollisionSpringForces(struct world * jello, struct point forces[numMass])
{
	int i;
	for (i = 0; i < numMass; i++) {
		checkForCollision(i, jello, forces);
	}
}

/*	Calculates the force each mass-point has on its bend neighbors, and
	sums the total force at each point. Returns forces on each mass
	in result array 'forces'*/
void computeBendSpringForces(struct world * jello, struct point forces[numMass])
{
	int i;
	for (i = 0; i < numMass - 1; i++) {
		applyBendForceToNeighbors(i, jello, forces);
	}
}

/*	Calculates the force each mass-point has on its structural neighbors, and
	sums the total force at each point. Returns forces on each mass
	in result array 'forces'*/
void computeStructuralSpringForces(struct world * jello, struct point forces[numMass])
{
	int i;
	for (i = 0; i < numMass; i++) {
		applyStructuralForceToNeighbors(i, jello, forces);
	}
}

/*Takes the total force applied to each mass and computes the acceleration by dividing each component
by the mass. Function should be called after all shear, bend, structural, collision, and external forces
have been calculated and summed*/
void computeEachMassAcceleration(struct world * jello, struct point forces[numMass], struct point a[numMass])
{
	int i;

	for (i = 0; i < numMass; i++) {
		a[i].x = forces[i].x / jello->mass;
		a[i].y = forces[i].y / jello->mass;
		a[i].z = forces[i].z / jello->mass;
	}
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[numMass])
{
  /* for you to implement ... */
	struct point forces[numMass];
	int i;
	double alpha = 0.0;

	for (i = 0; i < numMass; i++) {
		forces[i].x = 0.0;
		forces[i].y = 0.0;
		forces[i].z = 0.0; 
		if (up == 1) {
			forces[i].z += 1.0 * ((double)i);
		}
		else if (down == 1) {
			forces[i].z += -1.0 * ((double)i);
		}

		if (left == 1) {
			forces[i].y += -1.0 * ((double)i);
		} 
		else if (right == 1) {
			forces[i].y += 1.0 * ((double)i);
		}
	}

	computeStructuralSpringForces(jello, forces);
	springSmoothingFunction(jello, alpha);
	ComputeFrames(jello);
	ComputeReferenceVectors(jello);
	computeBendSpringForces(jello, forces);
	
	//computeCollisionSpringForces(jello, forces);
	computeExternalForces(jello, forces);

	computeEachMassAcceleration(jello, forces, a);

}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j = 0,k = 0;
  point a[numMass];

  computeAcceleration(jello, a);
  
  //keep first mass frozen
  jello->p[0].x = jello->p_init[0].x;
  jello->p[0].y = jello->p_init[0].y;
  jello->p[0].z = jello->p_init[0].z;
  jello->v[0].x = 0;
  jello->v[0].y = 0;
  jello->v[0].z = 0;

  double maxLength = 3.0 / 7.0;
  double maxVelocity = 30;

  for (i=1; i < numMass; i++)
  {
	jello->p[i].x += jello->dt * jello->v[i].x;
	jello->p[i].y += jello->dt * jello->v[i].y;
	jello->p[i].z += jello->dt * jello->v[i].z;
	jello->v[i].x += jello->dt * a[i].x;
	jello->v[i].y += jello->dt * a[i].y;
	jello->v[i].z += jello->dt * a[i].z;
  }

  //LIMIT VELOCITY
  for (i = 1; i < numMass; i++) {
	  double velocityMagnitude = computeVectorMagnitude(jello->v[i]);
	  if (velocityMagnitude > maxVelocity) {
		  double ratio = maxVelocity / velocityMagnitude;
		  pMULTIPLY(jello->v[i], ratio, jello->v[i]);
	  }
  }

  //LIMIT POSITION
  for (i = 0; i < numMass - 1; i++) {
	  struct point lengthVector;
	  pDIFFERENCE(jello->p[i + 1], jello->p[i], lengthVector);

	  if (computeVectorMagnitude(lengthVector) > maxLength) {
		  struct point newPosition;
		  struct point offset;
		  double ratio = maxLength / computeVectorMagnitude(lengthVector);

		  pMULTIPLY(lengthVector, ratio, newPosition);
		  pSUM(jello->p[i], newPosition, newPosition);
		  pDIFFERENCE(newPosition, jello->p[i + 1], offset);
		  pCPY(newPosition, jello->p[i + 1]);

		  for (int j = i + 1; j < numMass - 1; j++) {
			  pSUM(jello->p[j], offset, jello->p[j]);
		  }
	  }
  }
}

point computeVectorNormalized(point vector)
{
	point ret;
	double length = sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
	ret.x = vector.x / length;
	ret.y = vector.y / length;
	ret.z = vector.z / length;
	return ret;

}

void ComputeFramesInitial(struct world* jello)
{
	// ADDED: calculate initial local frames by parallel transport of root
	point up;
	up.x = 0.0;
	up.y = 0.0;
	up.z = 1.0;

	for (int i = 0; i < numMass - 1; i++)
	{
		point start;
		start.x = jello->p_smoothed[i].x;
		start.y = jello->p_smoothed[i].y;
		start.z = jello->p_smoothed[i].z;

		point end;
		end.x = jello->p_smoothed[i + 1].x;
		end.y = jello->p_smoothed[i + 1].y;
		end.z = jello->p_smoothed[i + 1].z;

		point currUp;
		currUp.x = up.x;
		currUp.y = up.y;
		currUp.z = up.z;

		point pos;
		pos.x = start.x;
		pos.y = start.y;
		pos.z = start.z;

		point aim;
		pDIFFERENCE(start, end, aim);
		aim = computeVectorNormalized(aim);

		point cross;
		CROSSPRODUCTp(aim, up, cross);
		cross = computeVectorNormalized(cross);

		CROSSPRODUCTp(cross, aim, currUp);
		currUp = computeVectorNormalized(currUp);

		jello->local_frames_init[i][0][0] = aim.x;
		jello->local_frames_init[i][1][0] = aim.y;
		jello->local_frames_init[i][2][0] = aim.z;
		jello->local_frames_init[i][0][1] = currUp.x;
		jello->local_frames_init[i][1][1] = currUp.y;
		jello->local_frames_init[i][2][1] = currUp.z;
		jello->local_frames_init[i][0][2] = cross.x;
		jello->local_frames_init[i][1][2] = cross.y;
		jello->local_frames_init[i][2][2] = cross.z;
		

		up.x = currUp.x;
		up.y = currUp.y;
		up.z = currUp.z;
	}
}

void ComputeFrames(struct world* jello)
{
	// ADDED: calculate initial local frames by parallel transport of root
	point up;
	up.x = 0.0;
	up.y = 0.0;
	up.z = 1.0;

	for (int i = 0; i < numMass - 1; i++)
	{
		point start;
		start.x = jello->p_smoothed[i].x;
		start.y = jello->p_smoothed[i].y;
		start.z = jello->p_smoothed[i].z;

		point end;
		end.x = jello->p_smoothed[i + 1].x;
		end.y = jello->p_smoothed[i + 1].y;
		end.z = jello->p_smoothed[i + 1].z;

		point currUp;
		currUp.x = up.x;
		currUp.y = up.y;
		currUp.z = up.z;

		point pos;
		pos.x = start.x;
		pos.y = start.y;
		pos.z = start.z;

		point aim;
		pDIFFERENCE(start, end, aim);
		aim = computeVectorNormalized(aim);

		point cross;
		CROSSPRODUCTp(aim, up, cross);
		cross = computeVectorNormalized(cross);

		CROSSPRODUCTp(cross, aim, currUp);
		currUp = computeVectorNormalized(currUp);

		jello->local_frames_curr[i][0][0] = aim.x;
		jello->local_frames_curr[i][1][0] = aim.y;
		jello->local_frames_curr[i][2][0] = aim.z;
		jello->local_frames_curr[i][0][1] = currUp.x;
		jello->local_frames_curr[i][1][1] = currUp.y;
		jello->local_frames_curr[i][2][1] = currUp.z;
		jello->local_frames_curr[i][0][2] = cross.x;
		jello->local_frames_curr[i][1][2] = cross.y;
		jello->local_frames_curr[i][2][2] = cross.z;


		up.x = currUp.x;
		up.y = currUp.y;
		up.z = currUp.z;
	}
}

void multiplyFrameByEdgeInit(double matrix[][3], double edge[][1], point &dst) {
	double matrixT[3][3];
	transposeMatrix(matrix, matrixT);

	double temp[3][1];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			temp[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				temp[i][j] += matrix[i][k] * edge[k][j];
			}
		}
	}

	dst.x = temp[0][0];
	dst.y = temp[1][0];
	dst.z = temp[2][0];
}

void multiplyFrameByEdge(double matrix[][3], double edge[][1], point &dst) {
	double matrixT[3][3];

	double temp[3][1];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 1; j++) {
			temp[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				temp[i][j] += matrix[i][k] * edge[k][j];
			}
		}
	}

	dst.x = temp[0][0];
	dst.y = temp[1][0];
	dst.z = temp[2][0];
}

void transposeMatrix(double src[][3], double dst[][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			dst[j][i] = src[i][j];
		}
	}
}

void ComputeInitialReferenceVectors(struct world* jello) {

	for (int i = 1; i < numMass; i++) {
		point edge;
		pDIFFERENCE(jello->p_init[i], jello->p_init[i + 1], edge);
		double edgeAsArray[3][1];
		edgeAsArray[0][0] = edge.x;
		edgeAsArray[1][0] = edge.y;
		edgeAsArray[2][0] = edge.z;

		multiplyFrameByEdgeInit(jello->local_frames_init[i], edgeAsArray, jello->t_init[i]);
	}
}

void ComputeReferenceVectors(struct world* jello) {

	for (int i = 1; i < numMass; i++) {
		double edgeAsArray[3][1];
		edgeAsArray[0][0] = jello->t_init[i].x;
		edgeAsArray[1][0] = jello->t_init[i].y;
		edgeAsArray[2][0] = jello->t_init[i].z;

		multiplyFrameByEdge(jello->local_frames_curr[i], edgeAsArray, jello->t[i]);
	}
}
