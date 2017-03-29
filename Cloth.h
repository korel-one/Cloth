#ifndef CLOTH_H_
#define CLOTH_H_

#include "bimap.h"

#include <Vector3.h>
#include <deque>
#include <vector>
#include <map>
#include <memory>


namespace _Math {
  class Vector3;
}

class _Constraint;
struct Plane3D;
class ConstraintStretchUV;

class Cloth
{
public:
  Cloth(int rows
      , int colomns
      , float size
      , const _Math::Vector3& pointToBeginWith);

private:
  void init(int columns
      , int rows
      , float size
      , const _Math::Vector3& pointToBeginWith);

public:
  virtual ~Cloth();

  void CalculateNormals(void);

  const _Math::Vector3* GetVertices(void) const;
  const _Math::Vector3* GetNormals(void) const;
  const int* GetIndices(void) const;

  int GetIndicesCount(void) const; //remove candidate

  bool SetMassInversedForPartilce(int i, float massInversed);
  bool SetIsUnderForceForPartilce(int i, bool particleUnderForce);

  bool SetTimeStep(float dt);

  void Calculate(void);

  void SetFreeFallAcceleration(const _Math::Vector3& freeFallAcceleration);
  bool SetIterationsCount(int iterationsCount);

  bool AddPositionToPoint(int i, const _Math::Vector3& addToPosition);

  bool PickUpPoints(const _Math::Vector3& nearestPoint, const _Math::Vector3& furthestPoint, _Math::Vector3& pickedUpPoint, int* indicesOfPickedUpPoints, int countOfPickedUpPoints);

  void CreateTear();

private:
  void CalculatePredictedPositions();
  void ApplyConstraints();
  void CorrectVelocities();

  void CalculateAlphaBeta(const _Math::Vector3& point1, const _Math::Vector3& point2, const _Math::Vector3& point3, float* alpha, float* beta);

  void TearCloth(void);

private:
  //---------------------------------------------------------------
  class Tear;
  typedef std::shared_ptr<Tear> TearPtr;
  typedef bool TApply;
  typedef std::pair<TApply, TearPtr> TAlgorithmTear;
  TAlgorithmTear __tear;

  int __count;
  int __indicesCount; //remove candidate

  std::vector<_Math::Vector3> __vertices;
  std::vector<int> __indices;

  static const float THRESHOLD;
  //---------------------------------------------------------------

  std::vector<_Math::Vector3> __normals;

  std::vector<float> __massInverseParticles;
  std::deque<bool> __underForceParticles;

  float __dt;
  float __dt_1;
  int __iterationsCount;
  float __iterationsCountInversed;

  std::vector<_Math::Vector3> __velocities;

  _Math::Vector3 __freeFallAcceleration;
  std::vector<_Math::Vector3> __predictedPositions;

  std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> > __constraints;
};

#endif /* CLOTH_H_ */
