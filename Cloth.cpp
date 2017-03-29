#include "Cloth.h"
#include "UtilGeometry.h"

#include <MathVector3.h>
#include <MathVector4.h>
#include <MathMatrix4.h>
#include <MathMatrix3.h>
#include <IntConstr_ConstraintStretch.h>
#include <ConstraintStretchUV.h>
#include <ConstraintBendUV.h>
#include <IntConstr_BendingConstraint.h>
#include <IntConstr_ConstraintInternalStateAux.h>
#include <MathAdapterFunctions.h>
#include <MathPlane.h>

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <set>

#include <iostream>

namespace {
#include <time.h>
#include <string>
using std::string;
using std::cout;
clock_t t1, t2;
void
ST()
{
        t1 = clock();
}
void
ET(std::string s)
{
        t2 = clock();
        std::cout << (double)(t2 - t1) / (CLOCKS_PER_SEC / 1000) << "\t" << s  << "\n";
        //cout << (float)(t2 - t1) / (CLOCKS_PER_SEC / 1000) << "\t";
}
}//namespace

namespace details
{
//-----------------------------------------------------------------------
//print sequence
//-----------------------------------------------------------------------
template <
template < class T,
class Alloc = std::allocator<T>>
class Container
, class TValue>
void print(const Container<TValue>& container)
{
    std::cout<<"print sequence:\n";
    for(auto it = std::begin(container);
            it != std::end(container); ++it)
        std::cout<<" "<<*it;

    std::cout<<std::endl;
}

//-----------------------------------------------------------------------
//print set
//-----------------------------------------------------------------------
template<
template<class T
, class Compare = std::less<T>
, class Alloc = std::allocator<T>>
class Container
, class TValue>
void print(const Container<TValue>& container)
{
    std::cout<<"print set:\n";
    for(auto it = std::begin(container);
            it != std::end(container); ++it)
        std::cout<<" "<<*it;

    std::cout<<std::endl;
}

//-----------------------------------------------------------------------
//print map
//-----------------------------------------------------------------------
template<
template<class Key
, class T
, class Compare = std::less<Key>
, class Alloc = std::allocator<std::pair<const Key,T>>>
class Container
, class TKey
, class TValue>
void print(const Container<TKey, TValue>& container)//map
{
  for(auto it = std::begin(container); it != std::end(container); ++it)
    {
      std::cout<<"{";
      std::cout<<(it->first);//<typename Container::key_type>
      std::cout<<" :";

      print(it->second);//<typename Container::mapped_type>
      std::cout<<"}\n";
    }

  std::cout<<std::endl;
}
}

using namespace _Math;

const float Cloth::THRESHOLD = 1.5;

//-------------------------------------------------------------------------------
//class Cloth::Tear
//-------------------------------------------------------------------------------
class Cloth::Tear
 {
 public:
   Tear(int count);
   void Init(const std::vector<int>& indices
       , const std::vector<Vector3>& vertices);

   void Detect(
       const std::vector<Vector3>& vertices
       , const std::vector<int>& indices);

   std::map<int, std::vector<int>> Update(const std::vector<int>& indices);

   //by traversing neighbor triangles
   int MaxStretchAdjacentTo(int point
       , const std::vector<Vector3>& vertices
       , const std::vector<int>& indices) const;

   //accessors:
   std::set<int> GetAdjacentPointsTo(int point, const std::vector<int>& indices) const;

   bool EdgeProcessed(int point1, int point2) const;
   void MarkEdgeProcessed(int point1, int point2);

 private:
   bool SplitByPlane(const RuptureCalculation::Utility::Plane3D& splitPlane
       , int point
       , const std::vector<int>& indices
       , const std::vector<Vector3>& vertices);

 private:
   std::vector<std::set<int>> __neighborTriangles;
   std::vector<Vector3> __verticesInit;

   //accumulated data
   typedef int TearPoint;
   typedef std::vector<int> TApperTriangles;
   std::map<TearPoint, TApperTriangles> __tearingData;


   typedef int TOrigin;
   typedef int TMaxStretchAdjacent;
   utils::bimap<TOrigin, TMaxStretchAdjacent> __stretchConstraints;

   //util container;
   std::map<int, std::vector<int>> __dividedInto;
 };

Cloth::Tear::Tear(int count)
        : __neighborTriangles(count)
{
}

//-------------------------------------------------------------------------------
void
Cloth::Tear::Init(const std::vector<int>& indices
    , const std::vector<Vector3>& vertices)
{
    __verticesInit = vertices;

    int N = __neighborTriangles.size();
    int indicesCount = indices.size();

    for(int i = 0; i != N; ++i)
    {
        for(int j = 0; j != indicesCount; j += 3)
        {
            if(i == indices[j] ||
                    i == indices[j + 1] ||
                    i == indices[j + 2])
            {
                __neighborTriangles[i].insert(j);
            }
        }
    }
}

//-------------------------------------------------------------------------------
std::set<int>
Cloth::Tear::GetAdjacentPointsTo(int point, const std::vector<int>& indices) const
{
    auto& neighbors = __neighborTriangles[point];
    std::set<int> adjacentPoints;

    for(auto it = std::begin(neighbors); it != std::end(neighbors); ++it)
    {
        Triangle triangle(indices[*it], indices[*it+1], indices[*it + 2]);
        Triangle::TCoupleAdjacent points = triangle.GetAdjacentTo(point);

        adjacentPoints.insert(indices[*it + points.first]);
        adjacentPoints.insert(indices[*it + points.second]);
    }

    return adjacentPoints;
}

//-------------------------------------------------------------------------------
int
Cloth::Tear::MaxStretchAdjacentTo(int point
    , const std::vector<Vector3>& vertices
    , const std::vector<int>& indices) const
{
  const std::set<int>& _adjacentPoints = GetAdjacentPointsTo(point, indices);

  float maxStretchRatio = 0.f;
  int maxStretchAdjacentPoint = Triangle::NOT_POINT;

  std::for_each(std::begin(_adjacentPoints), std::end(_adjacentPoints)
  , [&](int adjacentPoint)
  {
    float deformationDistance = vertices[point].Distance(vertices[adjacentPoint]);
    float initDistance = __verticesInit[point].Distance(__verticesInit[adjacentPoint]);

    float currRatio = deformationDistance/initDistance;
    if(currRatio > maxStretchRatio && currRatio > THRESHOLD)
      {
        maxStretchRatio = currRatio;
        maxStretchAdjacentPoint = adjacentPoint;
      }
  });

  return maxStretchAdjacentPoint;
}

//TEST------------------------------------------
int foo(const std::vector<int>& indices
    , const std::vector<Vector3>& vertices
    , const Plane3D& splitPlane
    , int startIndex, int pointFrom)
  {

  if(pointFrom == indices[startIndex])
    {
    //TCoupleAdjacent(1, 2)
    auto& p1 = vertices[indices[startIndex + 1]];
    auto& p2 = vertices[indices[startIndex + 2]];

    return splitPlane.SpaceRelationTo(p1, p2);
    }
  else if(pointFrom == indices[startIndex+1])
    {
    //TCoupleAdjacent(0, 2)
    auto& p1 = vertices[indices[startIndex]];
    auto& p2 = vertices[indices[startIndex + 2]];

    return splitPlane.SpaceRelationTo(p1, p2);
    }
  else //if(pointFrom == indices[startIndex+2])
    {
    //TCoupleAdjacent(0, 1)
    auto& p1 = vertices[indices[startIndex]];
    auto& p2 = vertices[indices[startIndex + 1]];

    return splitPlane.SpaceRelationTo(p1, p2);
    }
  }

bool
Cloth::Tear::SplitByPlane(const Plane3D& splitPlane
                          , int point
                          , const std::vector<int>& indices
                          , const std::vector<Vector3>& vertices)
{
    auto& neighborTriangles = __neighborTriangles[point];
    __dividedInto.clear();

    for(auto it = std::begin(neighborTriangles); it != std::end(neighborTriangles); ++it)
    {
        int spaceRelation = foo(indices, vertices, splitPlane, *it, point);

        //TODO: handle triangles on plane?
        __dividedInto[spaceRelation].push_back(*it);
    }

    auto it_above = __dividedInto.find(Plane3D::ABOVE_PLANE);
    auto it_below = __dividedInto.find(Plane3D::BELOW_PLANE);

    //accumulate only triangles above a plane
    if(it_above != std::end(__dividedInto) &&
            it_below != std::end(__dividedInto))
    {
        __tearingData[point] = it_above->second;
        return true;
    }

    return false;
}

//-------------------------------------------------------------------------------
void
Cloth::Tear::Detect(const std::vector<Vector3>& vertices
                    , const std::vector<int>& indices)
{
    __stretchConstraints.clear();
    __tearingData.clear();

    int N = vertices.size();
    for(int i = 0; i != N; ++i)
    {
        int maxStretchAdjacentPoint = MaxStretchAdjacentTo(i, vertices, indices);

        if(Triangle::NOT_POINT != maxStretchAdjacentPoint
                && !EdgeProcessed(i, maxStretchAdjacentPoint))
        {
            //calculate rupture points
            Plane3D plane(vertices[i], vertices[maxStretchAdjacentPoint] - vertices[i]);

            if (SplitByPlane(plane, i, indices, vertices))
            {
                //uniqueness constraint
                MarkEdgeProcessed(i, maxStretchAdjacentPoint);
            }

        }
    }
}

//-------------------------------------------------------------------------------
std::map<int, std::vector<int>>
Cloth::Tear::Update(const std::vector<int>& indices)
{
    int N = __neighborTriangles.size();
    int initVerticesCount = __verticesInit.size();

    __neighborTriangles.resize(N + __tearingData.size());

    std::map<int, std::vector<int>> delta;

    std::for_each(std::begin(__tearingData), std::end(__tearingData)
                  , [=, &N, &initVerticesCount, &delta](const typename decltype(__tearingData)::value_type& tear)
                  {
                      __verticesInit.resize(++initVerticesCount);
                      __verticesInit[initVerticesCount-1] = __verticesInit[tear.first];

                      //apper triangles
                      for(auto it = std::begin(tear.second); it != std::end(tear.second); ++it)
                      {
                          __neighborTriangles[tear.first].erase(*it);
                          __neighborTriangles[N].insert(*it);

                          Triangle triangle(indices[*it], indices[*it + 1], indices[*it + 2]);
                          int originShift = triangle.GetIndex(tear.first);
                          delta[tear.first].push_back(*it + originShift);
                      }
                      ++N;
                  }
                 );

    return delta;
}

//-------------------------------------------------------------------------------
bool
Cloth::Tear::EdgeProcessed(int point1, int point2) const
{
    try
    {
        return (__stretchConstraints.direct(point1) == point2) &&
        (__stretchConstraints.inverse(point1) == point2);
    }
    catch(...)//std::string or std::exception
    {
        return false;
    }
}

//-------------------------------------------------------------------------------
void
Cloth::Tear::MarkEdgeProcessed(int point1, int point2)
{
  if(!EdgeProcessed(point1, point2))
    {
      __stretchConstraints.insert(point1, point2);
    }
}

//-------------------------------------------------------------------------------
//class Cloth
//-------------------------------------------------------------------------------
Cloth::Cloth(int rows, int colomns, float size
             , const Vector3& pointToBeginWith)
        : __count(rows * colomns)
        , __indicesCount(2*3*(rows-1)*(colomns-1))
        , __vertices(__count)
        , __indices(__indicesCount)
        , __normals(__count)
        , __massInverseParticles(__count, 1.f)
        , __underForceParticles(__count, true)
        , __dt(0.01f)
        , __dt_1(1.f / __dt)
        , __iterationsCount(1)
        , __iterationsCountInversed(1.f)
        , __velocities(__count)
        , __freeFallAcceleration(0.f, 0.f, 1.f)
        , __predictedPositions(__count)
{
    init(colomns, rows, size, pointToBeginWith);
    //__tear = std::make_pair<TApply, TearPtr>(false, std::make_shared<Tear>(colomns, rows));
    //__tear.second->Init(__indices, __vertices);
}

//-------------------------------------------------------------------------------
void Cloth::CreateTear()
  {
  if(__tear.second)
    return;

  __tear = std::make_pair<TApply, TearPtr>(false, std::make_shared<Tear>(__count));
  __tear.second->Init(__indices, __predictedPositions);//__vertices);
  }

//-------------------------------------------------------------------------------
void
Cloth::init(int columns, int rows, float size, const Vector3& pointToBeginWith)
{
    //indices
    int index = 0;
    for(int i = 0; i != columns; ++i)
    {
        for(int j = 0; j != rows; ++j)
        {
            int currStartTriangleIndex = i*columns + j;

            __vertices[currStartTriangleIndex] = pointToBeginWith + Vector3(j * size, i * size, 0.f);
            __predictedPositions[currStartTriangleIndex] = __vertices[currStartTriangleIndex];


            if(i != columns-1 && j != rows-1)
            {
                __indices[index++] = currStartTriangleIndex;
                __indices[index++] = currStartTriangleIndex + 1;
                __indices[index++] = currStartTriangleIndex + columns;

                __indices[index++] = currStartTriangleIndex + 1;
                __indices[index++] = currStartTriangleIndex + columns + 1;
                __indices[index++] = currStartTriangleIndex + columns;
            }
        }
    }


    for (int i = 0; i < __indicesCount; i += 3)
    {
        Vector3 p1 = __predictedPositions[__indices[i]] + Vector3(10.f, 15.f, 23.f);
        Vector3 p2 = __predictedPositions[__indices[i + 1]] + Vector3(10.f, 15.f, 23.f);
        Vector3 p3 = __predictedPositions[__indices[i + 2]] + Vector3(10.f, 15.f, 23.f);

        float alpha[3];
        float beta[3];
        CalculateAlphaBeta(p1, p2, p3, alpha, beta);

        _ConstraintInternalStateAux::NumbersAux numbersAux;
        numbersAux.push_back(1.f);
        numbersAux.push_back(alpha[0]);
        numbersAux.push_back(alpha[1]);
        numbersAux.push_back(alpha[2]);
        numbersAux.push_back(__indices[i]);
        numbersAux.push_back(__indices[i + 1]);
        numbersAux.push_back(__indices[i + 2]);
        _ConstraintInternalStateAux internalStateAux(numbersAux);


        ConstraintStretchUV* pConstraintAlpha =
        new ConstraintStretchUV(3, 1.f, _ConstraintInternalStateAux(numbersAux));

        internalStateAux.numbersAux[1] = beta[0];
        internalStateAux.numbersAux[2] = beta[1];
        internalStateAux.numbersAux[3] = beta[2];
        ConstraintStretchUV* pConstraintBeta = new ConstraintStretchUV(3, 1.f, internalStateAux);

        __constraints.insert(std::make_pair(i, std::make_pair(pConstraintAlpha, pConstraintBeta)));

        internalStateAux.numbersAux.clear();
        internalStateAux.numbersAux.push_back(alpha[0]);
        internalStateAux.numbersAux.push_back(alpha[1]);
        internalStateAux.numbersAux.push_back(alpha[2]);
        internalStateAux.numbersAux.push_back(beta[0]);
        internalStateAux.numbersAux.push_back(beta[1]);
        internalStateAux.numbersAux.push_back(beta[2]);
    }

}

Cloth::~Cloth()
{
    std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator it = __constraints.begin();
    std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator itEnd = __constraints.end();
    while (it != itEnd)
    {
        delete it->second.first;
        delete it->second.second;
        ++it;
    }
}

const Vector3*
Cloth::GetVertices(void) const
{
    return &__vertices.at(0);
}

const int*
Cloth::GetIndices(void) const
{
    return &(__indices.at(0));
}

int
Cloth::GetIndicesCount(void) const
{
    return __indicesCount;
}

const Vector3*
Cloth::GetNormals(void) const
{
    return &__normals.at(0);
}

void
Cloth::CalculateNormals(void)
{
    for (int i = 0; i < __count; ++i)
    {
        __normals[i].Set(0.f);
    }

    //iterate over triangles
    for (int i = 0; i < __indicesCount; i += 3)
    {
        const Vector3& v1 = __vertices[__indices[i]];
        const Vector3& v2 = __vertices[__indices[i + 1]];
        const Vector3& v3 = __vertices[__indices[i + 2]];

        Vector3 normal = (v2 - v1).Cross(v3 - v1);
        normal.Normalize();

        __normals[__indices[i]] += normal;
        __normals[__indices[i + 1]] += normal;
        __normals[__indices[i + 2]] += normal;
    }

    for (int i = 0; i < __count; ++i)
        __normals[i].Normalize();
}

bool
Cloth::SetMassInversedForPartilce(int i, float massInversed)
{
    if (i < 0 || i >= __count || massInversed < 0.f)
    {
        return false;
    }
    __massInverseParticles[i] = massInversed;
    return true;
}

bool
Cloth::SetIsUnderForceForPartilce(int i, bool particleUnderForce)
{
    if (i < 0 || i >= __count)
    {
        return false;
    }
    __underForceParticles[i] = particleUnderForce;
    return true;
}

bool
Cloth::SetTimeStep(float dt)
{
    if (dt < 0.f || Equal(dt, 0.f))
    {
        return false;
    }
    __dt = dt;
    __dt_1 = 1.f / dt;
    return true;
}

void
Cloth::Calculate(void)
{
  START_PROFILER_AV(CalculateApplyVeloc, 100)
    CalculatePredictedPositions();
    ApplyConstraints();
    CorrectVelocities();
  END_PROFILER_AV(CalculateApplyVeloc)
    TearCloth();
}

void
Cloth::CalculatePredictedPositions()
{
    for (int i = 0; i < __count; ++i)
    {
        if (__underForceParticles[i] == false)
        {
            continue;
        }
        Vector3& position = __vertices[i];
        Vector3& velocity = __velocities[i];
        velocity += __freeFallAcceleration * __dt;
        velocity *= 0.993;
        __predictedPositions[i] = position + velocity * __dt;
    }
    return;
}

void
Cloth::ApplyConstraints()
{
    for (int i = 0; i < __iterationsCount; ++i)
    {
        std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator it = __constraints.begin();
        std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator itEnd = __constraints.end();
        while (it != itEnd)
        {
            it->second.first->Apply(&__predictedPositions[0], &__massInverseParticles[0]);
            it->second.second->Apply(&__predictedPositions[0], &__massInverseParticles[0]);
            ++it;
        }
    }
}

void
Cloth::CorrectVelocities()
{
    for (int i = 0; i < __count; ++i)
    {
        if (__underForceParticles[i] == false)
        {
            continue;
        }
        __velocities[i] = (__predictedPositions[i] - __vertices[i]) * __dt_1;
        __vertices[i] = __predictedPositions[i];
    }
    return;
}

void
Cloth::SetFreeFallAcceleration(const Vector3& freeFallAcceleration)
{
    __freeFallAcceleration = freeFallAcceleration;
}

bool
Cloth::SetIterationsCount(int iterationsCount)
{
    if (iterationsCount <= 0)
    {
        return false;
    }
    __iterationsCount = iterationsCount;
    __iterationsCountInversed = 1.f / iterationsCount;

    std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator it = __constraints.begin();
    std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator itEnd = __constraints.end();
    while (it != itEnd)
    {
        (it)->second.first->Prepare(__iterationsCountInversed);
        (it)->second.second->Prepare(__iterationsCountInversed);
        ++it;
    }

    return true;
}

bool
Cloth::AddPositionToPoint(int i, const Vector3& addToPosition)
{
    if (i < 0 || i >= __count)
    {
        return false;
    }
    __vertices[i] += addToPosition;
    __predictedPositions[i] += addToPosition;
    return true;
}

bool
Cloth::PickUpPoints(const Vector3& nearestPoint, const Vector3& furthestPoint, Vector3& pickedUpPoint, int* indicesOfPickedUpPoints, int countOfPickedUpPoints)
{
  bool foundIntersection = false;

  if (countOfPickedUpPoints != 1 && countOfPickedUpPoints != 3)
  {
      __tear.first = foundIntersection;
      return foundIntersection;
  }

    pickedUpPoint = furthestPoint;
    Vector3 vecDir = furthestPoint - nearestPoint;
    for (int i = 0; i < __indicesCount; i += 3)
    {
        const Vector3& v1 = __vertices[__indices[i]];
        const Vector3& v2 = __vertices[__indices[i + 1]];
        const Vector3& v3 = __vertices[__indices[i + 2]];

        Vector3 v(v2 - v1);
        Vector3 w(v3 - v1);
        Vector3 normal(v.GetCrossed(w));
        if (normal.LengthSqr() < fEPS)
        {
            continue;
        }

        float dPlane = -normal.Dot(v1);

        float nv = normal.Dot(vecDir);

        if (fabsf(nv) < fEPS)
        {
            continue;
        }

        float s = -(dPlane + normal.Dot(nearestPoint)) / nv;

        if (s < 0)
        {
            continue;
        }

        //2. determine does ray intersect triangle (find out variables c and d)
        Vector3 q(nearestPoint + s * vecDir);
        Vector3 t(q - v1);

        float c;
        float d;
        if (!Plane::FindCoefficientsFoPointInPlane(v, w, t, c, d))
        {
            continue;
        }
        if (c >= 0.f && d >= 0.f && c + d <= 1.f)
        {
            if ((pickedUpPoint - nearestPoint).Length() > (q - nearestPoint).Length())
            {
                foundIntersection = true;
                pickedUpPoint = q;
                if (countOfPickedUpPoints == 1)
                {
                    indicesOfPickedUpPoints[0] = __indices[i];
                }
                else if (countOfPickedUpPoints == 3)
                {
                    indicesOfPickedUpPoints[0] = __indices[i];
                    indicesOfPickedUpPoints[1] = __indices[i + 1];
                    indicesOfPickedUpPoints[2] = __indices[i + 2];
                }
            }
        }
    }

    __tear.first = foundIntersection;
    return foundIntersection;
}

//-------------------------------------------------------------------------------
void
Cloth::TearCloth()
{
  if(!__tear.first)
    return;

  START_PROFILER_AV(Detect, 100)
    __tear.second->Detect(__vertices, __indices);
  END_PROFILER_AV(Detect)

//  START_PROFILER_AV(Update, 100)
    auto delta = __tear.second->Update(__indices);
//  END_PROFILER_AV(Update)

    for(auto it = std::begin(delta); it != std::end(delta); ++it)
    {
//      START_PROFILER_AV(ApplyDelta, 100)
        __vertices.resize(__count + 1);
        __predictedPositions.resize(__count + 1);
        __normals.resize(__count + 1);
        __velocities.resize(__count + 1);
        __massInverseParticles.resize(__count + 1);
        __underForceParticles.resize(__count + 1);

        __vertices[__count] = __vertices[it->first];
        __predictedPositions[__count] = __predictedPositions[it->first];
        __velocities[__count] = __velocities[it->first];
        __massInverseParticles[__count] = 1.f; //default;
        __underForceParticles[__count] = true;
        //__pNormal[__count] = __pNormal[i];

        for(auto it_indices = std::begin(it->second); it_indices != std::end(it->second)
                ; ++it_indices)
        {
            int index = *it_indices;

            std::map<int, std::pair<ConstraintStretchUV*, ConstraintStretchUV*> >::iterator it_res = __constraints.find(index - index % 3);
            assert(it_res != __constraints.end());

            //          const int* pOldIndices = it_res->second.first->GetIndices();

            //          if (__pIndices[index] == pOldIndices[0])
            if (index % 3 == 0)
            {
                it_res->second.first->SetNewIndex(0, __count);
                it_res->second.second->SetNewIndex(0, __count);
            }
            //          else if (__pIndices[index] == pOldIndices[1])
            else if (index % 3 == 1)
            {
                it_res->second.first->SetNewIndex(1, __count);
                it_res->second.second->SetNewIndex(1, __count);
            }
            else if (index % 3 == 2)
            {
                it_res->second.first->SetNewIndex(2, __count);
                it_res->second.second->SetNewIndex(2, __count);
            }
            else
            {
                assert(false);
            }
            __indices[index] = __count;
        }



        ++__count;

//        END_PROFILER_AV(ApplyDelta)
    }
}

void
Cloth::CalculateAlphaBeta(const Vector3& point1, const Vector3& point2, const Vector3& point3, float* alpha, float* beta)
{
    float d1 = (point2 - point3).Dot(point1);
    float d2 = (point2 - point3).Dot(point2);
    float d3 = (point2 - point3).Dot(point3);

    float r = d3 * d3 * (point1 - point2).Dot(point1 - point2);
    r += 2 * d1 * d3 * (point1 - point2).Dot(point2 - point3);
    r += d2 * d2 * (point1 - point3).Dot(point1 - point3);
    r -= 2 * d2 * (d3 * (point1 - point2).Dot(point1 - point3) + d1 * (point1 - point3).Dot(point2 - point3));
    r += d1 * d1 * (point2 - point3).Dot(point2 - point3);
    r = sqrtf(r);

    //	float alpha[3], beta[3];
    alpha[0] = 0.f;
    alpha[1] = 1.f / (point2 - point3).Length();
    alpha[2] = -alpha[1];
    beta[0] = (d2 - d3) / r;
    beta[1] = (d3 - d1) / r;
    beta[2] = (d1 - d2) / r;


    Vector3 u = alpha[0] * point1 + alpha[1] * point2 + alpha[2] * point3;
    Vector3 v = beta[0] * point1 + beta[1] * point2 + beta[2] * point3;
    float maxMin = 0.f;
    //	float maxAngle = 0.f;
    for (int j = 30; j < 90; ++j)
    {
        float angle = j * M_PI / 180.f;
        Matrix4 rotation;
        rotation = Matrix4::CreateRotationAxis((point2 - point3).Cross(point1 - point3), angle);
        Vector3 u1 = rotation * u;
        Vector3 v1 = rotation * v;

        float delta = Matrix3(point1.x, point2.x, point3.x,
                              point1.y, point2.y, point3.y,
                              point1.z, point2.z, point3.z).Determinant();

        float deltaU1 = Matrix3(u1.x, point2.x, point3.x,
                                u1.y, point2.y, point3.y,
                                u1.z, point2.z, point3.z).Determinant();
        float deltaU2 = Matrix3(point1.x, u1.x, point3.x,
                                point1.y, u1.y, point3.y,
                                point1.z, u1.z, point3.z).Determinant();
        float deltaU3 = Matrix3(point1.x, point2.x, u1.x,
                                point1.y, point2.y, u1.y,
                                point1.z, point2.z, u1.z).Determinant();

        float deltaV1 = Matrix3(v1.x, point2.x, point3.x,
                                v1.y, point2.y, point3.y,
                                v1.z, point2.z, point3.z).Determinant();
        float deltaV2 = Matrix3(point1.x, v1.x, point3.x,
                                point1.y, v1.y, point3.y,
                                point1.z, v1.z, point3.z).Determinant();
        float deltaV3 = Matrix3(point1.x, point2.x, v1.x,
                                point1.y, point2.y, v1.y,
                                point1.z, point2.z, v1.z).Determinant();

        float alpha1[3];
        float beta1[3];
        alpha1[0] = deltaU1 / delta;
        alpha1[1] = deltaU2 / delta;
        alpha1[2] = deltaU3 / delta;
        beta1[0] = deltaV1 / delta;
        beta1[1] = deltaV2 / delta;
        beta1[2] = deltaV3 / delta;

        float minCur = fabsf(alpha1[0]);
        if (fabsf(alpha1[1]) < minCur)
            minCur = fabsf(alpha1[1]);
        if (fabsf(alpha1[2]) < minCur)
            minCur = fabsf(alpha1[2]);
        if (fabsf(beta1[0]) < minCur)
            minCur = fabsf(beta1[0]);
        if (fabsf(beta1[1]) < minCur)
            minCur = fabsf(beta1[1]);
        if (fabsf(beta1[2]) < minCur)
            minCur = fabsf(beta1[2]);

        if (minCur > maxMin)
        {
            //			maxAngle = angle;

            maxMin = minCur;

            alpha[0] = alpha1[0];
            alpha[1] = alpha1[1];
            alpha[2] = alpha1[2];

            beta[0] = beta1[0];
            beta[1] = beta1[1];
            beta[2] = beta1[2];
        }
    }
    return;
}
