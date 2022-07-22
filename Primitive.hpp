// Fall 2018

#pragma once

#include <glm/glm.hpp>
#include <vector>
//#include "Material.hpp"


struct Intersection 
{
  double t = 0;
  glm::vec3 p;
  glm::vec3 n;
  glm::vec3 kd;
  glm::vec3 ks;
  double kr;
  double s;
  int id = -1;
  bool burning = false;
  float temperature;
  int time;
};


struct Triangle
{
  size_t v1;
  size_t v2;
  size_t v3;
  bool burning = false;

  Triangle( size_t pv1, size_t pv2, size_t pv3 )
    : v1( pv1 )
    , v2( pv2 )
    , v3( pv3 )
  {}
};

class Primitive {
  
public:
  float epsilon = 0.00001;
  virtual ~Primitive();

  virtual Intersection Intersect(glm::vec3 orig, glm::vec3 dest);
  double IntersectTri(glm::vec3 orig, glm::vec3 dest, glm::vec3 V0, glm::vec3 V1, glm::vec3 V2);
  void setEp(float e) {epsilon = e;}
};

class Cylinder : public Primitive {
public:
  virtual ~Cylinder();
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;
};

class Cone : public Primitive {
public:
  virtual ~Cone();
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;
};

class Sphere : public Primitive {
public:
  virtual ~Sphere();
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;
};

class Cube : public Primitive {
public:
  virtual ~Cube();
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;
private:
  std::vector<glm::vec3> m_vertices = {};
  std::vector<Triangle> m_faces = {};
};

class NonhierSphere : public Primitive {
public:
  NonhierSphere(const glm::vec3& pos, double radius)
    : m_pos(pos), m_radius(radius)
  {
  }
  virtual ~NonhierSphere();
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;
 
private:
  glm::vec3 m_pos;
  double m_radius;
};

class NonhierBox : public Primitive {
public:
  NonhierBox(const glm::vec3& pos, double size)
    : m_pos(pos), m_size(size)
  {
  }
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;
  virtual ~NonhierBox();
  
private:
  glm::vec3 m_pos;
  double m_size;
  std::vector<glm::vec3> m_vertices = {};
  std::vector<Triangle> m_faces = {};
};
