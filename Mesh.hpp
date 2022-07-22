// Fall 2018

#pragma once

#include <vector>
#include <iosfwd>
#include <string>

#include <glm/glm.hpp>

#include "Primitive.hpp"
/*
struct Triangle
{
	size_t v1;
	size_t v2;
	size_t v3;

	Triangle( size_t pv1, size_t pv2, size_t pv3 )
		: v1( pv1 )
		, v2( pv2 )
		, v3( pv3 )
	{}
	};*/

struct Tuple
{
  int t1;
  int t2;

  Tuple( int pt1, int pt2 )
    : t1( pt1 )
    , t2( pt2 )
  {}

};


// A polygonal mesh.
class Mesh : public Primitive {
public:
  Mesh() {};
  Mesh( const std::string& fname );

  bool IntersectBound(glm::vec3 orig, glm::vec3 dest);
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;  
  int get_verts() {return m_vertices.size();};

  void IFS(int level, float scale);

  //private:
protected:
	std::vector<glm::vec3> m_vertices;
	std::vector<Triangle> m_faces;
  std::vector<glm::vec3> m_vnorms;
  std::vector<int> m_facecount;
  glm::vec3 m_min;
  glm::vec3 m_max;
 
    friend std::ostream& operator<<(std::ostream& out, const Mesh& mesh);
};
