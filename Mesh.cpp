// Fall 2018

#include <iostream>
#include <fstream>

#include <glm/ext.hpp>
#include <map>
#include <string>
#include <random>
#include <stdlib.h> 
#include <time.h>    
// #include "cs488-framework/ObjFileDecoder.hpp"
#include "Mesh.hpp"


Mesh::Mesh( const std::string& fname )
	: m_vertices()
	, m_faces()
{

  //std::cout << "Creating mesh" << std::endl;
  //std::cout << fname << std::endl;
	std::string code;
	double vx, vy, vz;
	size_t s1, s2, s3;

	std::ifstream ifs( fname.c_str() );
	while( ifs >> code ) {
	  //std::cout << "Reading File" << std::endl; 
		if( code == "v" ) {
			ifs >> vx >> vy >> vz;
			m_vertices.push_back( glm::vec3( vx, vy, vz ) );
		} else if( code == "f" ) {
			ifs >> s1 >> s2 >> s3;
			//std::cout << "Adding a face" << std::endl;
			m_faces.push_back( Triangle( s1 - 1, s2 - 1, s3 - 1 ) );
		}
	}

	if (m_vertices.size() != 0) {	
	glm::vec3 cur_max = m_vertices[0];
	glm::vec3 cur_min = m_vertices[0];
	for (glm::vec3 v : m_vertices) {
	  if (v[0] > cur_max[0]) {cur_max[0] = v[0];}
	  if (v[1] > cur_max[1]) {cur_max[1] = v[1];}
	  if (v[2] > cur_max[2]) {cur_max[2] = v[2];}
	  if (v[0] < cur_min[0]) {cur_min[0] = v[0];}
	  if (v[1] < cur_min[1]) {cur_min[1] = v[1];}
	  if (v[2] < cur_min[2]) {cur_min[2] = v[2];}
	}
	m_max = cur_max;
	m_min = cur_min;
	}
}

void Mesh::IFS(int level, float scale) {
  for (int i = 0; i < 3; ++i) {
    m_vnorms.push_back(glm::vec3(0));
    m_facecount.push_back(0);
  }

  for (Triangle tri : m_faces) {

  }

  for (int i = 0; i < level; ++i) {  
    std::vector<Triangle> new_faces;

    for (Triangle tri : m_faces) {
        
      glm::vec3 V1 = m_vertices[tri.v1];
      glm::vec3 V2 = m_vertices[tri.v2];
      glm::vec3 V3 = m_vertices[tri.v3];

      glm::vec3 V4;
      glm::vec3 V5;
      glm::vec3 V6;

      int cur_vert = get_verts();
      int t1 = cur_vert;
      int t2 = cur_vert + 1;
      int t3 = cur_vert + 2;
      

      V4 = (V1 + V2)/2;
      V4[1] += (float)(scale/(i+1));
      bool changed = false;

      if (V4[1] > m_max[1]) { m_max[1] = V4[1]; }
      for (int i = 0; i < m_vertices.size(); ++i) {
	if (V4[0] == m_vertices[i][0] && V4[2] == m_vertices[i][2]) {
	  t1 = i;
	  changed = true;
	  --t2; 
	  --t3;
	  break;
	}
      } 
      if (!changed) {//(t1 == cur_vert) {
      m_vertices.push_back(V4);
      m_vnorms.push_back(glm::vec3(0));
      m_facecount.push_back(0);  
      }

      changed = false;

      V5 = (V2 + V3)/2;
      V5[1] += (float)(scale/(i+1));
        
      if (V5[1] > m_max[1]) { m_max[1] =V5[1]; }
      for (int i = 0; i < m_vertices.size(); ++i) {
        if (V5[0] == m_vertices[i][0] && V5[2] == m_vertices[i][2]) {
          t2 = i;
	  changed = true;
	  --t3;
          //break;
        }
      }
      if (!changed) { //(t2 == cur_vert + 1) {
      m_vertices.push_back(V5);
      m_vnorms.push_back(glm::vec3(0));
      m_facecount.push_back(0);        
      }

      changed = false;

      V6 = (V3 + V1)/2;
      V6[1] += (float)(scale/(i+1));      

      if (V6[1] > m_max[1]) { m_max[1] =V6[1]; }
      for (int i = 0; i < m_vertices.size(); ++i) {
        if (V6[0] == m_vertices[i][0] && V6[2] == m_vertices[i][2]) {
          t3 = i;
	  changed = true;
          //break;
        }
      }
      if (!changed) { //(t3 == cur_vert+2) {
      m_vertices.push_back(V6);
      m_vnorms.push_back(glm::vec3(0));
      m_facecount.push_back(0);     
      }

      new_faces.push_back( Triangle( tri.v1, t1 , t3) );
      new_faces.push_back( Triangle( t1, tri.v2 , t2) );
      new_faces.push_back( Triangle( t3, t2 , tri.v3) );;
      new_faces.push_back( Triangle( t1, t2 , t3) );
    }
    m_faces = new_faces;
  }

  for (Triangle tri : m_faces) {
    glm::vec3 n = glm::normalize(glm::cross(m_vertices[tri.v2] - m_vertices[tri.v1]
					    , m_vertices[tri.v3] - m_vertices[tri.v1]));
    m_vnorms[tri.v1] += n;
    ++m_facecount[tri.v1]; 
    m_vnorms[tri.v2] += n;
    ++m_facecount[tri.v2];
    m_vnorms[tri.v3] += n;
    ++m_facecount[tri.v3];
  }

  for (int i = 0; i < m_vertices.size(); ++i) {

    m_vnorms[i] /= m_facecount[i];

  }

  //std::cout << m_vertices.size() << std::endl;
  //std::cout << m_vnorms.size() << std::endl;
  //std::cout << m_facecount.size() << std::endl;
}


std::ostream& operator<<(std::ostream& out, const Mesh& mesh)
{
  out << "mesh {";
  /*
  
  for( size_t idx = 0; idx < mesh.m_verts.size(); ++idx ) {
  	const MeshVertex& v = mesh.m_verts[idx];
  	out << glm::to_string( v.m_position );
	if( mesh.m_have_norm ) {
  	  out << " / " << glm::to_string( v.m_normal );
	}
	if( mesh.m_have_uv ) {
  	  out << " / " << glm::to_string( v.m_uv );
	}
  }

*/
  out << "}";
  return out;
}


bool Mesh::IntersectBound(glm::vec3 orig, glm::vec3 dest) {
  glm::vec3 tmin = (m_min - orig) / dest;
  glm::vec3 tmax = (m_max - orig) / dest;
   
  if (tmin[0] > tmax[0]) { 
    double temp = tmax[0];
    tmax[0] = tmin[0];
    tmin[0] = temp;
  }

  if (tmin[1] > tmax[1]) {
    double temp = tmax[1];
    tmax[1] = tmin[1];
    tmin[1] = temp;
  }

  if ((tmin[0] > tmax[1]) || (tmin[1] > tmax[0])) {return false;}

  if (tmin[1] > tmin[0]) {tmin[0] = tmin[1];}
  if (tmax[0] > tmax[1]) {tmax[0] = tmax[1];}

  if (tmin[2] > tmax[2]) {
    double temp = tmax[2];
    tmax[2] = tmin[2];
    tmin[2] = temp;
  }

  if ((tmin[0] > tmax[2]) || (tmin[2] > tmax[0])) {return false;}

  return true;
}


Intersection Mesh::Intersect(glm::vec3 orig, glm::vec3 dest) {
  Intersection I; //= Intersection();
  
  /*
  glm::vec3 m_pos = glm::vec3(0,0,0);
  double m_size = 1;
  std::vector<glm::vec3> m_vertices;
  std::vector<Triangle> m_faces;

  m_pos = m_min;
  m_size = m_max[0] - m_min[0];

  if (m_vertices.size() == 0) {
    m_vertices.push_back(m_pos);
    m_vertices.push_back(m_pos + glm::vec3(m_size, 0, 0));
    m_vertices.push_back(m_pos + glm::vec3(0, m_size, 0));
    m_vertices.push_back(m_pos + glm::vec3(m_size, m_size, 0));
    m_vertices.push_back(m_pos + glm::vec3(0, 0, m_size));
    m_vertices.push_back(m_pos + glm::vec3(m_size, 0, m_size));
    m_vertices.push_back(m_pos + glm::vec3(0, m_size, m_size));
    m_vertices.push_back(m_pos + glm::vec3(m_size, m_size, m_size));
  }

  if (m_faces.size() == 0) {
    m_faces.push_back(Triangle (0, 1, 2));
    m_faces.push_back(Triangle (2, 1, 3));

    m_faces.push_back(Triangle (5, 1, 7));
    m_faces.push_back(Triangle (7, 1, 3));

    m_faces.push_back(Triangle (4, 5, 6));
    m_faces.push_back(Triangle (6, 5, 7));

    m_faces.push_back(Triangle (0, 4, 2));
    m_faces.push_back(Triangle (2, 4, 6));


    m_faces.push_back(Triangle (6, 7, 2));
    m_faces.push_back(Triangle (2, 7, 3));

    m_faces.push_back(Triangle (4, 5, 0));
    m_faces.push_back(Triangle (0, 5, 1));
  }

  //Intersection I;

  glm::vec3 n;
  double cur_t = 0;
  for (Triangle tri : m_faces) {
    double t = IntersectTri(orig, dest,
                            m_vertices[tri.v1], m_vertices[tri.v2], m_vertices[tri.v3]);
    if (t != 0) {
      if ((cur_t == 0) || (t < cur_t)) {
        cur_t = t;
        n = glm::normalize(glm::cross(m_vertices[tri.v2] - m_vertices[tri.v1],
				      m_vertices[tri.v3]- m_vertices[tri.v1]));
      }
    }
  }
  I.t = cur_t;
  I.p = orig + (float)I.t * dest;
  I.n = n;
  return I;

  */





  if (!IntersectBound(orig, dest)) { return I;}
  //std::cout << "Calling Intersect" << std::endl;
  double cur_t = 0;
  glm::vec3 n;
  //std::cout << m_faces.size() << std::endl;
  for (Triangle tri : m_faces) {
    glm::vec3 BA = m_vertices[tri.v2] - m_vertices[tri.v1];
    glm::vec3 CA = m_vertices[tri.v3] - m_vertices[tri.v1];
    glm::vec3 N = glm::normalize(glm::cross(BA, CA));
    float denom = glm::dot(N, N);

    //std::cout << "Calling Triangle Intersect" << std::endl;
    double t = IntersectTri(orig, dest, 
			    m_vertices[tri.v1], m_vertices[tri.v2], m_vertices[tri.v3]);
    //std::cout << "\n\n\n" << t << std::endl;
    if (t != 0) {
      if ((cur_t == 0) || (t < cur_t)) {
	cur_t = t;
	if (m_vnorms.size() != 0) {
	  glm::vec3 p = orig + cur_t * dest;
	  /*
	  glm::vec3 PA = p - m_vertices[tri.v1];
	  float d00 = glm::dot(BA, BA);
	  float d01 = glm::dot(BA, CA);
	  float d02 = glm::dot(BA, PA);
	  float d11 = glm::dot(CA, CA);
	  float d12 = glm::dot(PA, PA);
	  */
	  //float denom = 1 / (d00 * d11 - d01 * d01);
	  glm::vec3 CB = m_vertices[tri.v3] - m_vertices[tri.v2];
	  glm::vec3 PB = p - m_vertices[tri.v2];
	  float u = glm::dot(N, glm::cross(CB, PB));

	  glm::vec3 AC = m_vertices[tri.v1] - m_vertices[tri.v3];
	  glm::vec3 PC = p - m_vertices[tri.v3];
	  float v = glm::dot(N,glm::cross(AC, PC));

	  u /= denom;
	  v /= denom;

	  // u + v + w = 1                                                              
	  //float v = (d11 * d02 - d01 * d12) * denom;
	  //float w = (d00 * d12 - d01 * d02) * denom;
	  float w = 1 - u - v;

	  n = glm::normalize(u * m_vnorms[tri.v1] + v * m_vnorms[tri.v2] + w * m_vnorms[tri.v3]);
	  //std::cout << m_vnorms[tri.v1][0] << ", " << m_vnorms[tri.v1][1] << ", " << m_vnorms[tri.v1][2] << ", "  << std::endl; 
	  //std::cout << u << ", " << v << ", " << w << std::endl; 
	  //std::cout << n[0] << ", " << n[1] << ", " << n[2] << ", "  << std::endl; 
	} else {

	n = glm::normalize(glm::cross(m_vertices[tri.v2] - m_vertices[tri.v1],
				      m_vertices[tri.v3]- m_vertices[tri.v1]));
	}
      }
    }
  }
  I.t = cur_t;
  I.p = orig + (float)I.t * dest;
  /*
  if (m_vnoms.size() != 0) {
  glm::vec3 PA = I.p - m_vertices[tri.v1];
  float d00 = glm::dot(BA, BA);
  float d01 = glm::dot(BA, CA);
  float d02 = glm::dot(BA, PA);
  float d11 = glm::dot(CA, CA);
  float d12 = glm::ot(PA, PA);
  float denom = 1 / (d00 * d11 - d01 * d01);

  // u + v + w = 1
  float v = (d11 * d02 - d01 * d12) * dnom;
  float w = (d00 * d12 - d01 * d02) * dnom;
  float u = 1 - v - w;
  */

  I.n = n;
  
  return I;
}
