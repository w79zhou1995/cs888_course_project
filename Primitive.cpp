// Fall 2018

#include "Primitive.hpp"
#include "polyroots.hpp"
#include <math.h>

//const double epsilon = 0.00000;

Primitive::~Primitive()
{
}

Sphere::~Sphere()
{
}

Cube::~Cube()
{
}

Cylinder::~Cylinder()
{
}

Cone::~Cone()
{
}

NonhierSphere::~NonhierSphere()
{
}

NonhierBox::~NonhierBox()
{
}


Intersection Primitive::Intersect(glm::vec3 orig, glm::vec3 dest) {
  Intersection I;
  return I;
}


// Moller and Trumbores' algorithm 
double Primitive::IntersectTri(glm::vec3 orig, glm::vec3 dest, 
			       glm::vec3 V0, glm::vec3 V1, glm::vec3 V2) {

  // calculate two edges from V0
  glm::vec3 E1 = V1 - V0;
  glm::vec3 E2 = V2 - V0;

  // calculate determinant
  glm::vec3 P = glm::cross(dest, E2);
  double det = glm::dot(E1, P);

  // if determinant is 0, ray and triangle are in the same plane, no intersection
  if ((det < epsilon) && (det > -epsilon)) {
    return 0;
  }

  // T(u, v) = (1-u-v)V0 + uV1 + vV2   
  //   where u,v >= 0, u+v <= 1

  glm::vec3 T = orig - V0;
  double u = (1 / det) * glm::dot(T, P);
  // out of bound
  if ((u < 0) || (u > 1)) {
    return 0;
  }

  glm::vec3 Q = glm::cross(T, E1);
  double v = (1 / det) * glm::dot(dest, Q);
  // out of bound                                                                                    
  if ((v < 0) || (u + v > 1)) {
    return 0;
  }

  double t = (1 / det) * glm::dot(E2, Q);
  if (t > epsilon) {
    return t;
  } 
  return 0;
}


Intersection Cube::Intersect(glm::vec3 orig, glm::vec3 dest) {
  glm::vec3 m_pos = glm::vec3(0,0,0);
  double m_size = 1;
  //std::vector<glm::vec3> m_vertices;
  //std::vector<Triangle> m_faces;

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

  Intersection I;
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
}


Intersection Sphere::Intersect(glm::vec3 orig, glm::vec3 dest) {
  Intersection I;
  glm::vec3 c = glm::vec3(0,0,0);
  double r = 1;
  double A = glm::dot(dest, dest);
  double B = 2 * glm::dot(dest, (orig - c));
  double C = glm::dot((orig - c), (orig - c)) - r * r;
  double roots[2];
  size_t n = quadraticRoots(A, B, C, roots);

  if (n == 0) {
    return I;
  }

  else if (n == 1) {
    if (roots[0] > epsilon) {
      I.t = roots[0];
      glm::vec3 p = orig + (float)I.t * dest;                                                         

      I.p = p;
      I.n = glm::normalize(p - c);
      return I;
    }
    return I;
  }

  else if (n == 2) {
    if ((roots[0] <  roots[1]) && (roots[0] > epsilon)) {
      I.t = roots[0];
      glm::vec3 p = orig + (float)I.t * dest;
      I.p = p;
      I.n = glm::normalize(p - c);
      return I;
    }
    else if (roots[1] > epsilon) {
      I.t = roots[1];
      glm::vec3 p = orig + (float)I.t * dest;
      I.p = p;
      I.n = glm::normalize(p - c);
      return I;
    }
    else {
      return I;
    }
  }
  return I;
}

Intersection Cylinder::Intersect(glm::vec3 orig, glm::vec3 dest) {
  Intersection I;
  glm::vec3 c = glm::vec3(0,0,0);
  double r = 1;
  double A = dest[0]*dest[0] + dest[1]*dest[1];
  double B = 2*(orig - c)[0]*dest[0] + 2*(orig - c)[1]*dest[1];
  double C = (orig - c)[0]*(orig - c)[0] + (orig - c)[1]*(orig - c)[1] - r*r;
  double roots[2];
  size_t n = quadraticRoots(A, B, C, roots);

  if (n == 0) {
    return I;
  }

  else if (n == 1) {
    
    double z = (orig - c)[2] + roots[0]*dest[2];
    
    if (roots[0] > epsilon) {
      if ((0 <= z) && (z <= 1)) {
	I.t = roots[0];
	glm::vec3 p = orig + (float)I.t * dest;	

	I.p = p;
	I.n = glm::normalize(p - glm::vec3(0,0,I.p[2]));
	return I;
      }
    }
    
    return I;
  }
  else if (n == 2) {
    double z1 = (orig - c)[2] +roots[0]*dest[2];
    double z2 = (orig - c)[2] +roots[1]*dest[2];
    double t1 = roots[0];
    double t2 = roots[1];
    if (t1 > t2) {
      std::swap(t1, t2);
      std::swap(z1, z2);
    }
   
    // choose root[0]
    //if ((roots[0] <  roots[1]) && (roots[0] > epsilon)) {
    
      if (z1 < 0) {
	if (z2 < 0) {
	  return I;
	} else {
	  double t = (double)((0 - (orig - c)[2])/dest[2]);                                      
	  // if (t > 0) {
	    I.t = t;   
	  glm::vec3 p = orig + (float)I.t * dest;
	  I.p = p;                                                                                        
	  I.n = glm::vec3(0,0,-1);
	  return I; 
	  // } else {
	  //  return I;
	  // }
	}
      }
      
      else if ((0 <= z1) && (z1 < 1)) {
	I.t = t1; //roots[0];
	glm::vec3 p = orig + (float)I.t * dest;
	I.p = p;
	I.n = glm::normalize(p - glm::vec3(0,0,I.p[2]));
	//I.n = glm::normalize(glm::vec3(I.p[0], I.p[1], 0)); 
	return I;
      } 
      else {
	if (z2 > 1) {
	  return I;
	} else {
	  double t = (double)((1 - (orig - c)[2])/dest[2]);
	  //if (t > 0) {
	    I.t = t;
	  glm::vec3 p = orig + (float)I.t * dest;
	  I.p = p;
	  I.n = glm::vec3(0,0,1);
          return I;
	  // } else {
	  // return I;
	  //}
	}
      }
  }
      
  

  /* choose root[1]
    else if (roots[1] > epsilon) {
    
      if (z2 < 0) {
	if (z1 < 0) {
	  return I;
	}
	
	else {
	  double t = (double)((0 - (orig - c)[2])/dest[2]);                                                             
	  // if (t > 0) {
	    I.t = t;
	  glm::vec3 p = orig + (float)I.t * dest;
	  I.p = p;                                                                                                
	  I.n = glm::vec3(0,0,-1);
	  return I;
	  // } else {
	  //  return I;
	  // }
	}
      }
      if ((0 <= z2) && (z2 < 1)) {
	I.t = roots[1];
	glm::vec3 p = orig + (float)I.t * dest;
	I.p = p;
	I.n = glm::normalize(p - glm::vec3(0,0,I.p[2]));
	//I.n = glm::normalize(glm::vec3(I.p[0], I.p[1], 0));
	return I;
      }
      else {
	if (z1 > 1) {
	  return I;
	} 
	
	else {
          double t = (double)((1 - (orig - c)[2])/dest[2]);
	  //double t = roots[0] + (roots[1] - roots[0]) * (z1 - 1) / (z1 - z2);
	  // if (t > 0) {
	  I.t = t;
	  glm::vec3 p = orig + (float)I.t * dest;
          I.p = p;
          I.n = glm::normalize(glm::vec3(0,0,1));
          return I;
	  // } else {
	  // return I;
	  //}
	  //}
	
      
	}      
      }
  }
    */
    else {
      return I;
    }

  return I;
}






Intersection Cone::Intersect(glm::vec3 orig, glm::vec3 dest) {
  Intersection I;
  glm::vec3 c = glm::vec3(0,0,0);
  double r = 1;
  double A = dest[0]*dest[0] + dest[1]*dest[1] - dest[2]*dest[2];
  double B = 2*(orig - c)[0]*dest[0] + 2*(orig - c)[1]*dest[1] - 2*(orig - c)[2]*dest[2];
  double C = (orig - c)[0]*(orig - c)[0] + (orig - c)[1]*(orig - c)[1] - (orig - c)[2]*(orig - c)[2];
  double roots[2];
  size_t n = quadraticRoots(A, B, C, roots);

  if (n == 0) {
    return I;
  }

  else if (n == 1) {
    
    double z = (orig - c)[2] + roots[0]*dest[2];

    if (roots[0] > epsilon) {
      if ((0 <= z) && (z <= 1)) {
        I.t = roots[0];
	glm::vec3 p = orig + (float)I.t * dest;

        I.p = p;
        I.n = glm::vec3(0,0,1);
        return I;
      }
    }
    return I;
  }
  else if (n == 2) {
    double z1 = (orig - c)[2] +roots[0]*dest[2];
    double z2 = (orig - c)[2] +roots[1]*dest[2];
    double t1 = roots[0];
    double t2 = roots[1];
    // choose root[0]    

    if (t1 > t2) {
      std::swap(t1, t2);
      std::swap(z1, z2);
    }                                                                                                                                                                                                               
    //if ((roots[0] <  roots[1]) && (roots[0] > epsilon)) {
      if (z1 < 0) {
        if (z2 < 0) {
          return I;
        } else {
	  /*
	  if (t2 > epsilon) {
	    I.t = t2;
	    glm::vec3 p = orig + (float)I.t * dest;
	    I.p = p;
	    glm::vec3 v = glm::vec3((I.p - c)[0], (I.p - c)[1], 0);
	    //double m = sqrt(I.p[0]*I.p[0] + I.p[1]*I.p[1]);                                               
	    double m = sqrt(v[0]*v[0] + v[1]*v[1]);
	    v[0] /= m;
	    v[1] /= m;
	    I.n = glm::normalize(glm::vec3(v[0]/r, v[1]/r, r));
	    return I;
	  } else {
	    return I;
	  } */
	  
          double t = (double)((0 - (orig - c)[2])/dest[2]);
	  glm::vec3 p = orig + (float)I.t * dest;

	  if (p[0]*p[0] + p[1]*p[1] <= 1) {
	    I.t = t;
          I.p = p;
          I.n = glm::vec3(0,0,-1);
          return I;
	  } else {
	  return I;
	  }
	  
        }
      }
      else if ((0 <= z1) && (z1 < 1)) {
        I.t = t1; //roots[0];
	glm::vec3 p = orig + (float)I.t * dest;
        I.p = p;
	glm::vec3 v = glm::vec3((I.p - c)[0], (I.p - c)[1], 0);
	double m = sqrt(I.p[0]*I.p[0] + I.p[1]*I.p[1]);
	//double m = sqrt(v[0]*v[0] + v[1]*v[1]);
	v[0] /= m;
	v[1] /= m;
	//I.n = glm::normalize(glm::vec3(v[0]/r, v[1]/r, r));
        I.n = glm::normalize(glm::cross(glm::vec3(0, 1, I.p[1]/m),
				glm::vec3(1, 0, I.p[0]/m)));
        return I;
      }
      else {
	//return I;


        if (z2 > 1) {
          return I;
        } else {
	  double t = (double)((1 - (orig - c)[2])/dest[2]);
	  glm::vec3 p = orig + (float)t * dest;

          if (p[0]*p[0] + p[1]*p[1] <= 1) {
            I.t = t;
            I.p = p;
            I.n = glm::vec3(0,0,1);
            return I;
	    } else {
            return I;
	    }
	    
	  /*
          I.t = (double)((1 - (orig - c)[2])/dest[2]);
	  glm::vec3 p = orig + (float)I.t * dest;
          I.p = p;
          I.n = glm::vec3(0,0,1);
          return I;
	  */
      }
      
      //}

      /* choose root[1]                                                                      
    else if (roots[1] > epsilon) {
      
      if (z2 < 0) {
        if (z1 < 0) {
          return I;
        }
        else {
          I.t = (double)((0 - (orig - c)[2])/dest[2]);
	  glm::vec3 p = orig + (float)I.t * dest;
          I.p = p;
          I.n = glm::vec3(0,0,-1);
          return I;
        }
      }
      
      if ((0 <= z2) && (z2 < 1)) {
        I.t = roots[1];
	glm::vec3 p = orig + (float)I.t * dest;
        I.p = p;
        double m = sqrt(I.p[0]*I.p[0] + I.p[1]*I.p[1]);
        I.n = glm::normalize(glm::cross(glm::vec3(0, 1, I.p[1]/m),
                                        glm::vec3(1, 0, I.p[0]/m)));
        return I;
      }
      
      else {
        if (z1 > 1) {
          return I;
        } 
	
	else {
          double t = (double)((1 - (orig - c)[2])/dest[2]);
	  glm::vec3 p = orig + (float)t * dest;
 
	  if (p[0]*p[0] + p[1]*p[1] <= 1) {
	    I.t = t;
	    I.p = p;
	    I.n = glm::vec3(0,0,-1);
	    return I;
	  } else {
	    return I;
	  }
	}
      }
      */
      }
  }
  else {
      return I;
  }
  
  return I;
}

Intersection NonhierSphere::Intersect(glm::vec3 orig, glm::vec3 dest) {
  Intersection I;
  double A = glm::dot(dest, dest);
  double B = 2 * glm::dot(dest, (orig - m_pos));
  double C = glm::dot((orig - m_pos), (orig - m_pos)) - m_radius * m_radius;
  double roots[2];
  size_t n = quadraticRoots(A, B, C, roots);

  if (n == 0) {
    return I;
  }
  else if (n == 1) {
    if (roots[0] > epsilon) {
      I.t = roots[0];
      glm::vec3 p = orig + (float)I.t * dest;                                                                            
      I.p = p;
      I.n = glm::normalize(p - m_pos);
      return I;
    }
    return I;
  }
  else if (n == 2) {
    if ((roots[0] <  roots[1]) && (roots[0] > epsilon)) {
      I.t = roots[0];
      glm::vec3 p = orig + (float)I.t * dest;
      I.p = p;
      I.n = glm::normalize(p - m_pos);
      return I;
    } 
    else if (roots[1] > epsilon) {
      I.t = roots[1];
      glm::vec3 p = orig + (float)I.t * dest;
      I.p = p;
      I.n = glm::normalize(p - m_pos);
      return I;
    }
    else {
      return I;
    }
  }
  return I;
}

Intersection NonhierBox::Intersect(glm::vec3 orig, glm::vec3 dest) {
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

  Intersection I;
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
}

