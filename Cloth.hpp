// Fall 2018                                                                                               

#pragma once

#include <vector>
#include <iosfwd>
#include <string>

#include <glm/glm.hpp>

#include "Mesh.hpp"

struct Constraint {
  int v1;
  int v2;
    float rest;
  float initial;
    char type;

  Constraint(int pv1, int pv2, float len, char c)
    : v1( pv1 )
    , v2( pv2 )
    , rest( len )
    , initial( len )
    , type( c )
  {}

};

class Cloth : public Mesh {
 public:
  Cloth(int N, std::string geometry, glm::vec3 center, float r, float softness, int time, int start);

  void add_force(glm::vec3 force, int index);

  void apply_force(int index);
  
  void move_vertex(int index);

  void check_constraints();
    
  void update_constraints();

  void pin(int index);
 
  void collision (std::string geometry, glm::vec3 center, float r, int index);
    
  //void burn (int index);
    
  void heat_transfer(int index1, int index2);
    
  void burn_cloth(int N);
    
  Intersection Intersect(glm::vec3 orig, glm::vec3 dest) override;

 private:
  std::vector<bool> m_burntfaces;
  std::vector<bool> m_burningfaces;

  std::vector<bool> m_burntvertices;
  std::vector<bool> m_burningvertices;
  std::vector<float> m_temperature;
  std::vector<float> m_lasttemp;
  std::vector<int> m_timer;
  std::vector<bool> m_movable;
  std::vector<glm::vec3> m_last;
  std::vector<float> m_masses;
  std::vector<glm::vec3> m_forces;
  std::vector<glm::vec3> m_acc;
    
  std::vector<bool> m_burntsprings;
  std::vector<Constraint> m_springs;
  //std::vector<glm::vec3> m_normals; //Phong Shade
};
