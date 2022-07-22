// Fall 2018                                                                                                    

#pragma once

#include <glm/glm.hpp>
#inlcude "SceneNode.hpp"

struct Intersection
{
  double t = 0;
  glm::vec3 p;
  glm::vec3 n;
  GeometryNode *gn;
  bool burning = false;

};
