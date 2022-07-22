// Fall 2018

#pragma once

#include <glm/glm.hpp>

class Material {
public:
  virtual ~Material();
  virtual glm::vec3 get_kd() {return glm::vec3(0,0,0);};
  virtual glm::vec3 get_ks() {return glm::vec3(0,0,0);};
  virtual double get_s() {return 0;};
  virtual double get_rf() {return 0;};

protected:
  Material();
};
