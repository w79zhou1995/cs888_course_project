// Fall 2018

#pragma once

#include <glm/glm.hpp>

#include "Material.hpp"

class PhongMaterial : public Material {
public:
  PhongMaterial(const glm::vec3& kd, const glm::vec3& ks, double shininess, double kr);
  virtual ~PhongMaterial();
  glm::vec3 get_kd() override {return m_kd;};  
  glm::vec3 get_ks() override {return m_ks;};
  double get_s() override {return m_shininess;};
  double get_rf() override {return m_kr;};

private:
  glm::vec3 m_kd;
  glm::vec3 m_ks;
  double m_kr;
  double m_shininess;
};
