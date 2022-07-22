// Fall 2018

#pragma once

#include "SceneNode.hpp"
//#include "Primitive.hpp"
#include "Material.hpp"

class GeometryNode : public SceneNode {
public:
	GeometryNode( const std::string & name, Primitive *prim, 
		Material *mat = nullptr );

	void setMaterial( Material *material );

  Material *get_material() const override {return m_material;};
  Primitive *get_primitive() const override {return m_primitive;};

  Material *m_material;
	Primitive *m_primitive;

  Intersection Intersect(glm::vec3 orig, glm::vec3 dest, int type) const override;
};
