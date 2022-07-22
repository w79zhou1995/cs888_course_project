// Fall 2018

#include "SceneNode.hpp"

#include "cs488-framework/MathUtils.hpp"

#include <iostream>
#include <sstream>
using namespace std;

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtx/transform.hpp>

using namespace glm;
//double e = 0.1;

// Static class variable
unsigned int SceneNode::nodeInstanceCount = 0;


//---------------------------------------------------------------------------------------
SceneNode::SceneNode(const std::string& name)
  : m_name(name),
	m_nodeType(NodeType::SceneNode),
	trans(mat4()),
	invtrans(mat4()),
	m_nodeId(nodeInstanceCount++)
{

}

//---------------------------------------------------------------------------------------
// Deep copy
SceneNode::SceneNode(const SceneNode & other)
	: m_nodeType(other.m_nodeType),
	  m_name(other.m_name),
	  trans(other.trans),
	  invtrans(other.invtrans)
{
	for(SceneNode * child : other.children) {
		this->children.push_front(new SceneNode(*child));
	}
}

//---------------------------------------------------------------------------------------
SceneNode::~SceneNode() {
	for(SceneNode * child : children) {
	  delete child;
	}
}

//---------------------------------------------------------------------------------------
void SceneNode::set_transform(const glm::mat4& m) {
	trans = m;
	invtrans = glm::inverse(m);
}

//---------------------------------------------------------------------------------------
const glm::mat4& SceneNode::get_transform() const {
	return trans;
}

//---------------------------------------------------------------------------------------
const glm::mat4& SceneNode::get_inverse() const {
	return invtrans;
}

//---------------------------------------------------------------------------------------
void SceneNode::add_child(SceneNode* child) {
	children.push_back(child);
}

//---------------------------------------------------------------------------------------
void SceneNode::remove_child(SceneNode* child) {
	children.remove(child);
}

//---------------------------------------------------------------------------------------
void SceneNode::rotate(char axis, float angle) {
	vec3 rot_axis;

	switch (axis) {
		case 'x':
			rot_axis = vec3(1,0,0);
			break;
		case 'y':
			rot_axis = vec3(0,1,0);
	        break;
		case 'z':
			rot_axis = vec3(0,0,1);
	        break;
		default:
			break;
	}
	mat4 rot_matrix = glm::rotate(degreesToRadians(angle), rot_axis);
	set_transform( rot_matrix * trans );
}

//---------------------------------------------------------------------------------------
void SceneNode::scale(const glm::vec3 & amount) {
	set_transform( glm::scale(amount) * trans );
}

//---------------------------------------------------------------------------------------
void SceneNode::translate(const glm::vec3& amount) {
	set_transform( glm::translate(amount) * trans );
}


//---------------------------------------------------------------------------------------
int SceneNode::totalSceneNodes() const {
	return nodeInstanceCount;
}

//---------------------------------------------------------------------------------------
std::ostream & operator << (std::ostream & os, const SceneNode & node) {

	//os << "SceneNode:[NodeType: ___, name: ____, id: ____, isSelected: ____, transform: ____"
	switch (node.m_nodeType) {
		case NodeType::SceneNode:
			os << "SceneNode";
			break;
		case NodeType::GeometryNode:
			os << "GeometryNode";
			break;
		case NodeType::JointNode:
			os << "JointNode";
			break;
	}
	os << ":[";

	os << "name:" << node.m_name << ", ";
	os << "id:" << node.m_nodeId;

	os << "]\n";
	return os;
}

Intersection SceneNode::Intersect(glm::vec3 orig, glm::vec3 dest, int type) const{
  glm::vec3 torig = glm::vec3(invtrans * glm::vec4(orig, 1));
  glm::vec3 tdest = glm::normalize(glm::vec3(invtrans * glm::vec4(dest, 0)));

  //cur_I.id = m_nodeId;    
  Intersection cur_I;  
  // cur_I.id = m_nodeId; 
  //double e = 0.01;


  //Intersection cur_I;
  //int count = 0;
  for (const SceneNode * node : this->children) {
    if (type == node->m_nodeId) {
      //e = 1;
      // continue;
     } 
    /*
    if (node->m_nodeType == NodeType::GeometryNode) {
      //std::cout << count << std::endl;
      //++count;
      Intersection I = (node->get_primitive())->Intersect(orig, dest);
      if (I.t > std::abs(e)) {
	if ((cur_I.t < std::abs(e)) || (I.t < cur_I.t)) {
	  cur_I = I;
	  cur_I.kd = node->get_material()->get_kd();
	  cur_I.ks = node->get_material()->get_ks();
	  cur_I.s = node->get_material()->get_s();
	}
      }
      }*/
    Intersection I ;
    //I.id = node->m_nodeId;  
    I = node->Intersect(torig, tdest, type);
    //I.id = node->m_nodeId;
    /*    
    if (I.id == m_nodeId) {
      e = 1;
    }
    */
    //if (I.t > std::abs(e)){
    if (I.t > 0) {                                                                                            
      //if ((cur_I.t < std::abs(e)) || (glm::distance(torig, I.p) < glm::distance(torig, cur_I.p))) {
	if ((cur_I.t == 0) || (glm::distance(torig, I.p) < glm::distance(torig, cur_I.p))) {
	cur_I = I;
	//cur_I.id = node->m_nodeId;
	}
    }
  }
  
  
  if (cur_I.t != 0) {
    cur_I.p = glm::vec3(trans * glm::vec4(cur_I.p, 1));
    cur_I.n = glm::normalize(glm::transpose(glm::mat3(invtrans)) * cur_I.n);
  }
  
  return cur_I;

}
