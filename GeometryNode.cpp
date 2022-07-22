// Fall 2018

#include "GeometryNode.hpp"

//---------------------------------------------------------------------------------------
GeometryNode::GeometryNode(
	const std::string & name, Primitive *prim, Material *mat )
	: SceneNode( name )
	, m_material( mat )
	, m_primitive( prim )
{
	m_nodeType = NodeType::GeometryNode;
}

void GeometryNode::setMaterial( Material *mat )
{
	// Obviously, there's a potential memory leak here.  A good solution
	// would be to use some kind of reference counting, as in the 
	// C++ shared_ptr.  But I'm going to punt on that problem here.
	// Why?  Two reasons:
	// (a) In practice we expect the scene to be constructed exactly
	//     once.  There's no reason to believe that materials will be
	//     repeatedly overwritten in a GeometryNode.
	// (b) A ray tracer is a program in which you compute once, and 
	//     throw away all your data.  A memory leak won't build up and
	//     crash the program.

	m_material = mat;
}


Intersection GeometryNode::Intersect(glm::vec3 orig, glm::vec3 dest, int type) const{
  glm::vec3 torig = glm::vec3(invtrans * glm::vec4(orig, 1));
  glm::vec3 tdest = glm::normalize(glm::vec3(invtrans * glm::vec4(dest, 0)));

  //double e = 0.01;
  //if (type == m_nodeId) {m_primitive->setEp(1);}
  
  Intersection cur_I;
  if (type == m_nodeId) {
    return cur_I;
    m_primitive->setEp(1);
  } else {
    m_primitive->setEp(0.00001);
  }
cur_I = m_primitive->Intersect(torig, tdest);
    
//cur_I.id = m_nodeId;
  //cur_I.id = m_nodeId;
  if (cur_I.t > 0) {    
    cur_I.kd = m_material->get_kd();
    cur_I.ks = m_material->get_ks();
    cur_I.s = m_material->get_s();
    cur_I.kr = m_material->get_rf();
    cur_I.id = m_nodeId;
  }


  for (const SceneNode * node : this->children) {
    //        if (type ==node->m_nodeId){
    //       continue;

    // }
    //if (node->m_nodeType == NodeType::GeometryNode) {                                                      
	    Intersection I;   
	    //I.id = node->m_nodeId;      
            I = node->Intersect(torig, tdest, type);
    //I.id = node->m_nodeId;
    /*    
    if (I.id == m_nodeId) {
      e = 1;
      }*/
    
    //if (I.t > std::abs(e)){
      if (I.t > 0) {
	//if ((cur_I.t < std::abs(e)) || (glm::distance(torig, I.p) < glm::distance(torig, cur_I.p))) {
      if ((cur_I.t == 0) || (glm::distance(torig, I.p) < glm::distance(torig, cur_I.p))) {
	cur_I = I;
	//cur_I.kd = node->get_material()->get_kd();
	//cur_I.ks = node->get_material()->get_ks();
	//cur_I.s = node->get_material()->get_s();
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

