// Fall 2018

#pragma once

#include <glm/glm.hpp>

#include "SceneNode.hpp"
#include "GeometryNode.hpp"
#include "Light.hpp"
#include "Image.hpp"

struct Reflection
{
  Intersection I;
  glm::vec3 color;
};

void A4_Render(
		// What to render
		SceneNode * root,

		// Image to write to, set to a given width and height
		Image & image,

		// Viewing parameters
		const glm::vec3 & eye,
		const glm::vec3 & view,
		const glm::vec3 & up,
		double fovy,

		// Lighting parameters
		const glm::vec3 & ambient,
		const std::list<Light *> & lights
);

void multi_thread(  // What to render                                                                                                      
		  SceneNode * root,

		  // Image to write to, set to a given width and height                                                                  
		  Image & image,

		  // Viewing parameters                                                                                                  
		  const glm::vec3 & eye,
		  const glm::vec3 & view,
		  const glm::vec3 & up,
		  double fovy,

		  // Lighting parameters                                                                                                 
		  const glm::vec3 & ambient,
		  const std::list<Light *> & lights,
		  size_t H,
		  size_t W,
		  std::vector<std::vector <glm::vec3>> & buffer
		    );


Reflection get_color(SceneNode * root,

                    // Image to write to, set to a given width and height                                                \
                                                                                                                          

                    // Viewing parameters                                                                                 

                    // Lighting parameters                                                                               \
                                                                                                                          
                    const glm::vec3 & ambient,
                    const std::list<Light *> & lights,
                    glm::vec3 orig, glm::vec3 dest,
                    int ref_num,
		    int y, int H, 
int type);
