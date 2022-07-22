// Fall 2018
//
#include <glm/ext.hpp>
#include <thread> 
#include <math.h>
#include <vector>
#include <mutex>
#include "A4.hpp"
#include <iostream>
#include <random>

using namespace glm;
using namespace std;

// Epsilon
const double EPSILON = 0.01;

// Tone Shading Coeffients
const double blue = 0.4;
const double yellow = 0.4;
const double alpha = 0.2;
const double beta = 0.6;
const vec3 EarlyWood = vec3(0.725, 0.612, 0.419);
const vec3 OldWood = vec3(0.628, 0.441, 0.237);
const vec3 EarlyWood_1 = vec3(0.505, 0.404, 0.230);
const vec3 OldWood_1 = vec3(0.411, 0.345, 0.267);
const vec3 BabyBlue = vec3(0.404, 0.38, 0.286);
const vec3 sky_b = vec3(0.467, 0.533, 0.612);
const vec3 sky_y = vec3(1, 0.91, 0.73); 
const vec3 sky = vec3(0.733, 0.765, 0.773);

// Rendering Options
long long int COUNT = 0;
long long int WORK = 0;
std::mutex WORK_mutex; 
const int DEPTH = 0;
const bool TONE_SHADE_1 = false;
const bool TONE_SHADE_2 = false;
const bool PHONG_SHADE = true;
const bool ANTI_ALIASING = true;
const bool FEATURE_LINE = false;
const bool FAKE_RENDER = false;
const bool SOFT_SHADOW = false;

std::vector<glm::vec3> SHADOW_RAYS;

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
) {

  // Fill in raytracing code here...  

  std::cout << "Calling A4_Render(\n" <<
		  "\t" << *root <<
          "\t" << "Image(width:" << image.width() << ", height:" << image.height() << ")\n"
          "\t" << "eye:  " << glm::to_string(eye) << std::endl <<
		  "\t" << "view: " << glm::to_string(view) << std::endl <<
		  "\t" << "up:   " << glm::to_string(up) << std::endl <<
		  "\t" << "fovy: " << fovy << std::endl <<
          "\t" << "ambient: " << glm::to_string(ambient) << std::endl <<
		  "\t" << "lights{" << std::endl;

	for(const Light * light : lights) {
		std::cout << "\t\t" <<  *light << std::endl;
	}
	std::cout << "\t}" << std::endl;
	std:: cout <<")" << std::endl;

	size_t H = image.height();
	size_t W = image.width();

	vector<vector <glm::vec3>> buffer(H+2, vector<glm::vec3>(W+2));

	for (int i = 0; i < 5000; ++i) {
	  SHADOW_RAYS.push_back(glm::vec3(glm::ballRand(50.0)));
	}
	
	thread threads[100];
	for (int i = 0; i < 100; ++i) {

	  threads[i] = thread(multi_thread, root, std::ref(image), std::ref(eye), std::ref(view), std::ref(up), fovy, std::ref(ambient), std::ref(lights), H, W, std::ref(buffer));
	 

	}
	
	  
	for (int i = 0; i < 100; ++i) {
	  threads[i].join();
	}
	
       

	/*
	
	// Compute u, v, w basis vectors                                                                    
	glm::vec3 w = glm::normalize(view);
	glm::vec3 u = glm::normalize(glm::cross(w, glm::normalize(up)));
	glm::vec3 v = glm::normalize(glm::cross(u, w));

	int count = 0;
	// for each pixel 
	for (uint y = 0; y < H + 2; ++y) {
		for (uint x = 0; x < W + 2; ++x) {
		  
		  // create a ray
		  glm::vec3 orig = eye; 

		  double d = 1;
		  double aspect = (double)W/(double)H;

		  double virtual_H = 2 * d * glm::tan(glm::radians(fovy / 2));
		  double virtual_W = virtual_H * aspect;

		  glm::mat4 trans_to_center = glm::translate(glm::mat4(), glm::vec3(-(double)W/2, -(double)H/2, d));
		  glm::mat4 scale_p = glm::scale(glm::mat4(), glm::vec3(virtual_W/(double)W, -virtual_H/(double)H, 1));
		  glm::mat4 rotate_p = glm::mat4(1.0);
		  rotate_p[0] = {u[0], v[0], w[0], 0};
		  rotate_p[1] = {u[1], v[1], w[1], 0};
		  rotate_p[2] = {u[2], v[2], w[2], 0};

		  glm::mat4 trans_to_eye = glm::translate(glm::mat4(), eye);
		  glm::vec3 s = glm::vec3(trans_to_eye * rotate_p * scale_p * trans_to_center * glm::vec4(x, y, 0, 1));
		  glm::vec3 dest = glm::normalize(s - eye);

		  Reflection ref = get_color(root, ambient, lights,
                                     orig, dest, 0, -5);
		  
		  if ((0 < x) && (x < W+1) && (0 < y) && (y < H+1)) {
		    image(x-1, y-1, 0) = ref.color[0];
		    image(x-1, y-1, 1) = ref.color[1];
		    image(x-1, y-1, 2) = ref.color[2];
		  }
		  
		  buffer[x][y] = ref.color;
		}

		float percent = ((float)y/(W-1)) * 100;
		std::cout << "\e[A\r\e[0K" << (int)percent << '%' << std::endl;
		
	}
	*/
	
	if (ANTI_ALIASING) {
	  for (int i = 1; i < H+1; ++i) {
	    for (int j = 1; j < W+1; ++ j) {
	      vec3 color = (buffer[i-1][j-1] + buffer[i-1][j] + buffer[i-1][j+1]
			    + buffer[i][j-1] + buffer[i][j] + buffer[i][j+1]
			    + buffer[i+1][j-1] + buffer[i+1][j] + buffer[i+1][j+1]) / 9;
	      image(i-1, j-1, 0) = color[0];
	      image(i-1, j-1, 1) = color[1];
	      image(i-1, j-1, 2) = color[2];
	    }
	  }
	}
        if (FEATURE_LINE) {
	  for (int i = 1; i < H; ++i) {
            for (int j = 1; j < W; ++ j) {
	      vec3 diff_1 = vec3(image(i, j, 0), image(i, j, 1), image(i, j, 2)) - 
		vec3(image(i-1, j, 0), image(i-1, j, 1), image(i-1, j, 2));
	      //cout <<  diff_1 << endl;
	      vec3 diff_2 = vec3(image(i, j, 0), image(i, j, 1), image(i, j, 2)) -
		vec3(image(i, j-1, 0), image(i, j-1, 1), image(i, j-1, 2));
	      //cout <<  diff_2 << endl; 
	      float t = 0.14;
	      float d = -0.1;
	      if ((diff_1[0] > t) || (diff_1[0] < -t) ||
		  (diff_1[1] > t) || (diff_1[1] < -t) || 
		  (diff_1[2] > t) || (diff_1[2] < -t)){
		//cout << 'T' << endl;
		image(i-1, j, 0) += d;
		image(i-1, j, 1) += d;
		image(i-1, j, 2) += d;
	      } 
	      if ((diff_2[0] > t) || (diff_2[0] < -t) ||
		  (diff_2[1] > t) || (diff_2[1] < -t) || 
		  (diff_2[2] > t) || (diff_2[2] < -t)) {
		image(i, j-1, 0) += d;
		image(i, j-1, 1) += d;
		image(i, j-1, 2) += d;
	      }
	    }
	  }
	}
}


void multi_thread(// What to render                                                                                                            
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
		 vector<vector <glm::vec3>> & buffer
		 ) 
{

  while(true) {

    //size_t H = image.height();
    //size_t W = image.width();
    int x;
    int y;

    //std::mutex WORK_mutex;
    {
      std::lock_guard<std::mutex> lock(WORK_mutex); 
      if (WORK < ((H+2) * (W+2))) {
      //std::lock_guard<std::mutex> lock(WORK_mutex);

	y = WORK / (H+2);
	x = WORK % (W+2);
	++WORK;
      //cout << x << " " << y << endl;
        float percent = ((float)WORK/((H+2) * (W+2) - 1)) * 100;
	//std::cout << '/r' << percent << '%' << flush;
	std::cout << "\e[A\r\e[0K" << (int)percent << '%' << std::endl;
      } else {
	return;
      }
    }


    //cout << x << " " << y << endl;
    glm::vec3 w = glm::normalize(view);
    glm::vec3 u = glm::normalize(glm::cross(w, glm::normalize(up)));
    glm::vec3 v = glm::normalize(glm::cross(u, w));

    int count = 0;
    // create a ray                                                                                                            
    glm::vec3 orig = eye;

    double d = 1;
    double aspect = (double)W/(double)H;

    double virtual_H = 2 * d * glm::tan(glm::radians(fovy / 2));
    double virtual_W = virtual_H * aspect;

    glm::mat4 trans_to_center = glm::translate(glm::mat4(), glm::vec3(-(double)W/2, -(double)H/2, d));
    glm::mat4 scale_p = glm::scale(glm::mat4(), glm::vec3(virtual_W/(double)W, -virtual_H/(double)H, 1));
    glm::mat4 rotate_p = glm::mat4(1.0);
    rotate_p[0] = {u[0], v[0], w[0], 0};
    rotate_p[1] = {u[1], v[1], w[1], 0};
    rotate_p[2] = {u[2], v[2], w[2], 0};

    glm::mat4 trans_to_eye = glm::translate(glm::mat4(), eye);
    glm::vec3 s = glm::vec3(trans_to_eye * rotate_p * scale_p * trans_to_center * glm::vec4(x, y, 0, 1));
    glm::vec3 dest = glm::normalize(s - eye);

    
    Reflection ref = get_color(root, ambient, lights,
			       orig, dest, 0, y, H,  -5);
    
 
    if ((0 < x) && (x < W+1) && (0 < y) && (y < H+1)) {
      image(x-1, y-1, 0) = ref.color[0];
      image(x-1, y-1, 1) = ref.color[1];
      image(x-1, y-1, 2) = ref.color[2];
    }
    
    /*
    image(x, y, 0) = ref.color[0];                                                         
    image(x, y, 1) = ref.color[1];                                                         
    image(x, y, 2) = ref.color[2];
    */

    buffer[x][y] = ref.color;
    
  }
}



Reflection get_color(SceneNode * root,
                     const glm::vec3 & ambient,
                     const std::list<Light *> & lights,
                     glm::vec3 orig, glm::vec3 dest,
                     int ref_num,
                     int y, int H, 
int type) {
  glm::vec3 col;
 
  //cout << "intersect" << endl;
    Intersection I = root->Intersect(orig, dest, type);
    
    if (I.t == 0) {
      //col = sky;
      //col = (1- (double)((double)y/H)) * sky + (double)((double)y/H) * sky_y;
      
      col[0] = 0.3;                                                                                                     
      col[1] = 0.4;                                                                                                     
      col[2] = 0.7;  
      
      return Reflection{I, col};
    }

    glm::vec3 kd = I.kd;
    glm::vec3 ks = I.ks;
    double shininess = I.s;
    double kr = I.kr;
    int id  = I.id;
    if (kd[0] == -1) {
      // Wood Solid Texture 1
      float f = fmod((((I.p)[0]+id)*((I.p)[0]+id) + ((I.p)[2]+id)*((I.p)[2])+id), 2);   
      kd = EarlyWood + f * (OldWood - EarlyWood);
    }

    if (kd[0] == -2) {
        //if (I.burning == true) {
            //std::cout << "burning" << std::endl;
            //kd = glm::vec3(0.325, 0.204, 0.047);
        //} else {
      // Stripe solid texture
      int r = (int)((I.p)[1]*50) % 5;
      //std::cout << (I.p)[2]+ 5 << " " << r << std::endl;
            //kd = glm::vec3(1, 1, 1);
    if (r == 0) {kd = glm::vec3(1, 1, 1);}
    else {kd = glm::vec3(0.72, 0.75, 0.81);}
      //  }
      //if ((I.p)[1] <= 0) {kd = glm::vec3(1, 0, 0);
        //}
    }

    if (kd[0] == -3) {
      float f = fmod((((I.p)[0]+id)*((I.p)[0]+id) + ((I.p)[2]+id)*((I.p)[2])+id), 2);
      kd = EarlyWood_1 + f * (OldWood_1 - EarlyWood_1);
    }
    
    if (kd[0] == -4) {
        kd = glm::vec3(1, 1, 1);
        //kd[0] -= (I.temperature - 100) /220;
        //kd[1] -= (I.temperature - 100) /200;
        //kd[2] -= (I.temperature - 100) /200;
        if (I.burning == true) {
            kd = glm::vec3(0.325, 0.204, 0.047);
        }
        
    }
    
    if (kd[0] == -5) {
        kd = glm::vec3(1, 1, 1);
        //kd[0] -= (I.temperature - 100) /220;
        //kd[1] -= (I.temperature/400);
        //kd[2] -= (I.temperature/400);
        
    }
    //glm::vec3 ks = I.ks;
    //kd = vec3(1.0);
    //double shininess = I.s;
    //int id = I.id;    

    //glm::vec3 col;
    vec3 k_blue = vec3(0, 0, blue);
    vec3 k_yellow = vec3(yellow, yellow, 0);
    vec3 k_cool = k_blue + alpha * kd;
    vec3 k_warm = k_yellow + beta * kd;
    if (PHONG_SHADE) {
      col[0] = ambient[0] * kd[0];
      col[1] = ambient[1] * kd[1];
      col[2] = ambient[2] * kd[2];
    }
    else if (TONE_SHADE_1) { 
      //col = (k_cool + k_warm)/2;
    }
    else if (TONE_SHADE_2) {
      col = (k_cool + k_warm)/2;
    }
  
    if (ref_num < DEPTH && kr > 0) {
    glm::vec3 r = glm::normalize(dest - 2 * glm::dot(dest, I.n) *I.n);
    Reflection col_ref = get_color(root, ambient, lights,
                                   I.p, r, ref_num+1,y, H,  I.id);

    if (col_ref.I.t >= 0) {
      //col = 0.5 * col + 0.5 * col_ref.color;
      //col += ks * col_ref.color;
      col = (1-kr) * col + kr * col_ref.color;
      //col = 0.5 * col + 0.5 * col_ref.color; 
    }
  }

  glm::vec3 v = -glm::normalize(dest);
  glm::vec3 sum = glm::vec3(0,0,0);
  
  //int i = 0;
  for (Light * light : lights) {
    //if (i == 0) {continue;}
    //++i;
    glm::vec3 l = light->position - I.p;
    glm::vec3 l2 = -l;
    double q = glm::length(l);
    l = glm::normalize(l);
    l2 = glm::normalize(l2);
    double n_dot_l = glm::dot(I.n, l);
    double n_dot_l2 = glm::dot(I.n, l2);
    glm::vec3 diffuse = kd * std::max(n_dot_l, 0.0);
    glm::vec3 diffuse2 = kd * std::max(n_dot_l2, 0.0);
    float tone_diff = std::max(n_dot_l, 0.0);
    
    bool in_shadow = false;
    Intersection shadow_I = root->Intersect(I.p, l, I.id);
    
    float blocked = 0;
    float casted = 0;
    
    
    if (SOFT_SHADOW) {
    for (int i = 0; i < 100; ++i) {
	  ++casted;
	  //static thread_local std::mt19937 generator;
	  //std::uniform_int_distribution<int> distribution(0, 4999);
	  //glm::vec3 l_disk = light->position+SHADOW_RAYS[distribution(generator)] - I.p;
	  glm::vec3 l_disk = light->position+SHADOW_RAYS[COUNT % 5000] - I.p;
	  Intersection soft_shadow_I = root->Intersect(I.p, l_disk, I.id);
	  ++COUNT;
	  if (soft_shadow_I.t != 0) {
	    ++blocked;
	  }
    }
    }
   
    
    else {
      if (shadow_I.t != 0) {
	in_shadow = true;
      }
      
      if (in_shadow && PHONG_SHADE) {
	continue;
      }
      
    }
    

    glm::vec3 r = glm::normalize(2 * (glm::dot(l, I.n)) * I.n - l);
    glm::vec3 r2 = glm::normalize(2 * (glm::dot(l2, I.n)) * I.n - l2);
    glm::vec3 h = glm::normalize(v + l);
    glm::vec3 h2 = glm::normalize(v + l2);
    double r_dot_v = glm::dot(r, v);
    double r_dot_v2 = glm::dot(r2, v);
    double n_dot_h = glm::dot(I.n, h);
    double n_dot_h2 = glm::dot(I.n, h2);
    glm::vec3 specular = ks * std::pow(std::max(r_dot_v, 0.0), shininess);
    glm::vec3 specular2 = ks * std::pow(std::max(r_dot_v2, 0.0), shininess);
    double attenuation = 1 / (light->falloff[0] + light->falloff[1]*q + light->falloff[2]*q*q);

    vec3 tone_shade = mix(k_cool, k_warm, tone_diff) + specular;

    if (TONE_SHADE_1) {
      col += tone_shade * attenuation;
      if (SOFT_SHADOW) { {col /= (1+blocked/casted);} }
      else { 
	if (in_shadow) {col /= 2;}
      }
    }

    else if (TONE_SHADE_2) {
      col += (k_warm - k_cool)/2 * (diffuse + specular);
      col += (k_cool - k_warm)/2 * (diffuse2 + specular2);
    }

    else if (PHONG_SHADE) {
      if (SOFT_SHADOW) {
        col += (light->colour * (diffuse + specular)) * attenuation * (1-blocked/casted);
      } else {
	col += (light->colour * (diffuse + specular)) * attenuation;
      }
      //if (SOFT_SHADOW) { col *= (1-(float)(blocked/casted)); }// * (1-blocked/casted);
    }
    //++i;
  }
  
  if (FAKE_RENDER) {
    int r = (id & 0x000000FF) >> 0;
    int g = (id & 0x0000FF00) >> 8;
    int b = (id & 0x00FF0000) >> 16;
    col = vec3(r/255, g/255, b/255);
  }
  else if (TONE_SHADE_1) {
    col /= lights.size(); // + 1;  
  }
  
  else if (TONE_SHADE_2) {
    //col += (k_cool + k_warm)/2;
  }
  /*
  if (I.t == 0) {
    col[0] = 0.3;
    col[1] = 0.4;
    col[2] = 0.7;
    }
  */
   return Reflection{I, col};
}


