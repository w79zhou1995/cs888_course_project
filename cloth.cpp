// Fall 2018                                                                                               

#include <iostream>
#include <fstream>

#include <glm/ext.hpp>
#include <map>
#include <string>
#include <random>
#include <stdlib.h>
#include <time.h>
// #include "cs488-framework/ObjFileDecoder.hpp"                                                           
#include "Cloth.hpp"
float TIME_STEP = 0.01;
int TIME = 2000;
int BURNING_START = 1700;
float k = 2;
float damping = 0.01;
float diffusion = 13.2;
float pyrolysis = 150;
float burningtime = 60;
//glm::vec3 GRAVITY = glm::vec3(0, -5, 0);

Cloth::Cloth (int N, std::string geometry, glm::vec3 center, float r, float softness, int time, int start) {
  m_max = glm::vec3(N, 0, N);
  m_min = glm::vec3(0, 0, 0);
  glm::vec3 GRAVITY = glm::vec3(0, -softness, 0);
  

  /* Init Cloth */
  // add vertices  
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      glm::vec3 V = glm::vec3(j, 0, i);
      m_vertices.push_back(V);
      m_last.push_back(V);
      m_vnorms.push_back(glm::vec3(0));
      m_facecount.push_back(0);
      m_masses.push_back(1);
      m_forces.push_back(glm::vec3(0));
      m_acc.push_back(glm::vec3(0));
      m_movable.push_back(true);
      m_temperature.push_back(100);
      m_lasttemp.push_back(100);
      m_timer.push_back(0);
      m_burntvertices.push_back(false);
      m_burningvertices.push_back(false);
    }
  }

  // add faces
  for (int i = 0; i < N-1; ++i) {
    for (int j = 0; j < N-1; ++j) {
        m_faces.push_back( Triangle(j+i*N, j+i*N+N, j+i*N+1) );
        m_burningfaces.push_back(false);
        m_burntfaces.push_back(false);
        m_faces.push_back( Triangle(j+i*N+N, j+i*N+N+1, j+i*N+1) );
        m_burningfaces.push_back(false);
        m_burntfaces.push_back(false);
    }
  }

  // add structural constraints
  for (int i = 0; i < N-1; ++i) {
    for (int j = 0; j < N-1; ++j) {
        m_springs.push_back( Constraint(j+i*N, j+i*N+1, glm::length(m_vertices[j+i*N] - m_vertices[j+i*N+1]), 'c') );
        m_burntsprings.push_back(false);
        m_springs.push_back( Constraint(j+i*N, j+i*N+N, glm::length(m_vertices[j+i*N] - m_vertices[j+i*N+N]), 'c') );
        m_burntsprings.push_back(false);
    }
  }
  
  // add shear constraints
  for (int i = 0; i < N-1; ++i) {
    for (int j = 0; j < N-1; ++j) {
        m_springs.push_back( Constraint(j+i*N, j+i*N+N+1, glm::length(m_vertices[j+i*N] - m_vertices[j+i*N+N+1]), 's') );
        m_burntsprings.push_back(false);
        m_springs.push_back( Constraint(j+i*N+N, j+i*N+1, glm::length(m_vertices[j+i*N+N] - m_vertices[j+i*N+1]), 's') );
        m_burntsprings.push_back(false);
    }
  }
  
  for (int i = 0; i < N-2; ++i) {
    for (int j = 0; j < N-2; ++j) {
        m_springs.push_back( Constraint(j+i*N, j+i*N+2, glm::length(m_vertices[j+i*N] - m_vertices[j+i*N+2]), 'b') );
        m_burntsprings.push_back(false);
      //std::cout << j+i*N << " " << j+i*N+2 <<std::endl;
        m_springs.push_back( Constraint(j+i*N, j+i*N+2*N, glm::length(m_vertices[j+i*N] - m_vertices[j+i*N+2*N]), 'b') );
        m_burntsprings.push_back(false);
      //std::cout << j+i*N << " " << j+i*N+2*N << std::endl;
        m_springs.push_back( Constraint(j+i*N, j+i*N+2*N+2, glm::length(m_vertices[j+i*N] - m_vertices[j+i*N+2*N+2]), 'b') );
        m_burntsprings.push_back(false);
      //std::cout << j+i*N << " " << j+i*N+2*N+2 << std::endl;
        m_springs.push_back( Constraint(j+i*N+2*N, j+i*N+2, glm::length(m_vertices[j+i*N+2*N] - m_vertices[j+i*N+2]), 'b') );
        m_burntsprings.push_back(false);
      //std::cout << j+i*N+2*N << " " << j+i*N+2 << std::endl;
    }
  }
    
  
    if (geometry == "cube") {
        for (int i = 0; i < N * N; ++i) {
            if (m_vertices[i][0] > center[0]-r && m_vertices[i][0] < center[0]+r && m_vertices[i][2] > center[2]-r && m_vertices[i][2] < center[2]+r ) {
                pin(i);
            }
        }
    }
   
    
    if (geometry == "hang" || geometry == "ball") {
        pin(0);
        //pin(N);
        //pin(N-2);
        //pin(N-1);
    }
    
    if (geometry == "cap") {
        pin(0);
        pin(N-1);
    }
    
    
  /* Set Ignition Point */
  m_temperature[N*N] = 20000;

  /* Time Integration */
  for (int t = 0; t < time; ++t) {
      
      /* Moving Time Step */
      for (int i = 0; i < N * N; ++i) {
          // apply gravity
          add_force(GRAVITY * m_masses[i], i);

          apply_force(i);

      }
      
      //update_constraints();
      
      check_constraints();
      

      for (int i = 0; i < N * N; ++i) {
        move_vertex(i);
      }
      
      /* Solid Collision Detection */
      for (int i = 0; i < N * N; ++i) {
          //if (!m_burntvertices[i]) {
          glm::vec3 circle_center;

          if (geometry == "cap") {
              circle_center = glm::vec3(center[0],  m_vertices[i][1], center[2]);
              
              if (glm::length(m_vertices[i] - circle_center) < r) {
                      float to_cap = center[1] - m_vertices[i][1];
                      if (to_cap < 0.2 && to_cap > 0 && m_movable[i]) {
                          m_vertices[i] = m_vertices[i] + glm::vec3(0,1,0) * (to_cap);
                      }
                      else if (to_cap > 0 && m_movable[i]) {
                          m_vertices[i] = m_vertices[i] + glm::normalize(m_vertices[i] - circle_center) * (r-glm::length(m_vertices[i] - circle_center));
                      }
              }
              
              continue;
          }
          
          /*
          else if (geometry == "cube") {
              if (m_vertices[i][0] > center[0]-r && m_vertices[i][0] < center[0]+r && m_vertices[i][2] > center[2]-r && m_vertices[i][2] < center[2]+r) {
                  float to_cap = center[1] - m_vertices[i][1];
                  if (to_cap < 1 && to_cap > 0 && m_movable[i]) {
                      m_vertices[i] = m_vertices[i] + glm::vec3(0,1,0) * (to_cap);
                  }
              }
          }
          */
          
         
          else if (geometry == "ball") {
              circle_center = center;
              glm::vec3 to_center = m_vertices[i] - circle_center;
              float distance = glm::length(to_center);
              if (distance < r && m_movable[i]) {
                  m_vertices[i] = m_vertices[i] + glm::normalize(to_center) * (r-distance);
              }
          }
          //}
          //else { continue; }
          
      }
      
      /* self collision detection 
      for (int i = 0; i < N*N; ++i) {
          for (int j = 0; j < N*N; ++j) {
              glm::vec3 between_particals = m_vertices[i] - m_vertices[j];
              if (glm::length(between_particals) < 0.2) {
                  add_force(between_particals * 3000, i);
              }
          }
      }
      */
      
      /* Burning Starts */
      if (t > start) {
          //std::cout << "buring start" << std::endl;
          //update_constraints();
          burn_cloth(N);
      }
    
      
  }

    
  //m_temperature[N*N/2 + N/2] = 20000;
  /*for (int t = 0; t < 600; ++t) {
      burn_cloth(N);
  }
  */
      /*
      for (int i = 0; i < N*N; ++i) {
          if (m_temperature[i] >= pyrolysis) {
              ++m_timer[i];
          }
      }
      
      int count = 0;
      for (Triangle tri : m_faces) {
       
          if (std::max(std::max(m_temperature[tri.v1], m_temperature[tri.v2]), m_temperature[tri.v3]) >= pyrolysis) {
              m_burning[count] = true;
          }
          if (std::max(std::max(m_timer[tri.v1], m_timer[tri.v2]), m_timer[tri.v3]) >= burningtime) {
              //remove the cell
              m_burning[count] = false;
          }
          ++count;
      }
      
      
      float temp;
      
      temp = m_temperature[0];
      heat_transfer(0, 1);
      heat_transfer(0, N);
      heat_transfer(0, N+1);
      m_lasttemp[0] = temp;
      
      temp = m_temperature[N-1];
      heat_transfer(N-1, N-2);
      heat_transfer(N-1, 2*N-1);
      heat_transfer(N-1, 2*N-2);
      m_lasttemp[N-1] = temp;
      
      temp = m_temperature[N*(N-1)];
      heat_transfer(N*(N-1), N*(N-2));
      heat_transfer(N*(N-1), N*(N-2)+1);
      heat_transfer(N*(N-1), N*(N-1)+1);
      m_lasttemp[N*(N-1)] = temp;
      
      temp = m_temperature[N*N-1];
      heat_transfer(N*N-1, N*N-2);
      heat_transfer(N*N-1, N*(N-1)-1);
      heat_transfer(N*N-1, N*(N-1)-2);
      m_lasttemp[N*N-1] = temp;
      
      
      for (int i = 1; i < N-1; ++i) {
          for (int j = 1; j < N-1; ++j) {
              
              temp = m_temperature[j+i*N];
              heat_transfer(j+i*N, j+i*N-1);
              heat_transfer(j+i*N, j+i*N+1);
              heat_transfer(j+i*N, j+i*N-N-1);
              heat_transfer(j+i*N, j+i*N-N);
              heat_transfer(j+i*N, j+i*N-N+1);
              heat_transfer(j+i*N, j+i*N+N-1);
              heat_transfer(j+i*N, j+i*N+N);
              heat_transfer(j+i*N, j+i*N+N+1);
              m_lasttemp[j+i*N] = temp;
              
          }
      }
      
      for (int i = 1; i < N-1; ++i) {
          
          temp = m_temperature[i];
          heat_transfer(i, i-1);
          heat_transfer(i, i+1);
          heat_transfer(i, i+N-1);
          heat_transfer(i, i+N);
          heat_transfer(i, i+N+1);
          m_lasttemp[i] = temp;
          
      }
      
      for (int i = 1; i < N-1; ++i) {
          
          temp = m_temperature[i*N];
          heat_transfer(i*N, i*N-N);
          heat_transfer(i*N, i*N+N);
          heat_transfer(i*N, i*N-N+1);
          heat_transfer(i*N, i*N+1);
          heat_transfer(i*N, i*N+N+1);
          m_lasttemp[i*N] = temp;
          
      }
      
      for (int i = 1; i < N-1; ++i) {
          
          temp = m_temperature[i*N+N-1];
          heat_transfer(i*N+N-1, i*N-1);
          heat_transfer(i*N+N-1, i*N+2*N-1);
          heat_transfer(i*N+N-1, i*N-2);
          heat_transfer(i*N+N-1, i*N+N-2);
          heat_transfer(i*N+N-1, i*N+2*N-2);
          m_lasttemp[i*N+N-1] = temp;
          
      }
      
      for (int i = 1; i < N-1; ++i) {
          
          temp = m_temperature[i+N*(N-1)];
          heat_transfer(i+N*(N-1), i+N*(N-1)-1);
          heat_transfer(i+N*(N-1), i+N*(N-1)+1);
          heat_transfer(i+N*(N-1), i+N*(N-1)-N-1);
          heat_transfer(i+N*(N-1), i+N*(N-1)-N);
          heat_transfer(i+N*(N-1), i+N*(N-1)-N+1);
          m_lasttemp[i+N*(N-1)] = temp;
          
      }
      
     */
  

    

  for (Triangle tri : m_faces) {
    glm::vec3 n = glm::normalize(glm::cross(m_vertices[tri.v2] - m_vertices[tri.v1],
					    m_vertices[tri.v3] - m_vertices[tri.v1]));
    m_vnorms[tri.v1] += n;
    ++m_facecount[tri.v1];
    m_vnorms[tri.v2] += n;
    ++m_facecount[tri.v2];
    m_vnorms[tri.v3] += n;
    ++m_facecount[tri.v3];
  }

  for (int i = 0; i < N*N; ++i) {
    m_vnorms[i] /= m_facecount[i];
  }

  glm::vec3 cur_max = m_vertices[0];
  glm::vec3 cur_min = m_vertices[0];
  for (glm::vec3 v : m_vertices) {
    if (v[0] > cur_max[0]) {cur_max[0] = v[0];}
    if (v[1] > cur_max[1]) {cur_max[1] = v[1];}
    if (v[2] > cur_max[2]) {cur_max[2] = v[2];}
    if (v[0] < cur_min[0]) {cur_min[0] = v[0];}
    if (v[1] < cur_min[1]) {cur_min[1] = v[1];}
    if (v[2] < cur_min[2]) {cur_min[2] = v[2];}
  }
  m_max = cur_max;
  m_min = cur_min;
  
}



void Cloth::apply_force(int index) {
    //if (!m_burntvertices[index]) {
        m_acc[index] = m_forces[index] / m_masses[index];
        m_forces[index] = glm::vec3(0);
    //}
}

void Cloth::add_force(glm::vec3 force,int index) {
    //if (!m_burntvertices[index]) {
        m_forces[index] += force;
    //}
}

void Cloth::move_vertex(int index) {
  //if (!m_burntvertices[index]) {
      if (!m_movable[index]) { return; }
      // x_t+1 = x_t + v * t + a * t * t
      glm::vec3 temp = m_vertices[index];
      m_vertices[index] = m_vertices[index] + (m_vertices[index] - m_last[index])*(1.0 - damping) + m_acc[index] * TIME_STEP;
      m_last[index] = temp;
  //}
}

void Cloth::heat_transfer(int index1, int index2) {
    //if (m_burntvertices[index1] || m_burntvertices[index2] ) { return; }
        glm::vec3 x = m_vertices[index1] - m_vertices[index2];
        float cur_x = glm::length(x);
        m_temperature[index1] = m_temperature[index1] + diffusion * TIME_STEP * ((m_lasttemp[index2] - m_lasttemp[index1]) / (cur_x * cur_x));
    
}

void Cloth::check_constraints() {
    for (Constraint con : m_springs) {
        //if (con.initial > con.rest) { std::cout << "wrong!" << std::endl; }
        //if ((m_burntvertices[con.v1] || m_burntvertices[con.v2])  && (con.type == 'c')) { continue; }
        
        glm::vec3 x = m_vertices[con.v2] - m_vertices[con.v1];
        float cur_x = glm::length(x);
        
        if (m_burntvertices[con.v2]) {
            m_vertices[con.v2] -= x * (1 - con.initial/cur_x);
        }
        else if (m_burntvertices[con.v1]) {
            m_vertices[con.v1] += x * (1 - con.initial/cur_x);
        }
        else {
            if (m_movable[con.v1]) {
                if (m_movable[con.v2]) {
                    m_vertices[con.v1] += x * (1 - con.initial/cur_x) / 2;
                    m_vertices[con.v2] -= x * (1 - con.initial/cur_x) / 2;
                }
                else {
                    m_vertices[con.v1] += x * (1 - con.initial/cur_x);
                }
            } else {
                if (m_movable[con.v2]) {
                    m_vertices[con.v2] -= x * (1 - con.initial/cur_x);
                }
            }
        }
        
    }
}

void Cloth::update_constraints() {
    
    int count = 0;
    for (Constraint con : m_springs) {
        
        if (m_burningvertices[con.v1] || m_burningvertices[con.v2]) { //|| m_temperature[con.v1] > 100 || m_temperature[con.v2] > 100) {
            //(m_temperature[con.v1] > 100 || m_temperature[con.v2] > 100) {

        
        //float temp_change1 = std::min(pyrolysis - m_temperature[con.v1], 0.0f);
        //float temp_change2 = std::min(pyrolysis - m_temperature[con.v2], 0.0f);
        //float rate = 1 + 0.00000000001 * (temp_change1 + temp_change2);
        //m_springs[count].initial = m_springs[count].rest * rate;
        //if (m_temperature[con.v1] < 0 || m_temperature[con.v2] < 0) { std::cout << "wrong!" << std::endl; }
            float rate =  1 - 0.00001; //* (m_temperature[con.v1] + m_temperature[con.v2]) / 2;
            m_springs[count].initial = rate * m_springs[count].initial;
            
        }
            
        ++ count;
    }
    
    
}

void Cloth::burn_cloth(int N) {
    
    for (int i = 0; i < N*N; ++i) {
        //if (!m_burntvertices[i]) {
        if (m_temperature[i] >= pyrolysis) {
            m_burningvertices[i] = true;
        }
        
        if (m_burningvertices[i]) {
            ++m_timer[i];
        }
        
        if (m_timer[i] >= burningtime) {
            m_burntvertices[i] = true;
        }
        
        //}
    }
    
    int count = 0;
    for (Triangle tri : m_faces) {
        //if (!m_burntfaces[count]) {
            //if (std::max(std::max(m_temperature[tri.v1], m_temperature[tri.v2]), m_temperature[tri.v3]) >= pyrolysis) {
            if (m_burningvertices[tri.v1] || m_burningvertices[tri.v2] || m_burningvertices[tri.v3]) {
                m_burningfaces[count] = true;
            }
            //if (std::max(std::max(m_timer[tri.v1], m_timer[tri.v2]), m_timer[tri.v3]) >= burningtime) {
            if (m_burntvertices[tri.v1] || m_burntvertices[tri.v2] || m_burntvertices[tri.v3]) {
            //remove the cell
                //m_burning[count] = false;
                m_burntfaces[count] = true;
            //std::vector<Triangle>::iterator it = m_faces.begin() + count;
            //m_faces.erase(it);
            }
        //}
        ++count;
    }
    
    count = 0;
    /*
    for (Constraint con : m_springs) {
        //if (!m_burntsprings[count]) {
            if (m_burntvertices[con.v1] && m_burntvertices[con.v2]) {
                //m_burntsprings[count] = true;
                //m_springs[count].initial *= 2;
            }
        //}
        ++count;
    }
     */
    
    float temp;
    
    temp = m_temperature[0];
    heat_transfer(0, 1);
    heat_transfer(0, N);
    heat_transfer(0, N+1);
    m_lasttemp[0] = temp;
    
    temp = m_temperature[N-1];
    heat_transfer(N-1, N-2);
    heat_transfer(N-1, 2*N-1);
    heat_transfer(N-1, 2*N-2);
    m_lasttemp[N-1] = temp;
    
    temp = m_temperature[N*(N-1)];
    heat_transfer(N*(N-1), N*(N-2));
    heat_transfer(N*(N-1), N*(N-2)+1);
    heat_transfer(N*(N-1), N*(N-1)+1);
    m_lasttemp[N*(N-1)] = temp;
    
    temp = m_temperature[N*N-1];
    heat_transfer(N*N-1, N*N-2);
    heat_transfer(N*N-1, N*(N-1)-1);
    heat_transfer(N*N-1, N*(N-1)-2);
    m_lasttemp[N*N-1] = temp;
    
    
    for (int i = 1; i < N-1; ++i) {
        for (int j = 1; j < N-1; ++j) {
            
            temp = m_temperature[j+i*N];
            heat_transfer(j+i*N, j+i*N-1);
            heat_transfer(j+i*N, j+i*N+1);
            heat_transfer(j+i*N, j+i*N-N-1);
            heat_transfer(j+i*N, j+i*N-N);
            heat_transfer(j+i*N, j+i*N-N+1);
            heat_transfer(j+i*N, j+i*N+N-1);
            heat_transfer(j+i*N, j+i*N+N);
            heat_transfer(j+i*N, j+i*N+N+1);
            m_lasttemp[j+i*N] = temp;
            
        }
    }
    
    for (int i = 1; i < N-1; ++i) {
        
        temp = m_temperature[i];
        heat_transfer(i, i-1);
        heat_transfer(i, i+1);
        heat_transfer(i, i+N-1);
        heat_transfer(i, i+N);
        heat_transfer(i, i+N+1);
        m_lasttemp[i] = temp;
        
    }
    
    for (int i = 1; i < N-1; ++i) {
        
        temp = m_temperature[i*N];
        heat_transfer(i*N, i*N-N);
        heat_transfer(i*N, i*N+N);
        heat_transfer(i*N, i*N-N+1);
        heat_transfer(i*N, i*N+1);
        heat_transfer(i*N, i*N+N+1);
        m_lasttemp[i*N] = temp;
        
    }
    
    for (int i = 1; i < N-1; ++i) {
        
        temp = m_temperature[i*N+N-1];
        heat_transfer(i*N+N-1, i*N-1);
        heat_transfer(i*N+N-1, i*N+2*N-1);
        heat_transfer(i*N+N-1, i*N-2);
        heat_transfer(i*N+N-1, i*N+N-2);
        heat_transfer(i*N+N-1, i*N+2*N-2);
        m_lasttemp[i*N+N-1] = temp;
        
    }
    
    for (int i = 1; i < N-1; ++i) {
        
        temp = m_temperature[i+N*(N-1)];
        heat_transfer(i+N*(N-1), i+N*(N-1)-1);
        heat_transfer(i+N*(N-1), i+N*(N-1)+1);
        heat_transfer(i+N*(N-1), i+N*(N-1)-N-1);
        heat_transfer(i+N*(N-1), i+N*(N-1)-N);
        heat_transfer(i+N*(N-1), i+N*(N-1)-N+1);
        m_lasttemp[i+N*(N-1)] = temp;
        
    }
}

void Cloth::pin(int index) {
  m_movable[index] = false;
}

//void Cloth::burn(int index) {
  //m_burningfaces[index] = true;
//}


/* For Ray Tracing */

Intersection Cloth::Intersect(glm::vec3 orig, glm::vec3 dest) {
    Intersection I;
    
    if (!IntersectBound(orig, dest)) { return I;}
    //std::cout << "Calling Intersect" << std::endl;
    double cur_t = 0;
    glm::vec3 n;
    //std::cout << m_faces.size() << std::endl;
    int count = 0;
    int found = -1;
    for (Triangle tri : m_faces) {
        
        //if (m_burntfaces[count]) {continue;}
        
        glm::vec3 BA = m_vertices[tri.v2] - m_vertices[tri.v1];
        glm::vec3 CA = m_vertices[tri.v3] - m_vertices[tri.v1];
        glm::vec3 N = glm::normalize(glm::cross(BA, CA));
        float denom = glm::dot(N, N);
        
        //std::cout << "Calling Triangle Intersect" << std::endl;
        
        double t = IntersectTri(orig, dest,
                                m_vertices[tri.v1], m_vertices[tri.v2], m_vertices[tri.v3]);
        //std::cout << "\n\n\n" << t << std::endl;
        if (t != 0 && !m_burntfaces[count]) {
            
            if ((cur_t == 0) || (t < cur_t)) {
                cur_t = t;
                /*
                if (m_burning[count] == true) {
                    //std::cout << "burning" << std::endl;
                    I.burning = true;
                }
                */
                found = count;
                if (m_vnorms.size() != 0) {
                    glm::vec3 p = orig + cur_t * dest;
                    /*
                     glm::vec3 PA = p - m_vertices[tri.v1];
                     float d00 = glm::dot(BA, BA);
                     float d01 = glm::dot(BA, CA);
                     float d02 = glm::dot(BA, PA);
                     float d11 = glm::dot(CA, CA);
                     float d12 = glm::dot(PA, PA);
                     */
                    //float denom = 1 / (d00 * d11 - d01 * d01);
                    glm::vec3 CB = m_vertices[tri.v3] - m_vertices[tri.v2];
                    glm::vec3 PB = p - m_vertices[tri.v2];
                    float u = glm::dot(N, glm::cross(CB, PB));
                    
                    glm::vec3 AC = m_vertices[tri.v1] - m_vertices[tri.v3];
                    glm::vec3 PC = p - m_vertices[tri.v3];
                    float v = glm::dot(N,glm::cross(AC, PC));
                    
                    u /= denom;
                    v /= denom;
                    
                    // u + v + w = 1
                    //float v = (d11 * d02 - d01 * d12) * denom;
                    //float w = (d00 * d12 - d01 * d02) * denom;
                    float w = 1 - u - v;
                    
                    n = glm::normalize(u * m_vnorms[tri.v1] + v * m_vnorms[tri.v2] + w * m_vnorms[tri.v3]);
                    //std::cout << m_vnorms[tri.v1][0] << ", " << m_vnorms[tri.v1][1] << ", " << m_vnorms[tri.v1][2] << ", "  << std::endl;
                    //std::cout << u << ", " << v << ", " << w << std::endl;
                    //std::cout << n[0] << ", " << n[1] << ", " << n[2] << ", "  << std::endl;
                } else {
                    
                    n = glm::normalize(glm::cross(m_vertices[tri.v2] - m_vertices[tri.v1],
                                                  m_vertices[tri.v3]- m_vertices[tri.v1]));
                }
            }
        }
        ++count;
    }
    //I.burning = true;
    if (found != -1) {
        I.temperature = std::max(std::max(m_temperature[m_faces[found].v1], m_temperature[m_faces[found].v2]), m_temperature[m_faces[found].v3]);
        if (m_burningfaces[found] == true) {
            //std::cout << "burning" << std::endl;
            I.burning = true;
        }
    }
    I.t = cur_t;
    I.p = orig + (float)I.t * dest;
    /*
     if (m_vnoms.size() != 0) {
     glm::vec3 PA = I.p - m_vertices[tri.v1];
     float d00 = glm::dot(BA, BA);
     float d01 = glm::dot(BA, CA);
     float d02 = glm::dot(BA, PA);
     float d11 = glm::dot(CA, CA);
     float d12 = glm::ot(PA, PA);
     float denom = 1 / (d00 * d11 - d01 * d01);
     
     // u + v + w = 1
     float v = (d11 * d02 - d01 * d12) * dnom;
     float w = (d00 * d12 - d01 * d02) * dnom;
     float u = 1 - v - w;
     */
    
    I.n = n;
    I.kd = glm::vec3(0);
    I.ks = glm::vec3(0);
    
    return I;
}
