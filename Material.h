//
//  Material.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef Material_h
#define Material_h

#include "glm/glm.hpp"
#include "Textures.h"

/**
 Structure describing a material of an object
 */
struct Material{
    glm::vec3 ambient = glm::vec3(0.0);
    glm::vec3 diffuse = glm::vec3(1.0);
    glm::vec3 specular = glm::vec3(0.0);
    float shininess = 0.0;
    float reflectiveness = 0.0;
    float refractiveness = 0.0;
    Texture * texture = nullptr;
};

#endif /* Material_h */
