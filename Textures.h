//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h

#define N_SQUARES 16
#define N_STRIPES 18

#include "glm/glm.hpp"

glm::vec3 checkerboardTexture(glm::vec2 uv);
glm::vec3 rainbowTexture(glm::vec2 uv);

#endif /* Textures_h */
