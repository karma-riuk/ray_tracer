//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h


#include "glm/glm.hpp"

glm::vec3 checkerboardTexture(glm::vec2 uv){
    float n = 20;
    float value = int(floor(n*uv.s) + floor(2*n*uv.t)) % 2;
    return glm::vec3(value);
}
glm::vec3 rainbowTexture(glm::vec2 uv){
    float n = 40;
    int value = int(floor(n*uv.t + 0.5*n*uv.s )) % 3;
    switch(value){
        case 0: return glm::vec3(1.0, 0.0, 0.0);
            break;
        case 1: return glm::vec3(0.0, 1.0, 0.0);
            break;
        default: return glm::vec3(0.0, 0.0, 1.0);
    }
}

#endif /* Textures_h */
