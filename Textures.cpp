#include "Textures.h"
#include "glm/gtx/string_cast.hpp"
#include <iostream>

#define RED glm::vec3(1, 0, 0)
#define GREEN glm::vec3(0, 1, 0)
#define BLUE glm::vec3(0, 0, 1)

glm::vec3 checkerboardTexture(glm::vec2 uv){
    return glm::vec3(
            (int) (glm::floor(N_SQUARES * uv.x) + glm::floor(N_SQUARES * uv.y)) % 2
            );
}

glm::vec3 rainbowTexture(glm::vec2 uv){
    switch ( (int) glm::floor(N_STRIPES * (uv.x + uv.y)) % 3) {
        case 0:
            return RED;
        case 1:
            return GREEN;
        case 2:
            return BLUE;
        default:
            return glm::vec3(0);
    }
}
