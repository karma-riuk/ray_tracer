#include "Textures.h"

glm::vec3 checkerboardTexture(glm::vec2 uv){
    return glm::vec3(
            (int) (glm::floor(N_SQUARES * uv.x) + glm::floor(N_SQUARES * uv.y)) % 2
            );
}
glm::vec3 rainbowTexture(glm::vec2 uv){
    return glm::vec3(0);
}
