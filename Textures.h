//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h

#include <vector>
#define N_SQUARES 16
#define N_STRIPES 18
#define TEXTURE_REPETITION 5

#include "glm/glm.hpp"

typedef struct PNG_Image {
    unsigned width, height;
    std::vector<unsigned char> data;
} PNG_Image_t;


glm::vec3 checkerboardTexture(glm::vec2 uv, PNG_Image_t * image);
glm::vec3 rainbowTexture(glm::vec2 uv, PNG_Image_t * image);
glm::vec3 imageTexture( glm::vec2 uv, PNG_Image_t * image);


#endif /* Textures_h */
