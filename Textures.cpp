#include "Textures.h"
#include "glm/gtx/string_cast.hpp"
#include <iostream>

#define RED glm::vec3(1, 0, 0)
#define GREEN glm::vec3(0, 1, 0)
#define BLUE glm::vec3(0, 0, 1)

glm::vec3 checkerboardTexture(glm::vec2 uv, PNG_Image_t * image){
    return glm::vec3(
            (int) (glm::floor(N_SQUARES * uv.x) + glm::floor(N_SQUARES * uv.y)) % 2
            );
}

glm::vec3 rainbowTexture(glm::vec2 uv, PNG_Image_t * image){
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

// image is the result of reading the png file, it contains 3 fields:
// - the width and height of the image (number of pixels)
// - a one dimentional vector containing the color of the pixels in format RBGA RGBA ...
// (hence the length of the vector is 4 * width * height)
glm::vec3 imageTexture(glm::vec2 uv, PNG_Image_t * image){
    int x = (int) (round(TEXTURE_REPETITION * uv.x * image->width)) % ( image->width );
    int y = (int) (round(TEXTURE_REPETITION * uv.y * image->height)) % (image->height);
    int index = (y * image->width + x) * 4;
    // printf("x: %d, y: %d\n", x/4, y/4);

    glm::vec3 ret = glm::vec3(
            (float) image->data[index + 0]/255,
            (float) image->data[index + 1]/255,
            (float) image->data[index + 2]/255
            );
    return ret;
}
