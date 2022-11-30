#include "Textures.h"
#include "glm/gtx/string_cast.hpp"
#include <iostream>

#define RED glm::vec3(1, 0, 0)
#define GREEN glm::vec3(0, 1, 0)
#define BLUE glm::vec3(0, 0, 1)

glm::vec3 CheckerBoardTexture::texture(glm::vec2 uv) {
    return glm::vec3((int)(glm::floor(N_SQUARES * uv.x) + glm::floor(N_SQUARES * uv.y)) % 2);
}

glm::vec3 RainbowTexture::texture(glm::vec2 uv) {
    switch ((int)glm::floor(N_STRIPES * (uv.x + uv.y)) % 3) {
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
glm::vec3 ImageTexture::texture(glm::vec2 uv) {
    int x =
        (int)(round(TEXTURE_REPETITION * uv.x * this->base_color.width)) % (this->base_color.width);
    int y = (int)(round(TEXTURE_REPETITION * uv.y * this->base_color.height)) %
            (this->base_color.height);
    int index = (y * this->base_color.width + x) * 4;

    glm::vec3 ret = glm::vec3((float)this->base_color.data[index + 0] / 255,
                              (float)this->base_color.data[index + 1] / 255,
                              (float)this->base_color.data[index + 2] / 255);
    return ret;
}
