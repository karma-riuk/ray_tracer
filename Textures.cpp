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

// returns a vec4 of numbers in [0; 1] representing RGBA values
glm::vec4 ImageTexture::getRGBAat(glm::vec2 uv, PNG_Image_t & image) {
    int x = (int)(round(TEXTURE_REPETITION * uv.x * image.width)) % (image.width);
    int y = (int)(round(TEXTURE_REPETITION * uv.y * image.height)) % (image.height);
    int index = (y * image.width + x) * 4;

    // printf("%d, %d, %d, %d\n", image.data[index + 0], image.data[index + 1], image.data[index + 2], image.data[index + 3]);
    return glm::vec4((float)image.data[index + 0], (float)image.data[index + 1],
                     (float)image.data[index + 2], (float)image.data[index + 3]) / 255.f;
}

// image is the result of reading the png file, it contains 3 fields:
// - the width and height of the image (number of pixels)
// - a one dimentional vector containing the color of the pixels in format RBGA RGBA ...
// (hence the length of the vector is 4 * width * height)
glm::vec3 ImageTexture::texture(glm::vec2 uv) {
    return glm::pow(ImageTexture::getRGBAat(uv, this->base_color), glm::vec4(2.2f));
}

glm::vec3 ImageTexture::normal(glm::vec2 uv) {
    return ImageTexture::getRGBAat(uv, this->normals_map);
}

glm::vec3 ImageTexture::occlusion(glm::vec2 uv) {
    return ImageTexture::getRGBAat(uv, this->ambient_occlusion_map);
}

float ImageTexture::roughness(glm::vec2 uv) {
    return ImageTexture::getRGBAat(uv, this->roughness_map).x;
}
