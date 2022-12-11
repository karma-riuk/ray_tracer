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
#define TEXTURE_REPETITION 2

#include "glm/glm.hpp"

typedef struct PNG_Image {
    unsigned width, height;
    std::vector<unsigned char> data;
} PNG_Image_t;

class Texture {

  public:
    Texture(){};
    virtual glm::vec3 texture(glm::vec2 uv) = 0;
};

class ImageTexture : public Texture {
  public:
    PNG_Image_t &base_color, &height_map, &normals_map, &ambient_occlusion_map, &roughness_map;
    ImageTexture(PNG_Image_t & base_color, PNG_Image_t & height_map, PNG_Image_t & normals,
                 PNG_Image_t & ambient_occlusion, PNG_Image_t & roughness)
        : base_color(base_color), height_map(height_map), normals_map(normals),
          ambient_occlusion_map(ambient_occlusion), roughness_map(roughness) {}
    static glm::vec4 getRGBAat(glm::vec2 uv, PNG_Image_t & image);
    virtual glm::vec3 texture(glm::vec2 uv);
    virtual glm::vec3 normal(glm::vec2 uv);
    virtual glm::vec3 occlusion(glm::vec2 uv);
    virtual float roughness(glm::vec2 uv);
};

class RainbowTexture : public Texture {
    virtual glm::vec3 texture(glm::vec2 uv);
};

class CheckerBoardTexture : public Texture {
    virtual glm::vec3 texture(glm::vec2 uv);
};
#endif /* Textures_h */
