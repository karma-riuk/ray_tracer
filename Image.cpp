#include "Image.h"
#include <fstream>

Image::Image(int width, int height) : width(width), height(height) {
    data = new int[3 * width * height];
}

void Image::writeImage(const char * path) {
    std::ofstream file;
    file.open(path);
    file << "P3" << std::endl;
    file << width << " " << height << std::endl;
    file << 255 << std::endl;
    for (int h = 0; h < height; h++) {
        for (int w = 0; w < width; w++) {
            file << data[3 * (h * width + w)] << " ";
            file << data[3 * (h * width + w) + 1] << " ";
            file << data[3 * (h * width + w) + 2] << "  ";
        }
        file << std::endl;
    }
    file.close();
}
void Image::setPixel(int x, int y, int r, int g, int b) {
    data[3 * (y * width + x)] = r;
    data[3 * (y * width + x) + 1] = g;
    data[3 * (y * width + x) + 2] = b;
}

void Image::setPixel(int x, int y, float r, float g, float b) {
    data[3 * (y * width + x)] = (float)(255 * r);
    data[3 * (y * width + x) + 1] = (float)(255 * g);
    data[3 * (y * width + x) + 2] = (float)(255 * b);
}
void Image::setPixel(int x, int y, glm::vec3 color){
        data[3 * (y * width + x)] = (float)(255 * color.r);
        data[3 * (y * width + x) + 1] = (float)(255 * color.g);
        data[3 * (y * width + x) + 2] = (float)(255 * color.b);
}
