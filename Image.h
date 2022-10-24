//
/**
@file Image.h
*/


//  Image.h
//  Raytracer
//
//  Created by Piotr Didyk on 14.07.21.
//

#ifndef Image_h
#define Image_h
#include "glm/glm.hpp"

using namespace std;

/**
 Class allowing for creating an image and writing it to a file
 */
class Image{

private:
    int width, height; ///< width and height of the image
    int *data; ///< a pointer to the data representing the images
    
public:
    /**
     @param width with of the image
     @param height height of the image
     */
    Image(int width, int height);
    
    /**
     Writes and image to a file in ppm format
     @param path the path where to the target image
     */
    void writeImage(const char *path);
    
    /**
     Set a value for one pixel
     @param x x coordinate of the pixel - index of the column counting from left to right
     @param y y coordinate of the pixel - index of the row counting from top to bottom
     @param r red chanel value in range from 0 to 255
     @param g green chanel value in range from 0 to 255
     @param b blue chanel value in range from 0 to 255
     */
    void setPixel(int x, int y, int r, int g, int b);
    
    /**
     Set a value for one pixel
     @param x x coordinate of the pixel - index of the column counting from left to right
     @param y y coordinate of the pixel - index of the row counting from top to bottom
     @param r red chanel value in range from 0 to 1
     @param g green chanel value in range from 0 to 1
     @param b blue chanel value in range from 0 to 1
     */
    void setPixel(int x, int y, float r, float g, float b);
    
    /**
     Set a value for one pixel
     @param x x coordinate of the pixel - index of the column counting from left to right
     @param y y coordinate of the pixel - index of the row counting from top to bottom
     @param color color of the pixel expressed as vec3 of RGB values in range from 0 to 1
     */
    void setPixel(int x, int y, glm::vec3 color);
};

#endif /* Image_h */
