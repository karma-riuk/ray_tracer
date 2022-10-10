/**
@file main.cpp
*/

#include "glm/ext/scalar_constants.hpp"
#include "glm/glm.hpp"
#include "glm/gtx/string_cast.hpp"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "Image.h"
#include "Material.h"
#include "glm/gtx/projection.hpp"

using namespace std;

/**
 Class representing a single ray.
 */
class Ray {
  public:
    glm::vec3 origin;    ///< Origin of the ray
    glm::vec3 direction; ///< Direction of the ray
                         /**
                          Contructor of the ray
                          @param origin Origin of the ray
                          @param direction Direction of the ray
                          */
    Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction) {}
};

class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit {
    bool hit; ///< Boolean indicating whether there was or there was no intersection with an object
    glm::vec3 normal;       ///< Normal vector of the intersected object at the intersection point
    glm::vec3 intersection; ///< Point of Intersection
    float distance;         ///< Distance from the origin of the ray to the intersection point
    Object * object;        ///< A pointer to the intersected object
    glm::vec2 uv;           ///< Coordinates for computing the texture (texture coordinates)
};

/**
 General class for the object
 */
class Object {
  public:
    glm::vec3 color;   ///< Color of the object
    Material material; ///< Structure describing the material of the object
                       /** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;
};

float norm(glm::vec3 v) { return glm::distance(glm::vec3(0), v); }

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object {
  private:
    float radius;     ///< Radius of the sphere
    glm::vec3 center; ///< Center of the sphere

  public:
    /**
     The constructor of the sphere
     @param radius Radius of the sphere
     @param center Center of the sphere
     @param color Color of the sphere
     */
    Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center) {
        this->color = color;
    }
    Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center) {
        this->material = material;
    }
    /** Implementation of the intersection function*/
    Hit intersect(Ray ray) {

        glm::vec3 c = center - ray.origin;

        float cdotc = glm::dot(c, c);
        float cdotd = glm::dot(c, ray.direction);

        Hit hit;

        float D = 0;
        if (cdotc > cdotd * cdotd) {
            D = sqrt(cdotc - cdotd * cdotd);
        }
        if (D <= radius) {
            hit.hit = true;
            float t1 = cdotd - sqrt(radius * radius - D * D);
            float t2 = cdotd + sqrt(radius * radius - D * D);

            float t = t1;
            if (t < 0)
                t = t2;
            if (t < 0) {
                hit.hit = false;
                return hit;
            }

            hit.intersection = ray.origin + t * ray.direction;
            hit.normal = glm::normalize(hit.intersection - center);
            hit.distance = glm::distance(ray.origin, hit.intersection);
            hit.object = this;

            hit.uv = glm::vec2((glm::asin(hit.normal.y) + M_PI_2) / M_PI,
                               (glm::atan(hit.normal.z, hit.normal.x) + M_PI) / M_PI);


        } else {
            hit.hit = false;
        }
        return hit;
    }
};

class Plane : public Object {

  private:
    glm::vec3 point;
    glm::vec3 normal;

  public:
    Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal) {}
    Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal) {
        this->material = material;
    }

    Hit intersect(Ray ray) {
        Hit hit;
        hit.hit = false;

        float d_dot_N = glm::dot(ray.direction, normal);
        if (d_dot_N == 0)
            return hit;

        float po_dot_N = glm::dot(point - ray.origin, normal);

        float t = po_dot_N / d_dot_N;
        if (t < 0)
            return hit;

        hit.hit = true;
        hit.intersection = ray.origin + ray.direction * t;
        hit.normal = (d_dot_N > 0 ? -1.f : 1.f) * normal;
        hit.distance = glm::distance(ray.origin, hit.intersection);
        hit.object = this;
        glm::vec3 p1 = point;
        glm::vec3 p2 = glm::vec3(-hit.normal.y, hit.normal.x, hit.normal.z);
        glm::vec3 p3 = glm::cross(hit.normal, p2 - p1);

        glm::vec3 v1 = glm::normalize(p2-p1), v2 = glm::normalize(p3-p1);

        hit.uv = glm::vec2(
                norm(glm::proj(hit.intersection, v1)),
                norm(glm::proj(hit.intersection, v2))
                );
        return hit;
    }
};

/**
 Light class
 */
class Light {
  public:
    glm::vec3 position; ///< Position of the light source
    glm::vec3 color;    ///< Color/intentisty of the light source
    Light(glm::vec3 position) : position(position) { color = glm::vec3(1.0); }
    Light(glm::vec3 position, glm::vec3 color) : position(position), color(color) {}
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(1.0, 1.0, 1.0);
vector<Object *> objects; ///< A list of all objects in the scene

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction,
                     Material material) {

    glm::vec3 color(0.0);
    for (size_t light_num = 0; light_num < lights.size(); light_num++) {

        glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
        glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

        float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
        float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

        // cout << glm::to_string(uv) << endl;
        glm::vec3 diffuse;
        if (material.texture)
             diffuse = material.texture(uv) * glm::vec3(NdotL);
        else 
            diffuse = material.diffuse * glm::vec3(NdotL);
        glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));

        /*


         Excercise 3 - Modify the code by adding attenuation of the light due to distance from the
         intersection point to the light source



         */

        color += lights[light_num]->color * (diffuse + specular);
    }
    color += ambient_light * material.ambient;

    color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
    return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray) {

    Hit closest_hit;

    closest_hit.hit = false;
    closest_hit.distance = INFINITY;
    closest_hit.object = NULL;

    for (size_t k = 0; k < objects.size(); k++) {
        Hit hit = objects[k]->intersect(ray);
        if (hit.hit == true && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }

    glm::vec3 color(0.0);

    if (closest_hit.hit) {
        color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv,
                           glm::normalize(-ray.direction), closest_hit.object->material);
    } else {
        color = glm::vec3(0.0, 0.0, 0.0);
    }
    return color;
}
/**
 Function defining the scene
 */
void sceneDefinition() {

    Material green_diffuse;
    green_diffuse.ambient = glm::vec3(0.07f, 0.09f, 0.07f);
    green_diffuse.diffuse = glm::vec3(0.7f, 0.9f, 0.7f);

    Material red_specular;
    red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
    red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
    red_specular.specular = glm::vec3(0.5);
    red_specular.shininess = 10.0;

    Material blue_specular;
    blue_specular.ambient = glm::vec3(0.07f, 0.07f, 0.1f);
    blue_specular.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);
    blue_specular.specular = glm::vec3(0.6);
    blue_specular.shininess = 100.0;

    objects.push_back(new Sphere(1.0, glm::vec3(1, -2, 8), blue_specular));
    objects.push_back(new Sphere(0.5, glm::vec3(-1, -2.5, 6), red_specular));
    objects.push_back(new Sphere(1.0, glm::vec3(3, -2, 6), green_diffuse));

    // Excercise 2 - Textured sphere
    Material textured;
    textured.texture = &rainbowTexture;
    objects.push_back(new Sphere(7.0, glm::vec3(-6, 4, 23), textured));

    objects.push_back(new Plane(glm::vec3(0, 0, -.01), glm::vec3(0, 0, 1))); // plane behind
    objects.push_back(new Plane(glm::vec3(-15, 0, 0), glm::vec3(1, 0, 0)));  // plane left
    objects.push_back(new Plane(glm::vec3(15, 0, 0), glm::vec3(1, 0, 0)));   // plane right
    objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0, 1, 0)));   // plane below
    objects.push_back(new Plane(glm::vec3(0, 17, 0), glm::vec3(0, 1, 0)));   // plane above
    objects.push_back(new Plane(glm::vec3(0, 0, 30), glm::vec3(0, 0, 1)));   // plane in front

    lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.4)));
    lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.4)));
    lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity) {

    glm::vec3 tonemapped = intensity; // tonemapped intensity

    /*


     Excercise 3 - Tone mapping



     */

    return glm::clamp(tonemapped, glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char * argv[]) {

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; // width of the image
    int height = 768; // height of the image
    float fov = 90;   // field of view

    sceneDefinition(); // Let's define a scene

    Image image(width, height); // Create an image where we will store the result

    float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
    float X = -s * width / 2;
    float Y = s * height / 2;

    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++) {

            float dx = X + i * s + s / 2;
            float dy = Y - j * s - s / 2;
            float dz = 1;

            glm::vec3 origin(0, 0, 0);
            glm::vec3 direction(dx, dy, dz);
            direction = glm::normalize(direction);

            Ray ray(origin, direction);

            image.setPixel(i, j, toneMapping(trace_ray(ray)));
        }

    t = clock() - t;
    cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
    cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second."
         << endl;

    // Writing the final results of the rendering
    if (argc == 2) {
        image.writeImage(argv[1]);
    } else {
        image.writeImage("./result.ppm");
    }

    return 0;
}
