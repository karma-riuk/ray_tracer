/**
@file main.cpp
*/
#include "glm/glm.hpp"
#include "glm/gtx/string_cast.hpp"
#include "glm/gtx/transform.hpp"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Image.h"
#include "Material.h"
#include "lib/lodepng/lodepng.h"
// #include "obj_reader.h"

using std::stringstream;

using namespace std;
#define EPS 0.0001f

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
    bool from_outside = true;
};

/**
 General class for the object
 */
class Object {

  protected:
    glm::mat4 transformationMatrix; ///< Matrix representing the transformation from the local to
                                    ///< the global coordinate system
    glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the
                                           ///< global to the local coordinate system
    glm::mat4 normalMatrix; ///< Matrix for transforming normal vectors from the local to the global
                            ///< coordinate system

  public:
    glm::vec3 color;   ///< Color of the object
    Material material; ///< Structure describing the material of the object
                       /** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;

    virtual glm::vec3 getMinCoords() = 0;
    virtual glm::vec3 getMaxCoords() = 0;
    /** Function that returns the material struct of the object*/
    Material getMaterial() { return material; }
    /** Function that set the material
     @param material A structure describing the material of the object
    */
    virtual void setMaterial(Material material) { this->material = material; }
    /** Functions for setting up all the transformation matrices
     @param matrix The matrix representing the transformation of the object in the global
     coordinates */
    void setTransformation(glm::mat4 matrix) {

        transformationMatrix = matrix;

        inverseTransformationMatrix = glm::inverse(matrix);
        normalMatrix = (glm::abs(glm::determinant(matrix)) == 1)
                           ? matrix
                           : glm::transpose(inverseTransformationMatrix);
    }
};

template <class DstType, class SrcType> bool isType(const SrcType * src) {
    return dynamic_cast<const DstType *>(src) != nullptr;
}

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
    Sphere(glm::vec3 color) : radius(1), center(glm::vec3(0)) { this->color = color; }
    Sphere(Material material) : radius(1), center(glm::vec3(0)) { this->material = material; }

    glm::vec3 getMinCoords() { return transformationMatrix * glm::vec4(-1, -1, -1, 1); }
    glm::vec3 getMaxCoords() { return transformationMatrix * glm::vec4(1, 1, 1, 1); }

    /** Implementation of the intersection function*/
    Hit intersect(Ray ray) {

        glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1);
        glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0);

        Ray local_ray(o, d);
        Hit hit;
        hit.hit = false;

        float a, b, c;
        a = pow(d.x, 2) + pow(d.y, 2) + pow(d.z, 2);
        b = 2 * (d.x * o.x + d.y * o.y + d.z * o.z);
        c = pow(o.x, 2) + pow(o.y, 2) + pow(o.z, 2) - 1;

        float delta = pow(b, 2) - 4 * a * c;
        if (delta < 0)
            return hit;

        float t_1, t_2;
        t_1 = (-b + sqrt(delta)) / (2 * a);
        t_2 = (-b - sqrt(delta)) / (2 * a);

        if (t_1 < 0 && t_2 < 0)
            return hit;
        glm::vec3 i1, i2;
        i1 = o + t_1 * d;
        i2 = o + t_2 * d;

        float t;
        t = min(t_1 < 0 ? INFINITY : t_1, t_2 < 0 ? INFINITY : t_2); // find lowest positive number

        glm::vec3 i = o + t * d;
        hit.intersection = i;
        hit.normal = i;
        hit.uv.x = (atan2(i.z, i.x) + M_PI) / M_PI;
        hit.uv.y = (asin(i.y) + M_PI / 2) / M_PI;

        bool from_inside = glm::all(glm::lessThan(glm::abs(o), glm::vec3(1)));
        hit.from_outside = !from_inside;

        // Retransform in global coordinates
        hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1);

        if (isType<ImageTexture>(material.texture)) {
            glm::vec3 tangent = glm::vec3(sin(hit.uv.x), 0, cos(hit.uv.x));
            // glm::vec3 tangent = glm::vec3(sin(hit.intersection.x), 0, cos(hit.intersection.z));
            glm::vec3 bitangent = glm::cross(hit.normal, tangent);

            glm::mat3 tbn(glm::normalize(tangent), glm::normalize(bitangent),
                          glm::normalize(hit.normal));

            ImageTexture * texture = (ImageTexture *)material.texture;
            glm::vec3 tmp_normal = texture->normal(hit.uv);
            tmp_normal = (2.f * tmp_normal) - glm::vec3(1);
            tmp_normal = tbn * tmp_normal;
            hit.normal = normalMatrix * glm::vec4(tmp_normal, 0);
        } else {
            hit.normal = normalMatrix * glm::vec4(hit.normal, 0);
        }
        hit.normal = glm::normalize(hit.normal);

        hit.distance = glm::distance(ray.origin, hit.intersection);

        hit.object = this;
        hit.hit = true;

        return hit;
    }
};

class Plane : public Object {

  private:
    glm::vec3 point;
    glm::vec3 b1, b2;
    glm::vec3 normal;
    static glm::vec3 computeB1(glm::vec3 n) {
        glm::vec3 a = glm::cross(n, glm::vec3(1, 0, 0));
        glm::vec3 b = glm::cross(n, glm::vec3(0, 1, 0));
        glm::vec3 max_a_b = glm::dot(a, a) > glm::dot(b, b) ? a : b;
        glm::vec3 c = glm::cross(n, glm::vec3(0, 0, 1));
        return glm::normalize(glm::dot(max_a_b, max_a_b) < glm::dot(c, c) ? c : max_a_b);
    }

  public:
    void computeBaseVectors() {
        b1 = .05f * computeB1(normal);
        b2 = .1f * glm::cross(normal, b1);
    }

    Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal) {
        computeBaseVectors();
    }
    Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal) {
        computeBaseVectors();
        this->material = material;
    }

    glm::vec3 getMaxCoords() { return point; }
    glm::vec3 getMinCoords() { return point; }
    Hit intersect(Ray ray) {

        Hit hit;
        hit.hit = false;
        float DdotN = glm::dot(ray.direction, normal);
        if (DdotN < 0) {

            float PdotN = glm::dot(point - ray.origin, normal);
            float t = PdotN / DdotN;

            if (t > 0) {
                hit.hit = true;
                hit.normal = normal;
                hit.distance = t;
                hit.object = this;
                hit.intersection = t * ray.direction + ray.origin;
                hit.uv.x = glm::dot(b1, hit.intersection);
                hit.uv.y = glm::dot(b2, hit.intersection);

                if (isType<ImageTexture>(material.texture)) {
                    glm::mat3 tbn(b1, b2, normal);
                    ImageTexture * texture = (ImageTexture *)material.texture;
                    glm::vec3 tmp_normal = texture->normal(hit.uv);
                    tmp_normal = (2.f * tmp_normal) - glm::vec3(1);
                    tmp_normal = tbn * tmp_normal;
                    hit.normal = tmp_normal;
                }
                hit.normal = glm::normalize(hit.normal);
                // glm::mat3 base(b1, b2, normal);
                // glm::vec3 ts = glm::inverse(base) * (hit.intersection - point);
                // hit.uv.x = ts.x;
                // hit.uv.y = ts.y;
            }
        }
        return hit;
    }
};

class Cone : public Object {
    Plane base;

  public:
    Cone(Material material) : base(glm::vec3(0, 1, 0), glm::vec3(0, 1, 0), material) {
        this->material = material;
    }
    glm::vec3 getMinCoords() { return transformationMatrix * glm::vec4(-1, 0, -1, 1); }
    glm::vec3 getMaxCoords() { return transformationMatrix * glm::vec4(1, 1, 1, 1); }
    Hit intersect(Ray ray) {

        Hit hit;
        hit.hit = false;

        glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0);
        glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1);
        Ray local_ray(o, d);

        float a, b, c; // constant to solve at^2 + bt + c = 0
        a = pow(d.x, 2) - pow(d.y, 2) + pow(d.z, 2);
        b = 2 * (d.x * o.x - d.y * o.y + d.z * o.z);
        c = pow(o.x, 2) - pow(o.y, 2) + pow(o.z, 2);

        float delta = pow(b, 2) - 4 * a * c;
        if (delta < 0)
            return hit;

        float t_1, t_2;
        t_1 = (-b + sqrt(delta)) / (2 * a);
        t_2 = (-b - sqrt(delta)) / (2 * a);

        if (t_1 < 0 && t_2 < 0)
            return hit;
        glm::vec3 i1, i2;
        i1 = o + t_1 * d;
        i2 = o + t_2 * d;
        if ((i1.y < 0 && i2.y < 0) || (i1.y > 1 && i2.y > 1))
            return hit;

        float t;
        t = min(t_1 < 0 ? INFINITY : t_1, t_2 < 0 ? INFINITY : t_2); // find lowest positive number

        glm::vec3 i;
        i = o + t * d;
        if (i.y > 1) {
            hit = base.intersect(local_ray);
        } else {
            hit.hit = true;
            hit.object = this;
            hit.intersection = i;
            float y = (pow(i.x, 2) + pow(i.z, 2)) / i.y + i.y;

            glm::vec3 normal = i - glm::vec3(0, y, 0);
            hit.normal = normal;
            hit.uv.x = glm::atan(i.x, i.z);
            hit.uv.y = i.y;
        }

        // Retransform in global coordinates
        hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1);
        hit.normal = normalMatrix * glm::vec4(hit.normal, 0);
        hit.normal = glm::normalize(hit.normal);
        hit.distance = glm::distance(ray.origin, hit.intersection);

        return hit;
    }
};

class Fragment : public Object {};

class Triangle : public Fragment {

  private:
    glm::vec3 p1, p2, p3;
    glm::vec3 n1, n2, n3;
    glm::vec2 uv1, uv2, uv3;

  public:
    Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) : p1(p1), p2(p2), p3(p3) {
        n1 = n2 = n3 = glm::normalize(cross(p3 - p1, p2 - p1));
    }

    Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 n1, glm::vec3 n2, glm::vec3 n3)
        : p1(p1), p2(p2), p3(p3), n1(n1), n2(n2), n3(n3) {}
    Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 n1, glm::vec3 n2, glm::vec3 n3,
             glm::vec2 uv1, glm::vec2 uv2, glm::vec2 uv3)
        : p1(p1), p2(p2), p3(p3), n1(n1), n2(n2), n3(n3), uv1(uv1), uv2(uv2), uv3(uv3) {}
    Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 n1, glm::vec3 n2, glm::vec3 n3,
             Material material)
        : p1(p1), p2(p2), p3(p3), n1(n1), n2(n2), n3(n3) {
        ((Fragment *)this)->material = material;
    }

    glm::vec3 getMinCoords() {
        float min_x = min(p1.x, min(p2.x, p3.x));
        float min_y = min(p1.y, min(p2.y, p3.y));
        float min_z = min(p1.z, min(p2.z, p3.z));
        return glm::vec3(min_x, min_y, min_z);
    }
    glm::vec3 getMaxCoords() {
        float max_x = max(p1.x, max(p2.x, p3.x));
        float max_y = max(p1.y, max(p2.y, p3.y));
        float max_z = max(p1.z, max(p2.z, p3.z));
        return glm::vec3(max_x, max_y, max_z);
    }

    void set_normals(glm::vec3 n1, glm::vec3 n2, glm::vec3 n3) {
        this->n1 = n1;
        this->n2 = n2;
        this->n3 = n3;
    }

    Hit intersect(Ray ray) {
        Hit hit;
        hit.distance = INFINITY;
        hit.hit = false;

        glm::vec3 sn = glm::cross(p2 - p1, p3 - p1);
        Plane plane_from_triangle(p1, sn);
        Hit plane_hit = plane_from_triangle.intersect(ray);
        if (!plane_hit.hit)
            return hit;

        glm::vec3 p = plane_hit.intersection;
        glm::vec3 sn1, sn2, sn3; // normals for barycentric coords
        sn1 = glm::cross(p2 - p, p3 - p);
        sn2 = glm::cross(p3 - p, p1 - p);
        sn3 = glm::cross(p1 - p, p2 - p);

        float sign1, sign2, sign3;
        sign1 = glm::sign(glm::dot(sn1, sn));
        sign2 = glm::sign(glm::dot(sn2, sn));
        sign3 = glm::sign(glm::dot(sn3, sn));
        if (sign1 < 0 || sign2 < 0 || sign3 < 0) {
            return hit; // not hit
        }

        hit.hit = true;
        hit.intersection = p;

        hit.normal = glm::normalize(interpolate(sn1, sn2, sn3, n1, n2, n3));
        // hit.normal = glm::normalize(n1);
        hit.distance = plane_hit.distance;
        hit.object = (Fragment *)this;
        float W = length(sn);

        if (glm::length(uv1) > 0) {
            glm::vec3 interpolated = interpolate(sn1 / W, sn2 / W, sn3 / W, glm::vec3(uv1, 0),
                                                 glm::vec3(uv2, 0), glm::vec3(uv3, 0));
            hit.uv = interpolated;
            // hit.uv = uv1;
        } else
            hit.uv = glm::vec2(length(sn1) / W, length(sn2) / W);
        hit.from_outside = dot(hit.normal, ray.direction) < 0;
        // hit.from_outside = true;

        return hit;
    }

    glm::vec3 interpolate(glm::vec3 sn1, glm::vec3 sn2, glm::vec3 sn3, glm::vec3 v1, glm::vec3 v2,
                          glm::vec3 v3) {
        return length(sn1) * v1 + length(sn2) * v2 + length(sn3) * v3;
    }
};

class Diamond : public Fragment {
  private:
    Triangle *t1, *t2;

  public:
    Diamond(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 p4) {
        t1 = new Triangle(p1, p2, p4);
        t2 = new Triangle(p2, p3, p4);
    }
    Diamond(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 p4, glm::vec3 n1, glm::vec3 n2,
            glm::vec3 n3, glm::vec3 n4) {
        t1 = new Triangle(p1, p2, p4, n1, n2, n4);
        t2 = new Triangle(p2, p3, p4, n2, n3, n4);
    }
    Diamond(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 p4, glm::vec3 n1, glm::vec3 n2,
            glm::vec3 n3, glm::vec3 n4, glm::vec2 uv1, glm::vec2 uv2, glm::vec2 uv3,
            glm::vec2 uv4) {
        t1 = new Triangle(p1, p2, p4, n1, n2, n4, uv1, uv2, uv4);
        t2 = new Triangle(p2, p3, p4, n2, n3, n4, uv2, uv3, uv4);
    }

    glm::vec3 getMinCoords() override {
        glm::vec3 min1 = t1->getMinCoords();
        glm::vec3 min2 = t2->getMinCoords();
        float min_x = min(min1.x, min2.x);
        float min_y = min(min1.y, min2.y);
        float min_z = min(min1.z, min2.z);
        return glm::vec3(min_x, min_y, min_z);
    }
    glm::vec3 getMaxCoords() override {
        glm::vec3 max1 = t1->getMaxCoords();
        glm::vec3 max2 = t2->getMaxCoords();
        float max_x = max(max1.x, max2.x);
        float max_y = max(max1.y, max2.y);
        float max_z = max(max1.z, max2.z);
        return glm::vec3(max_x, max_y, max_z);
    }

    void setMaterial(Material material) override {
        t1->setMaterial(material);
        t2->setMaterial(material);
    }

    Hit intersect(Ray ray) override {
        Hit hit;
        hit.hit = false;
        Hit t1_hit = t1->intersect(ray);
        Hit t2_hit = t2->intersect(ray);
        if (!(t1_hit.hit || t2_hit.hit)) {
            return hit;
        }
        hit = (t1_hit.distance < t2_hit.distance) ? t1_hit : t2_hit;
        return hit;
    }
};

class Mesh : public Object {

  private:
    std::vector<Fragment *> fragments;

  public:
    Mesh(std::vector<Fragment *> fragments) : fragments(fragments) {}

    void setMaterial(Material material) {
        for (auto fragment : fragments) {
            fragment->setMaterial(material);
        }
    }

    glm::vec3 getMinCoords() {
        float min_x = INFINITY;
        float min_y = INFINITY;
        float min_z = INFINITY;
        for (auto fragment : fragments) {
            // triangle->setTransformation(transformationMatrix);
            glm::vec3 fragment_min = fragment->getMinCoords();
            // triangle->setTransformation(glm::translate(glm::vec3(0)));
            min_x = min(min_x, fragment_min.x);
            min_y = min(min_y, fragment_min.y);
            min_z = min(min_z, fragment_min.z);
        }

        return glm::vec3(min_x, min_y, min_z);
    }
    glm::vec3 getMaxCoords() {
        float max_x = -INFINITY;
        float max_y = -INFINITY;
        float max_z = -INFINITY;
        for (auto fragment : fragments) {
            // triangle->setTransformation(transformationMatrix);
            glm::vec3 fragment_max = fragment->getMaxCoords();
            // triangle->setTransformation(glm::translate(glm::vec3(0)));
            max_x = max(max_x, fragment_max.x);
            max_y = max(max_y, fragment_max.y);
            max_z = max(max_z, fragment_max.z);
        }
        return glm::vec3(max_x, max_y, max_z);
    }

    Hit intersect(Ray ray) {
        Hit hit;
        hit.hit = false;
        hit.distance = INFINITY;

        glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0);
        glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1);
        Ray local_ray(o, d);

        for (auto fragment : fragments) {
            Hit tr_hit = fragment->intersect(local_ray);
            if (tr_hit.hit && tr_hit.distance < hit.distance) {
                hit = tr_hit;
            }
        }

        // Retransform in global coordinates
        hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1);
        hit.normal = normalMatrix * glm::vec4(hit.normal, 0);
        hit.normal = glm::normalize(hit.normal);
        hit.distance = glm::distance(ray.origin, hit.intersection);
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

class Box : public Object {
  private:
    glm::vec3 v1, v2;
    Object * obj;
    std::vector<Triangle> triangles;

  public:
    Box(glm::vec3 v1, glm::vec3 v2, Object * obj) : v1(v1), v2(v2), obj(obj) {
        glm::vec3 A(v1);
        glm::vec3 B(v2.x, v1.y, v1.z);
        glm::vec3 C(v2.x, v2.y, v1.z);
        glm::vec3 D(v1.x, v2.y, v1.z);

        glm::vec3 G(v2);
        glm::vec3 E(v1.x, v1.y, v2.z);
        glm::vec3 F(v2.x, v1.y, v2.z);
        glm::vec3 H(v1.x, v2.y, v2.z);
        triangles.push_back(Triangle(B, A, D));
        triangles.push_back(Triangle(C, B, D));

        triangles.push_back(Triangle(F, B, C));
        triangles.push_back(Triangle(G, F, C));

        triangles.push_back(Triangle(E, F, G));
        triangles.push_back(Triangle(H, E, G));

        triangles.push_back(Triangle(A, E, H));
        triangles.push_back(Triangle(D, A, H));

        triangles.push_back(Triangle(C, D, H));
        triangles.push_back(Triangle(G, C, H));

        triangles.push_back(Triangle(F, E, A));
        triangles.push_back(Triangle(B, F, A));
    }

    glm::vec3 getMinCoords() { return transformationMatrix * glm::vec4(v1, 1); }
    glm::vec3 getMaxCoords() { return transformationMatrix * glm::vec4(v2, 1); }
    Hit intersect(Ray ray) {
        Hit hit;
        hit.hit = false;
        hit.distance = INFINITY;

        glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0);
        glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1);
        Ray local_ray(o, d);

        Hit tmp_hit;
        for (auto triangle : triangles) {
            // triangle.setMaterial(material);
            tmp_hit = triangle.intersect(local_ray);
            // if (tmp_hit.hit && tmp_hit.distance < hit.distance)
            if (tmp_hit.hit)
                hit = tmp_hit;
        }
        // return hit;
        return hit.hit ? obj->intersect(ray) : hit;
    };
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.001, 0.001, 0.001);
vector<Object *> objects; ///< A list of all objects in the scene
glm::vec3 tray_ray(Ray ray);
Hit find_closest_hit(Ray ray);

glm::vec3 refract(glm::vec3 i, glm::vec3 n, float index) {
    glm::vec3 a = n * glm::dot(i, n);
    glm::vec3 b = i - a;
    float beta = 1 / index;
    float alpha = sqrt(1 + (1 - pow(beta, 2)) * (glm::dot(b, b) / glm::dot(a, a)));
    return alpha * a + beta * b;
}

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction,
                     Material material) {
    glm::vec3 refractive_component(0);

    glm::vec3 color = ambient_light * material.ambient;

    for (size_t light_num = 0; light_num < lights.size(); light_num++) {

        glm::vec3 light_position = lights[light_num]->position;
        glm::vec3 light_direction = glm::normalize(light_position - point);
        Ray light_ray(point + (EPS * light_direction), light_direction);
        Hit closest_hit = find_closest_hit(light_ray);
        // TODO: missing the fact that refracting objects can let light through,
        // so there wouldn't be a complete shadow
        if (closest_hit.distance < glm::distance(point, light_position))
            continue;
        glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

        float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
        float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

        glm::vec3 diffuse_color = material.diffuse;
        if (material.texture) {
            diffuse_color = material.texture->texture(uv);
        }

        glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
        float shininess = material.shininess;

        if (isType<ImageTexture>(material.texture)) {
            ImageTexture * texture = (ImageTexture *)material.texture;
            diffuse *= 5.f;
            shininess = .5 / pow(texture->roughness(uv), 4) - .5;
            // shininess = glm::clamp(shininess, 0.f, 1.f);
        }

        glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, shininess));
        // glm::vec3 specular = glm::vec3(0);

        // distance to the light
        float r = glm::distance(point, lights[light_num]->position);
        r = max(r, 0.1f);

        color += lights[light_num]->color * (diffuse + specular) * 1.3f / r / r;
    }

    if (isType<ImageTexture>(material.texture))
        color *= ((ImageTexture *)material.texture)->occlusion(uv);
    else
        color *= 4;
    color *= 2;

    return color;
}

Hit find_closest_hit(Ray ray) {
    Hit closest_hit;

    closest_hit.hit = false;
    closest_hit.distance = INFINITY;

    for (size_t k = 0; k < objects.size(); k++) {
        Hit hit = objects[k]->intersect(ray);
        if (hit.hit == true && hit.from_outside && hit.distance > 0.01f && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }
    return closest_hit;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray) {

    Hit closest_hit = find_closest_hit(ray);

    glm::vec3 color(0.0);

    if (!closest_hit.hit)
        return color;

    color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv,
                       glm::normalize(-ray.direction), closest_hit.object->getMaterial());
    // cout << glm::to_string(color) << endl;

    Material material = closest_hit.object->getMaterial();
    glm::vec3 view_direction = glm::vec3(-ray.direction);
    glm::vec3 normal = closest_hit.normal;
    glm::vec3 point = closest_hit.intersection;

    if (material.refractiveness > 0) {
        float index =
            closest_hit.from_outside ? material.refractiveness : 1 / material.refractiveness;

        // glm::vec3 refraction_direction = glm::refract(ray.direction, closest_hit.normal, index);
        glm::vec3 refraction_direction = refract(-view_direction, normal, index);

        Ray refraction_ray(point + (EPS * refraction_direction), refraction_direction);
        glm::vec3 refractive_component = trace_ray(refraction_ray);
        color = (material.refractive_ratio) * color +
                (1.f - material.refractive_ratio) * refractive_component;
        // color = refractive_component;
    }

    if (material.reflectiveness > 0) {
        glm::vec3 relfective_component(0);
        glm::vec3 reflective_direction = glm::reflect(-view_direction, normal);
        Ray reflective_ray(point + (EPS * reflective_direction), reflective_direction);

        if (closest_hit.from_outside)
            relfective_component = trace_ray(reflective_ray);
        color = (1 - material.reflectiveness) * color +
                (material.reflectiveness * relfective_component);
    }
    // color *= 2;

    return color;
}

// Decode from disk to raw pixels with a single function call
PNG_Image_t * decodeOneStep(const char * filename) {
    std::vector<unsigned char> data; // the raw pixels
    unsigned width, height;

    // decode
    unsigned error = lodepng::decode(data, width, height, filename);

    // if there's an error, display it
    if (error) {
        std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << "(file: '"
                  << filename << "')" << std::endl;
        return 0;
    }

    PNG_Image_t * image = (PNG_Image_t *)malloc(sizeof(PNG_Image_t));
    assert(image);
    // printf("png image: %dx%d\n", width, height);
    image->width = width;
    image->height = height;
    image->data = std::vector<unsigned char>(data);

    // *image = {
    //     .width = width,
    //     .height = height,
    //     .data = data,
    // };

    // the pixels are now in the vector "image", 4 bytes per pixel, ordered RGBARGBA..., use it as
    // texture, draw it, ...
    return image;
}

glm::vec3 get_vertex(std::string & line) {
    std::istringstream stream(line);
    float p1, p2, p3;
    std::string v;
    stream >> v;
    stream >> p1;
    stream >> p2;
    stream >> p3;
    return glm::vec3(p1, p2, p3);
}

glm::vec2 get_vt(std::string line) {
    std::istringstream stream(line);
    float p1, p2;
    std::string vt;
    stream >> vt;
    stream >> p1;
    stream >> p2;
    return glm::vec2(p1, p2);
}

glm::vec3 get_normals(std::string & line) { return glm::normalize(get_vertex(line)); }

std::vector<std::string> split(std::string s, char c) {
    std::stringstream test(s);
    std::string segment;
    std::vector<std::string> seglist;

    while (std::getline(test, segment, c)) {
        seglist.push_back(segment);
    }
    return seglist;
}

Fragment * get_face(std::string line, std::vector<glm::vec3> vertices,
                    std::vector<glm::vec3> normals, std::vector<glm::vec2> vertex_textures) {
    std::istringstream stream(line);
    string s1, s2, s3, s4;
    std::string f;
    if (std::count(line.begin(), line.end(), ' ') == 3) {
        stream >> f;
        stream >> s1;
        stream >> s2;
        stream >> s3;
        int i1, i2, i3;
        i1 = atoi(s1.substr(0, s1.find('/')).c_str());
        i2 = atoi(s2.substr(0, s2.find('/')).c_str());
        i3 = atoi(s3.substr(0, s3.find('/')).c_str());
        return normals.size() > 0
                   ? new Triangle(vertices[i1 - 1], vertices[i2 - 1], vertices[i3 - 1],
                                  normals[i1 - 1], normals[i2 - 1], normals[i3 - 1])
                   : new Triangle(vertices[i1 - 1], vertices[i2 - 1], vertices[i3 - 1]);
    } else {
        stream >> f;
        stream >> s1;
        stream >> s2;
        stream >> s3;
        stream >> s4;
        int v1, v2, v3, v4;
        int vt1, vt2, vt3, vt4;
        int n1, n2, n3, n4;

        std::vector<std::string> s1_split = split(s1, '/');
        v1 = atoi(s1_split[0].c_str()) - 1;
        vt1 = atoi(s1_split[1].c_str()) - 1;
        n1 = atoi(s1_split[2].c_str()) - 1;

        std::vector<std::string> s2_split = split(s2, '/');
        v2 = atoi(s2_split[0].c_str()) - 1;
        vt2 = atoi(s2_split[1].c_str()) - 1;
        n2 = atoi(s2_split[2].c_str()) - 1;

        std::vector<std::string> s3_split = split(s3, '/');
        v3 = atoi(s3_split[0].c_str()) - 1;
        vt3 = atoi(s3_split[1].c_str()) - 1;
        n3 = atoi(s3_split[2].c_str()) - 1;

        std::vector<std::string> s4_split = split(s4, '/');
        v4 = atoi(s4_split[0].c_str()) - 1;
        vt4 = atoi(s4_split[1].c_str()) - 1;
        n4 = atoi(s4_split[2].c_str()) - 1;
        return normals.size() > 0
                   ? new Diamond(vertices[v1], vertices[v2], vertices[v3], vertices[v4],
                                 normals[n1], normals[n2], normals[n3], normals[n4],
                                 vertex_textures[vt1], vertex_textures[vt2], vertex_textures[vt3],
                                 vertex_textures[vt4])
                   : new Diamond(vertices[v1], vertices[v2], vertices[v3],
                                 vertices[v4]); // does not work
    }
}

Mesh * getMeshFromOBJ(std::string filename) {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec2> vertex_textures;
    std::vector<glm::vec3> normals;
    std::vector<Fragment *> faces;

    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '#')
            continue;

        if (line[0] == 'v')
            vertices.push_back(get_vertex(line));

        if (line.substr(0, 2) == "vt")
            vertex_textures.push_back(get_vt(line));

        if (line.substr(0, 2) == "vn")
            normals.push_back(get_normals(line));

        if (line[0] == 'f')
            faces.push_back(get_face(line, vertices, normals, vertex_textures));
    }

    return new Mesh(faces);
}

/**
 Function defining the scene
 */
void sceneDefinition() {

    Material green_diffuse;
    green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
    green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

    Material red_specular;
    red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
    red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
    red_specular.specular = glm::vec3(0.5);
    red_specular.shininess = 10.0;

    Material blue_specular;
    blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
    blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);
    blue_specular.specular = glm::vec3(0.6);
    blue_specular.shininess = 100.0;
    // blue_specular.reflectiveness = .9f;

    Material highly_specular_yellow;
    highly_specular_yellow.ambient = glm::vec3(0.1f, 0.1f, 0.03f);
    highly_specular_yellow.diffuse = glm::vec3(.6f, .6f, 0.1f);
    highly_specular_yellow.specular = glm::vec3(.6f);
    highly_specular_yellow.shininess = 100;

    glm::vec3 p1(0, 0, -1), p2(-1, 0, 0), p3(1, 0, 0), p4(0, 1, 0);
    glm::vec3 n1(cross(p3 - p1, p2 - p1)), n2(cross(p2 - p1, p4 - p1)), n3(cross(p4 - p1, p3 - p1)),
        n4(cross(p4 - p3, p2 - p3));
    n1 = glm::normalize(n1);
    n2 = glm::normalize(n2);
    n3 = glm::normalize(n3);
    n4 = glm::normalize(n4);

    Triangle *t1 = new Triangle(p1, p3, p2, n1, n1, n1, green_diffuse),
             *t2 = new Triangle(p2, p4, p1, n2, n2, n2, red_specular),
             *t3 = new Triangle(p1, p4, p3, n3, n3, n3, blue_specular),
             *t4 = new Triangle(p3, p4, p2, n4, n4, n4, highly_specular_yellow);

    std::vector<Fragment *> triangles;
    triangles.push_back(t1);
    triangles.push_back(t2);
    triangles.push_back(t3);
    triangles.push_back(t4);
    Mesh * pyramid = new Mesh(triangles);

    glm::mat4 pyr_translation = glm::translate(glm::vec3(0, -10, 15));
    glm::mat4 pyr_rotate = glm::rotate(.7f, glm::vec3(0, 1, 0));
    glm::mat4 pyr_scale = glm::scale(glm::vec3(1, 1, 1));
    pyramid->setTransformation(pyr_translation * pyr_rotate * pyr_scale);
    // objects.push_back(pyramid);

    Mesh * teapot = getMeshFromOBJ("teapot.obj");
    teapot->setMaterial(red_specular);
    glm::mat4 teapotTransformation = glm::translate(glm::vec3(0, 1, 10));
    teapot->setTransformation(teapotTransformation);
    Box * teapotBox = new Box(teapot->getMinCoords(), teapot->getMaxCoords(), teapot);
    // Box * teapotBox = new Box(glm::vec3(-3, 0, -2), glm::vec3(5, 5, 3), teapot);
    teapotBox->setMaterial(blue_specular);
    teapotBox->setTransformation(teapotTransformation);
    // objects.push_back(teapotBox);

    // Mesh * pyramid2 = new Mesh(triangles);
    // pyramid2->setTransformation(seashellTransformation);
    // objects.push_back(pyramid2);
    // pyramid->setTransformation(seashellTransformation);
    // objects.push_back(pyramid);

    Material shell_textured;
    // shell_textured.specular = glm::vec3(.9f);
    shell_textured.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
    shell_textured.diffuse = glm::vec3(0.8f, 0.2f, 0.2f);
    // shell_textured.diffuse = glm::vec3(.6f, 0.4f, 0.3f);
    shell_textured.specular = glm::vec3(0.5);
    shell_textured.shininess = 10.0;
    float angle = atan(-3);
    Mesh * seashell = getMeshFromOBJ("seashell_obj.obj");
    // Mesh * seashell = getMeshFromOBJ("shell.obj");

    seashell->setMaterial(shell_textured);

    glm::mat4 seashellTransformation = glm::translate(glm::vec3(-10, 5, 24.7)) *
                                       glm::rotate(angle, glm::vec3(1, 0, 0)) *
                                       glm::scale(glm::vec3(0.15));
                                       // glm::scale(glm::vec3(6));
    Box * seashellBox = new Box(seashell->getMinCoords(), seashell->getMaxCoords(), seashell);
    seashell->setTransformation(seashellTransformation);
    seashellBox->setTransformation(seashellTransformation);

    //
    // Mesh * seashell2 = getMeshFromOBJ("shell.obj");
    // Box * seashellBox2 = new Box(seashell2->getMinCoords(), seashell2->getMaxCoords(), seashell2);
    // glm::mat4 seashellTransformation2 = glm::translate(glm::vec3(10, -8, 16.7)) *
    //                          glm::rotate(angle, glm::vec3(1, -.5, 0)) *
    //                          glm::scale(glm::vec3(1));
    //                          // glm::scale(glm::vec3(6));
    // seashell2->setMaterial(shell_textured);
    // seashell2->setTransformation(seashellTransformation2);
    // seashellBox2->setTransformation(seashellTransformation2);

    Mesh * seashell3 = getMeshFromOBJ("shell.obj");
    Box * seashellBox3 = new Box(seashell3->getMinCoords(), seashell3->getMaxCoords(), seashell3);
    glm::mat4 seashellTransformation3 =  glm::translate(glm::vec3(9, -8, 16.7)) *
                             glm::rotate(angle, glm::vec3(1, .5, 0)) *
                             glm::scale(glm::vec3(1));
                             // glm::scale(glm::vec3(6));
    Material shell3_textured;
    shell3_textured.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
    shell3_textured.diffuse = glm::vec3(148,106,154)/255.f;
    // shell3_textured.diffuse = glm::vec3(.9f, 0.4f, 0.3f);
    shell3_textured.specular = glm::vec3(0.5);
    shell3_textured.shininess = 10.0;
    seashell3->setMaterial(shell3_textured);
    seashell3->setTransformation(seashellTransformation3);
    seashellBox3->setTransformation(seashellTransformation3);
    // pyramid->setTransformation(seashellTransformation3);


    Mesh * coral1 = getMeshFromOBJ("coral2.obj");
    Box * coralBox1 = new Box(coral1->getMinCoords(), coral1->getMaxCoords(), coral1);
    glm::mat4 coralTransformation1 = glm::translate(glm::vec3(-20, 9.5, 26.7)) *
                             glm::rotate((float) M_PI - .5f * angle, glm::vec3(1, 0, 0)) *
                             glm::scale(glm::vec3(0.1));
                             // glm::scale(glm::vec3(6));
    Material coral1_textured;
    coral1_textured.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
    // coral1_textured.diffuse = glm::vec3(0.7f, 0.2f, 0.2f);
    coral1_textured.diffuse = glm::vec3(.9f, 0.4f, 0.3f);
    coral1_textured.specular = glm::vec3(0.5);
    coral1_textured.shininess = 10.0;
    coral1->setMaterial(coral1_textured);
    coral1->setTransformation(coralTransformation1);
    coralBox1->setTransformation(coralTransformation1);
    // pyramid->setTransformation(seashellTransformation3);
    
    objects.push_back(seashellBox);
    // objects.push_back(seashellBox2);
    objects.push_back(seashellBox3);
    // objects.push_back(coralBox1);

    Material refractive;
    refractive.reflectiveness = 0.1f;
    refractive.refractiveness = 2.0f;

    Sphere * blue_sphere = new Sphere(blue_specular);
    blue_sphere->setTransformation(glm::translate(glm::vec3(1, -2, 8)));
    // objects.push_back(blue_sphere);
    Sphere * red_sphere = new Sphere(red_specular);
    red_sphere->setTransformation(
        glm::scale(glm::translate(glm::vec3(-1, -2.5, 6)), glm::vec3(.5)));
    // objects.push_back(red_sphere);

    Sphere * refractive_sphere = new Sphere(refractive);
    refractive_sphere->setTransformation(glm::translate(glm::vec3(-3, -1, 8)) *
                                         glm::scale(glm::vec3(2)));
    // objects.push_back(refractive_sphere);

    // Textured sphere
    Material stone_textured;
    stone_textured.specular = glm::vec3(.6f);
    stone_textured.shininess = 100;
    Material waffle_textured;
    waffle_textured.specular = glm::vec3(.6f);
    waffle_textured.shininess = 100;

    Material water_textured;
    water_textured.specular = glm::vec3(.6f);
    water_textured.refractiveness = .825f;
    water_textured.refractive_ratio = .2f;

    Material sand_textured;
    sand_textured.specular = glm::vec3(.2f);
    sand_textured.shininess = 100;
    // textured.texture = &rainbowTexture;
    Texture * sandImageTexture =
        new ImageTexture(*decodeOneStep("./textures/png/Sand_005_baseColor.png"),
                         *decodeOneStep("./textures/png/Sand_005_height.png"),
                         *decodeOneStep("./textures/png/Sand_005_normal.png"),
                         *decodeOneStep("./textures/png/Sand_005_ambientOcclusion.png"),
                         *decodeOneStep("./textures/png/Sand_005_roughness.png"));
    Texture * waterImageTexture =
        new ImageTexture(*decodeOneStep("./textures/png/Water_002_COLOR.png"),
                         *decodeOneStep("./textures/png/Water_002_DISP.png"),
                         *decodeOneStep("./textures/png/Water_002_NORM.png"),
                         *decodeOneStep("./textures/png/Water_002_OCC.png"),
                         *decodeOneStep("./textures/png/Water_002_ROUGH.png"));
    Texture * stoneImageTexture = new ImageTexture(
        *decodeOneStep("./textures/png/Stylized_Stone_Floor_005_basecolor.png"),
        *decodeOneStep("./textures/png/Stylized_Stone_Floor_005_height.png"),
        *decodeOneStep("./textures/png/Stylized_Stone_Floor_005_normal.png"),
        *decodeOneStep("./textures/png/Stylized_Stone_Floor_005_ambientOcclusion.png"),
        *decodeOneStep("./textures/png/Stylized_Stone_Floor_005_roughness.png"));
    Texture * waffleImageTexture =
        new ImageTexture(*decodeOneStep("./textures/png/Waffle_001_basecolor.png"),
                         *decodeOneStep("./textures/png/Waffle_001_height.png"),
                         *decodeOneStep("./textures/png/Waffle_001_normal.png"),
                         *decodeOneStep("./textures/png/Waffle_001_ambientOcclusion.png"),
                         *decodeOneStep("./textures/png/Waffle_001_roughness.png"));
    stone_textured.texture = stoneImageTexture;
    waffle_textured.texture = waffleImageTexture;
    water_textured.texture = waterImageTexture;
    sand_textured.texture = sandImageTexture;
    Sphere * waffle_sphere = new Sphere(waffle_textured);
    Sphere * water_sphere = new Sphere(water_textured);
    Sphere * stone_sphere = new Sphere(stone_textured);
    // rainbow_sphere->setTransformation(glm::rotate(glm::translate(glm::vec3(-6,4,23)), .2f,
    // glm::vec3(0, 1, 0)));
    glm::mat4 translation = glm::translate(glm::vec3(0, 1.5, 10));
    glm::mat4 rotation = glm::rotate(.7f, glm::vec3(1, 0, 0));
    glm::mat4 scaling = glm::scale(glm::vec3(1.5));
    waffle_sphere->setTransformation(translation * rotation * scaling);
    water_sphere->setTransformation(translation * rotation * scaling);
    // objects.push_back(water_sphere);

    translation = glm::translate(glm::vec3(4, 2, 10));
    rotation = glm::rotate(.2f, glm::vec3(0, 1, 0));
    scaling = glm::scale(glm::vec3(3));
    stone_sphere->setTransformation(translation * rotation * scaling);
    // objects.push_back(stone_sphere);

    // Planes
    Material red_diffuse;
    red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
    red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);

    Material blue_diffuse;
    blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
    blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
    objects.push_back(
        new Plane(glm::vec3(0, 0, 10), glm::normalize(glm::vec3(0.0, 1, -2)), water_textured));
    objects.push_back(
        new Plane(glm::vec3(0, 0, 25), glm::normalize(glm::vec3(0, 1, -3)), sand_textured));

    // objects.push_back(new Plane(glm::vec3(0, 1, 30), glm::vec3(0.0, 0.0, -1.0), stone_textured));
    // objects.push_back(new Plane(glm::vec3(-15, 1, 0), glm::vec3(1.0, 0.0, 0.0), red_diffuse));
    // objects.push_back(new Plane(glm::vec3(15, 1, 0), glm::vec3(-1.0, 0.0, 0.0), blue_diffuse));
    // objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0.0, -1, 0)));
    // objects.push_back(new Plane(glm::vec3(0, 1, -0.01), glm::vec3(0.0, 0.0, 1.0),
    // green_diffuse));

    glm::mat4 green_cone_trans = glm::scale(
        glm::translate(glm::vec3(6, -3, 7)) * glm::rotate((float)glm::atan(3), glm::vec3(0, 0, 1)),
        glm::vec3(1, 3, 1));

    Cone * cone1 = new Cone(green_diffuse);
    cone1->setTransformation(green_cone_trans);
    // objects.push_back(cone1);

    glm::mat4 yellow_cone_trans =
        glm::scale(glm::translate(glm::vec3(-5, 9, 14)) * glm::rotate(3.1415f, glm::vec3(0, 0, 1)),
                   glm::vec3(3, 12, 3));

    Cone * cone2 = new Cone(highly_specular_yellow);
    cone2->setTransformation(yellow_cone_trans);
    // objects.push_back(cone2);

    lights.push_back(new Light(glm::vec3(10, 16, -7), glm::vec3(.5f)));
    // lights.push_back(new Light(glm::vec3(2, 1, -15), glm::vec3(.5f)));
    // lights.push_back(new Light(glm::vec3(6, 1, 17), glm::vec3(0.3)));
    // lights.push_back(new Light(glm::vec3(2, 7, 1), glm::vec3(0.2)));
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity) {
    float gamma = 1.0 / 1.8;
    float alpha = 18.0f;
    // return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), 0.0f, 1.0f);
    return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), 0.0f, 1.0f);
}

long N_PIXELS;

void print_progress_bar(int pixel) {
    float progress = (float)pixel / N_PIXELS;
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    // std::cout << "] " << int(progress * 100.0) << " %\r";
    printf("] %.2f %%\r", progress * 100);
    std::cout.flush();
}

int main(int argc, const char * argv[]) {

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; // width of the image
    int height = 768; // height of the image
    float fov = 90;   // field of view
    N_PIXELS = width * height;

    cout << "Defining the scene..." << endl;
    sceneDefinition(); // Let's define a scene

    Image image(width, height); // Create an image where we will store the result
#ifdef ANTI_ALIASING
    cout << "Anti aliasing" << endl;
    // N_PIXELS *= 9;
    uint subpixel_weights[] = {1, 2, 1, 2, 4, 2, 1, 2, 1};
#endif

    float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
    float X = -s * width / 2;
    float Y = s * height / 2;

    cout << "Coloratin each pixel on the screne..." << endl;
    glm::vec3 origin(0, 0, 0);
    for (int i = 0; i < width; i++) {

        for (int j = 0; j < height; j++) {

            float dz = 1;
            // top left of the pixel
#ifdef ANTI_ALIASING
            float dx = X + i * s;
            float dy = Y - j * s;
            glm::vec3 color(0);
            for (int k = 0; k < 9; ++k) {
                float tmp_dx = dx + (k % 3) * s / 2;
                float tmp_dy = dy - ((int)(k / 3)) * s / 2;

                glm::vec3 direction(tmp_dx, tmp_dy, dz);
                // direction += glm::vec3(0, -1.5, 0);
                direction = glm::normalize(direction);

                Ray ray(origin, direction);
                color += (float)subpixel_weights[k] * trace_ray(ray);
                // print_progress_bar((i * height + j) * 9 + k);
            }

            color /= 16.f;
#else
            float dx = X + i * s + s / 2;
            float dy = Y - j * s + s / 2;
            glm::vec3 direction(dx, dy, dz);
            // direction += glm::vec3(0, -1.5, 0);
            Ray ray(origin, direction);
            glm::vec3 color = trace_ray(ray);
            // print_progress_bar(i * height + j);
#endif
            image.setPixel(i, j, toneMapping(color));
        }
        print_progress_bar((i + 1) * height);
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
