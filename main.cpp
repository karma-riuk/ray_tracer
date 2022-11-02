/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"
#include "glm/gtx/string_cast.hpp"

#include "Image.h"
#include "Material.h"

using namespace std;
#define EPS 0.0001f

/**
 Class representing a single ray.
 */
class Ray{
public:
    glm::vec3 origin; ///< Origin of the ray
    glm::vec3 direction; ///< Direction of the ray
	/**
	 Contructor of the ray
	 @param origin Origin of the ray
	 @param direction Direction of the ray
	 */
    Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction){
    }
};


class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit{
    bool hit; ///< Boolean indicating whether there was or there was no intersection with an object
    glm::vec3 normal; ///< Normal vector of the intersected object at the intersection point
    glm::vec3 intersection; ///< Point of Intersection
    float distance; ///< Distance from the origin of the ray to the intersection point
    Object *object; ///< A pointer to the intersected object
	glm::vec2 uv; ///< Coordinates for computing the texture (texture coordinates)
    bool from_outside = true;
};

/**
 General class for the object
 */
class Object{
	
protected:
	glm::mat4 transformationMatrix; ///< Matrix representing the transformation from the local to the global coordinate system
	glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
	glm::mat4 normalMatrix; ///< Matrix for transforming normal vectors from the local to the global coordinate system
	
public:
	glm::vec3 color; ///< Color of the object
	Material material; ///< Structure describing the material of the object
	/** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;

	/** Function that returns the material struct of the object*/
	Material getMaterial(){
		return material;
	}
	/** Function that set the material
	 @param material A structure describing the material of the object
	*/
	void setMaterial(Material material){
		this->material = material;
	}
	/** Functions for setting up all the transformation matrices
	 @param matrix The matrix representing the transformation of the object in the global coordinates */
	void setTransformation(glm::mat4 matrix){
		
		transformationMatrix = matrix;
		
		inverseTransformationMatrix = glm::inverse(matrix);
        normalMatrix = (glm::abs(glm::determinant(matrix)) == 1) 
            ? matrix 
            : glm::transpose(inverseTransformationMatrix);
	}
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object{
private:
    float radius; ///< Radius of the sphere
    glm::vec3 center; ///< Center of the sphere

public:
	/**
	 The constructor of the sphere
	 @param radius Radius of the sphere
	 @param center Center of the sphere
	 @param color Color of the sphere
	 */
    Sphere(glm::vec3 color) : radius(1), center(glm::vec3(0)){
		this->color = color;
    }
	Sphere(Material material) : radius(1), center(glm::vec3(0)){
		this->material = material;
	}
	/** Implementation of the intersection function*/
    Hit intersect(Ray ray){

        glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1);
        glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0);

        Ray local_ray(o, d);
        Hit hit;
        hit.hit = false;

        float a, b, c;
        a = pow(d.x, 2) + pow(d.y, 2) + pow(d.z, 2);
        b = 2 * (d.x * o.x + d.y * o.y + d.z * o.z);
        c = pow(o.x, 2) + pow(o.y, 2) + pow(o.z, 2) - 1;

        float delta  = pow(b, 2) - 4 * a * c;
        if ( delta < 0 )
            return hit;

        float t_1, t_2;
        t_1 = (-b + sqrt(delta)) / ( 2 * a );
        t_2 = (-b - sqrt(delta)) / ( 2 * a );


        if (t_1 < 0 && t_2 < 0)
            return hit;
        glm::vec3 i1, i2;
        i1 = o + t_1 * d;
        i2 = o + t_2 * d;

        float t;
        t = min(t_1 < 0 ? INFINITY: t_1, t_2 < 0 ? INFINITY: t_2); // find lowest positive number 
                                                                   //
        glm::vec3 i = o + t * d;
        hit.intersection = i;
        hit.normal = i;
        hit.uv.x = (asin(i.y) + M_PI/2)/M_PI;
        hit.uv.y = (atan2(i.z, i.x) + M_PI)/(2*M_PI);

        bool from_inside = glm::all(glm::lessThan(glm::abs(o), glm::vec3(1)));
        hit.from_outside = !from_inside;

        // Retransform in global coordinates
        hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1);
        hit.normal = normalMatrix * glm::vec4(hit.normal, 0);
        hit.normal = glm::normalize(hit.normal);
		hit.distance = glm::distance(ray.origin, hit.intersection);

        hit.object = this;
        hit.hit = true;

        return hit;
    }
};


class Plane : public Object{

private:
	glm::vec3 point;
	glm::vec3 normal;

public:
	Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal){
	}
	Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal){
		this->material = material;
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;
		float DdotN = glm::dot(ray.direction, normal);
		if(DdotN < 0){
			
			float PdotN = glm::dot (point-ray.origin, normal);
			float t = PdotN/DdotN;
			
			if(t > 0){
				hit.hit = true;
				hit.normal = normal;
				hit.distance = t;
				hit.object = this;
				hit.intersection = t * ray.direction + ray.origin;
			}
		}
		return hit;
	}
};

class Cone : public Object{
    Plane base;
public:
	Cone(Material material): base(glm::vec3(0, 1, 0), glm::vec3(0, 1, 0), material){
		this->material = material;
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;

        glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0);
        glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1);
        Ray local_ray(o, d);


        float a, b, c; // constant to solve at^2 + bt + c = 0
        a = pow(d.x, 2) - pow(d.y, 2) + pow(d.z, 2);
        b = 2 * (d.x * o.x - d.y * o.y + d.z * o.z);
        c = pow(o.x, 2) - pow(o.y, 2) + pow(o.z, 2);

        float delta  = pow(b, 2) - 4 * a * c;
        if ( delta < 0 )
            return hit;

        float t_1, t_2;
        t_1 = (-b + sqrt(delta)) / ( 2 * a );
        t_2 = (-b - sqrt(delta)) / ( 2 * a );


        if (t_1 < 0 && t_2 < 0)
            return hit;
        glm::vec3 i1, i2;
        i1 = o + t_1 * d;
        i2 = o + t_2 * d;
        if ((i1.y < 0 && i2.y < 0) || (i1.y > 1 && i2.y > 1))
            return hit;

        float t;
        t = min(t_1 < 0 ? INFINITY: t_1, t_2 < 0 ? INFINITY: t_2); // find lowest positive number 
		 
        glm::vec3 i;
        i = o + t * d;
        if (i.y > 1) {
            hit = base.intersect(local_ray);
        }
        else{
            hit.hit = true;
            hit.object = this;
            hit.intersection = i;
            float y = (pow(i.x, 2) + pow(i.z, 2))/i.y + i.y;

            glm::vec3 normal = i - glm::vec3(0, y, 0) ;
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


/**
 Light class
 */
class Light{
public:
	glm::vec3 position; ///< Position of the light source
	glm::vec3 color; ///< Color/intentisty of the light source
	Light(glm::vec3 position): position(position){
		color = glm::vec3(1.0);
	}
	Light(glm::vec3 position, glm::vec3 color): position(position), color(color){
	}
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(0.001,0.001,0.001);
vector<Object *> objects; ///< A list of all objects in the scene
glm::vec3 tray_ray(Ray ray);
Hit find_closest_hit(Ray ray);

glm::vec3 refract(glm::vec3 i, glm::vec3 n, float index){
    glm::vec3 a = n * glm::dot(i, n);
    glm::vec3 b = i - a;
    float beta = 1/index;
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
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction, Material material){
    glm::vec3 refractive_component(0);

    glm::vec3 color = ambient_light * material.ambient;

    for(size_t light_num = 0; light_num < lights.size(); light_num++){

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
        if(material.texture){
            diffuse_color = material.texture(uv);
        }

        glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
        glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));


        // distance to the light
        float r = glm::distance(point,lights[light_num]->position);
        r = max(r, 0.1f);

        color += lights[light_num]->color * (diffuse + specular) / r/r;
    }

    color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));

    return color;
}

Hit find_closest_hit(Ray ray) {
	Hit closest_hit;

	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	for(size_t k = 0; k<objects.size(); k++){
		Hit hit = objects[k]->intersect(ray);
		if(hit.hit == true && hit.distance > 0.01f && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}
	return closest_hit;
}


/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray){

    Hit closest_hit = find_closest_hit(ray);

    glm::vec3 color(0.0);

    if(!closest_hit.hit)
        return color;

    color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv, glm::normalize(-ray.direction), closest_hit.object->getMaterial());

    Material material = closest_hit.object->getMaterial();
    glm::vec3 view_direction = glm::vec3(-ray.direction);
    glm::vec3 normal = closest_hit.normal;
    glm::vec3 point = closest_hit.intersection;


    if (material.refractiveness > 0){
        float index = closest_hit.from_outside ? material.refractiveness : 1 / material.refractiveness;

        // glm::vec3 refraction_direction = glm::refract(ray.direction, closest_hit.normal, index);
        glm::vec3 refraction_direction = refract(-view_direction, normal, index);


        Ray refraction_ray(point + (EPS * refraction_direction), refraction_direction);
        color = trace_ray(refraction_ray);
    }

    if (material.reflectiveness > 0){
        glm::vec3 relfective_component(0);
        glm::vec3 reflective_direction = glm::reflect(-view_direction, normal);
        Ray reflective_ray(point + (EPS * reflective_direction), reflective_direction);


        if (closest_hit.from_outside)
            relfective_component = trace_ray(reflective_ray);
        color = (1 - material.reflectiveness) * color + (material.reflectiveness * relfective_component);
    }

    return color;
}
/**
 Function defining the scene
 */
void sceneDefinition (){

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
	green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
	red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
	blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;
    blue_specular.reflectiveness = .9f;

    Material refractive{
        .reflectiveness = 0.1f,
        .refractiveness = 2.0f
    };



    Sphere * blue_sphere = new Sphere(blue_specular);
    blue_sphere->setTransformation(glm::translate(glm::vec3(1,-2,8)));
	objects.push_back(blue_sphere);
    Sphere * red_sphere = new Sphere(red_specular);
    red_sphere->setTransformation(glm::scale(glm::translate(glm::vec3(-1,-2.5,6)), glm::vec3(.5)));
	objects.push_back(red_sphere);

    Sphere * refractive_sphere = new Sphere(refractive);
    refractive_sphere->setTransformation(glm::translate(glm::vec3(-3, -1, 8)) * glm::scale(glm::vec3(2)));
	objects.push_back(refractive_sphere);
	
	
	
	
	
	//Textured sphere
	Material textured;
	textured.texture = &rainbowTexture;
    Sphere * rainbow_sphere = new Sphere(textured);
    // rainbow_sphere->setTransformation(glm::rotate(glm::translate(glm::vec3(-6,4,23)), .2f, glm::vec3(0, 1, 0)));
    glm::mat4 translation = glm::translate(glm::vec3(-6,4,23));
    glm::mat4 rotation = glm::rotate(.2f, glm::vec3(0, 1, 0));
    glm::mat4 scaling = glm::scale(glm::vec3(7));
    rainbow_sphere->setTransformation(translation * rotation * scaling);
	objects.push_back(rainbow_sphere);
	
	
	//Planes
	Material red_diffuse;
	red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
	red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);
		
	Material blue_diffuse;
	blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
	blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
	objects.push_back(new Plane(glm::vec3(0,-3,0), glm::vec3(0.0,1,0)));
	objects.push_back(new Plane(glm::vec3(0,1,30), glm::vec3(0.0,0.0,-1.0), green_diffuse));
	objects.push_back(new Plane(glm::vec3(-15,1,0), glm::vec3(1.0,0.0,0.0), red_diffuse));
	objects.push_back(new Plane(glm::vec3(15,1,0), glm::vec3(-1.0,0.0,0.0), blue_diffuse));
	objects.push_back(new Plane(glm::vec3(0,27,0), glm::vec3(0.0,-1,0)));
	objects.push_back(new Plane(glm::vec3(0,1,-0.01), glm::vec3(0.0,0.0,1.0), green_diffuse));
	

    glm::mat4 green_cone_trans = 
        glm::scale(
        glm::translate(glm::vec3(6, -3, 7)) *
        glm::rotate((float) glm::atan(3), glm::vec3(0, 0, 1)), glm::vec3(1, 3, 1));
	
	
	Cone *cone1 = new Cone(green_diffuse);
	cone1->setTransformation(green_cone_trans);
	objects.push_back(cone1);

    Material highly_specular_yellow{
        .ambient = glm::vec3(0.1f, 0.1f, 0.03f),
        .diffuse = glm::vec3(.6f, .6f, 0.1f),
        .specular = glm::vec3(.6f),
        .shininess = 100,
    };

    glm::mat4 yellow_cone_trans = 
        glm::scale(
        glm::translate(glm::vec3(5, 9, 14)) *
        glm::rotate(3.1415f, glm::vec3(0, 0, 1)), glm::vec3(3, 12, 3));
	
	Cone *cone2 = new Cone(highly_specular_yellow);
	cone2->setTransformation(yellow_cone_trans);
	objects.push_back(cone2);
	
	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(.4f)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.2)));
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity){
	float gamma = 1.0/2.0;
	float alpha = 12.0f;
	return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char * argv[]) {

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; //width of the image
    int height = 768; // height of the image
    float fov = 90; // field of view

	sceneDefinition(); // Let's define a scene

	Image image(width,height); // Create an image where we will store the result

    float s = 2*tan(0.5*fov/180*M_PI)/width;
    float X = -s * width / 2;
    float Y = s * height / 2;

    for(int i = 0; i < width ; i++)
        for(int j = 0; j < height ; j++){

			float dx = X + i*s + s/2;
            float dy = Y - j*s - s/2;
            float dz = 1;

			glm::vec3 origin(0, 0, 0);
            glm::vec3 direction(dx, dy, dz);
            direction = glm::normalize(direction);

            Ray ray(origin, direction);

			image.setPixel(i, j, toneMapping(trace_ray(ray)));

        }

    t = clock() - t;
    cout<<"It took " << ((float)t)/CLOCKS_PER_SEC<< " seconds to render the image."<< endl;
    cout<<"I could render at "<< (float)CLOCKS_PER_SEC/((float)t) << " frames per second."<<endl;

	// Writing the final results of the rendering
	if (argc == 2){
		image.writeImage(argv[1]);
	}else{
		image.writeImage("./result.ppm");
	}

    return 0;
}
