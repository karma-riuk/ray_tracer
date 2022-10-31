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
#define EPS 0.001f

using namespace std;

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
};

Hit find_closest_hit(Ray ray);
glm::vec3 trace_ray(Ray ray);

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
	
	void setTransformation(glm::mat4 matrix){
		transformationMatrix = matrix;
		inverseTransformationMatrix = glm::inverse(matrix);
		normalMatrix = glm::transpose(inverseTransformationMatrix);
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
    Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center){
		this->color = color;
    }
	Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center){
		this->material = material;
	}
	/** Implementation of the intersection function*/
    Hit intersect(Ray ray){

        glm::vec3 c = center - ray.origin;

        float cdotc = glm::dot(c,c);
        float cdotd = glm::dot(c, ray.direction);

        Hit hit;

        float D = 0;
		if (cdotc > cdotd*cdotd){
			D =  sqrt(cdotc - cdotd*cdotd);
		}
        if(D<=radius){
            hit.hit = true;
            float t1 = cdotd - sqrt(radius*radius - D*D);
            float t2 = cdotd + sqrt(radius*radius - D*D);

            float t = t1;
            if(t<0) t = t2;
            if(t<0){
                hit.hit = false;
                return hit;
            }

			hit.intersection = ray.origin + t * ray.direction;
			hit.normal = glm::normalize(hit.intersection - center);
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
			
			hit.uv.s = (asin(hit.normal.y) + M_PI/2)/M_PI;
			hit.uv.t = (atan2(hit.normal.z,hit.normal.x) + M_PI) / (2*M_PI);
        }
		else{
            hit.hit = false;
		}
		return hit;
    }
};


class Plane : public Object{

private:
	glm::vec3 normal;
	glm::vec3 point;

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
private:
	Plane *plane;
public:
	Cone(Material material){
		this->material = material;
		plane = new Plane(glm::vec3(0,1,0), glm::vec3(0.0,1,0));
	}
	Hit intersect(Ray ray){
		
		Hit hit;
		hit.hit = false;
		
		glm::vec3 d = inverseTransformationMatrix * glm::vec4(ray.direction, 0.0); //implicit cast to vec3
		glm::vec3 o = inverseTransformationMatrix * glm::vec4(ray.origin, 1.0); //implicit cast to vec3
		d = glm::normalize(d);
		
		
		float a = d.x*d.x + d.z*d.z - d.y*d.y;
		float b = 2 * (d.x * o.x + d.z * o.z - d.y * o.y);
		float c = o.x * o.x + o.z * o.z - o.y * o.y;
		
		float delta = b*b - 4 * a * c;
		
		if(delta < 0){
			return hit;
		}
		
		float t1 = (-b-sqrt(delta)) / (2*a);
		float t2 = (-b+sqrt(delta)) / (2*a);
		
		float t = t1;
		hit.intersection = o + t*d;
		if(t<0 || hit.intersection.y>1 || hit.intersection.y<0){
			t = t2;
			hit.intersection = o + t*d;
			if(t<0 || hit.intersection.y>1 || hit.intersection.y<0){
				return hit;
			}
		};
	
		hit.normal = glm::vec3(hit.intersection.x, -hit.intersection.y, hit.intersection.z);
		hit.normal = glm::normalize(hit.normal);
	
		
		Ray new_ray(o,d);
		Hit hit_plane = plane->intersect(new_ray);
		if(hit_plane.hit && hit_plane.distance < t && length(hit_plane.intersection - glm::vec3(0,1,0)) <= 1.0 ){
			hit.intersection = hit_plane.intersection;
			hit.normal = hit_plane.normal;
		}
		
		hit.hit = true;
		hit.object = this;
		hit.intersection = transformationMatrix * glm::vec4(hit.intersection, 1.0); //implicit cast to vec3
		hit.normal = (normalMatrix * glm::vec4(hit.normal, 0.0)); //implicit cast to vec3
		hit.normal = glm::normalize(hit.normal);
		hit.distance = glm::length(hit.intersection - ray.origin);
		
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
	for(int light_num = 0; light_num < lights.size(); light_num++){
		
		glm::vec3 light_position = lights[light_num]->position;
		glm::vec3 light_direction = glm::normalize(light_position - point);
		Ray light_ray(point + (EPS * light_direction), light_direction);
		Hit closest_hit = find_closest_hit(light_ray);
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
	
	
	color += refractive_component;
	color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
	glm::vec3 reflective_direction = glm::reflect(-view_direction, normal);
	Ray reflective_ray(point + (EPS * reflective_direction), reflective_direction);
	glm::vec3 relfective_component = material.reflectiveness > 0 ? trace_ray(reflective_ray) : glm::vec3(0);
	
	return (1 - material.reflectiveness) * color + (material.reflectiveness * relfective_component) ;
}

Hit find_closest_hit(Ray ray) {
	Hit closest_hit;

	closest_hit.hit = false;
	closest_hit.distance = INFINITY;

	for(int k = 0; k<objects.size(); k++){
		Hit hit = objects[k]->intersect(ray);
		if(hit.hit == true && hit.distance < closest_hit.distance)
			closest_hit = hit;
	}
	return closest_hit;
}


glm::vec3 refract(glm::vec3 i, glm::vec3 n, float index){
    glm::vec3 a = n * glm::dot(i, n);
    glm::vec3 b = i - a;
    float beta = 1/index;
    float alpha = sqrt(1 + (1 - pow(beta, 2)) * (glm::dot(b, b) / glm::dot(a, a)));
    return alpha * a + beta * b;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray){

	Hit closest_hit = find_closest_hit(ray);

	glm::vec3 color(0.0);

	if(closest_hit.hit){
		Material material = closest_hit.object->getMaterial();
		if (material.refractiveness > 0){
			bool from_outside = glm::dot(ray.direction, closest_hit.normal) <= 0 ;
			float index = from_outside ? material.refractiveness : 1 / material.refractiveness;
			cout << "ray direction " << glm::to_string(ray.direction) << endl;
			cout << "normal " << glm::to_string(closest_hit.normal) << endl;
			
			if (from_outside)
				cout << "from outside" << endl;
			else
				cout << "from inside" << endl;

			// glm::vec3 refraction_direction = glm::refract(ray.direction, closest_hit.normal, index);
			glm::vec3 refraction_direction = refract(ray.direction, closest_hit.normal, index);

			Ray refraction_ray(closest_hit.intersection + ( refraction_direction), refraction_direction);
            cout << "refraction direction " << glm::to_string(refraction_direction) << endl;
            color = trace_ray(refraction_ray);
		} else {
			// cout << "omg the end" << endl;
			color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
		}
	}else{
		color = glm::vec3(0.0, 0.0, 0.0);
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
	blue_specular.reflectiveness = 1.00f;
	blue_specular.shininess = 100.0;

	Material refractive;
	refractive.reflectiveness = .2f;
	refractive.refractiveness = 2.0f;


	objects.push_back(new Sphere(1.0, glm::vec3(1,-2,8), blue_specular));
	objects.push_back(new Sphere(0.5, glm::vec3(-1,-2.5,6), red_specular));
	objects.push_back(new Sphere(2.0, glm::vec3(-3,-1, 8), refractive));
	//objects.push_back(new Sphere(1.0, glm::vec3(3,-2,6), green_diffuse));
	
	
	//Textured sphere
	
	Material textured;
	textured.texture = &rainbowTexture;
	objects.push_back(new Sphere(7.0, glm::vec3(-6,4,23), textured));
	
	
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
	
	
	// Cone
	
	Material yellow_specular;
	yellow_specular.ambient = glm::vec3(0.1f, 0.10f, 0.0f);
	yellow_specular.diffuse = glm::vec3(0.4f, 0.4f, 0.0f);
	yellow_specular.specular = glm::vec3(1.0);
	yellow_specular.shininess = 100.0;
	
	Cone *cone = new Cone(yellow_specular);
	glm::mat4 translationMatrix = glm::translate(glm::vec3(5,9,14));
	glm::mat4 scalingMatrix = glm::scale(glm::vec3(3.0f, 12.0f, 3.0f));
	glm::mat4 rotationMatrix = glm::rotate(glm::radians(180.0f) , glm::vec3(1,0,0));
	cone->setTransformation(translationMatrix*scalingMatrix*rotationMatrix);
	objects.push_back(cone);
	
	Cone *cone2 = new Cone(green_diffuse);
	translationMatrix = glm::translate(glm::vec3(6,-3,7));
	scalingMatrix = glm::scale(glm::vec3(1.0f, 3.0f, 1.0f));
	rotationMatrix = glm::rotate(glm::atan(3.0f), glm::vec3(0,0,1));
	cone2->setTransformation(translationMatrix* rotationMatrix*scalingMatrix);
	objects.push_back(cone2);
	
	
	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
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
 cout << endl;
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
