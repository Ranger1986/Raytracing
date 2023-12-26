#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>

#define MESHE 0
#define SPHERE 1
#define SQUARE 2

enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }




    RaySceneIntersection computeIntersection(Ray const & ray) {
        RaySceneIntersection result;
        float currentDistance = MAXFLOAT;
        float tempDistance = 0;
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        for (int i = 0; i < spheres.size(); i++)
        {
            RaySphereIntersection temp = spheres[i].intersect(ray);
            if (temp.intersectionExists)
            {  
                tempDistance = sqrt((temp.intersection-ray.origin()).norm());
                if (tempDistance < currentDistance)
                {
                    result.raySphereIntersection = temp;
                    result.typeOfIntersectedObject = SPHERE;
                    result.objectIndex=i;
                    result.intersectionExists=temp.intersectionExists;
                    currentDistance=tempDistance;
                }
            }
        }
        for (int i = 0; i < squares.size(); i++)
        {
            RaySquareIntersection temp = squares[i].intersect(ray);
            if (temp.intersectionExists)
            {  
                tempDistance = sqrt((temp.intersection-ray.origin()).norm());
                if (tempDistance < currentDistance)
                {
                    result.raySquareIntersection = temp;
                    result.typeOfIntersectedObject = SQUARE;
                    result.objectIndex=i;
                    result.intersectionExists=temp.intersectionExists;
                    currentDistance=tempDistance;
                }
            }
        }
        return result;
    }


    Vec3 Phong(Ray ray, RaySceneIntersection raySceneIntersection){
        Vec3 color;
        if (raySceneIntersection.intersectionExists)
        {   
            Material material;
            Vec3 normal;
            Vec3 intersection;

            switch (raySceneIntersection.typeOfIntersectedObject)
            {
                case SPHERE:
                    material = spheres[raySceneIntersection.objectIndex].material;
                    normal = raySceneIntersection.raySphereIntersection.normal;
                    intersection = raySceneIntersection.raySphereIntersection.intersection;
                    break;
                case SQUARE:
                    material = squares[raySceneIntersection.objectIndex].material;
                    normal = raySceneIntersection.raySquareIntersection.normal;
                    intersection = raySceneIntersection.raySquareIntersection.intersection;
                    break;
                default:
                    break;
            }

            Vec3 Ia, Id, Is;
            Ia = Vec3(0.1,0.1,0.1);
            Id = Vec3(0.f,0.f,0.f);
            Is = Vec3(0.f,0.f,0.f);

            for (int i = 0; i < lights.size(); i++)
            {
                Vec3 L = lights[i].pos-intersection;
                L.normalize();
                Vec3 N = normal;
                N.normalize();
                Vec3 R = 2 * Vec3::dot(N,L)*N-L;
                R.normalize();
                Vec3 V = ray.origin()-intersection;
                V.normalize();
                //std::cout << "i : " << i << std::endl;
                for (int j = 0; j < 3; j++)
                {
                    Id[j]+=lights[i].material[j]*material.diffuse_material[j]*fmax(0.,Vec3::dot(L,N));
                    Is[j]+=lights[i].material[j]*material.specular_material[j]*pow(fmax(0.,Vec3::dot(R,V)),50);
                    /*std::cout << "j : " << j << std::endl;
                    std::cout << Is[j] << std::endl;
                    std::cout << pow(fmax(0.,Vec3::dot(R,V)),material.shininess) << std::endl;
                    std::cout << Is[j] << std::endl;*/
                }
                
            }
            /*
            if (intersection[1] >=1-1e-2)
            {
                std::cout << "J'ai" << std::endl;
            }
            */
            
            color=Ia+Id+Is;
            //std::cout << Is[0] << "|" << Is[1] << "|" << Is[2] << std::endl;
        }
        return color;
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces ) {
        //TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 color;
        color=Phong(ray, raySceneIntersection);
        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
        //TODO appeler la fonction recursive
        Vec3 color=rayTraceRecursive(rayStart, 0);
        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(0,5,0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 20;
        }
    }

    void setup_double_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1 , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1 , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,0.0,0.0 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 200;
        }
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }
        /*
        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }
        */


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

};



#endif
