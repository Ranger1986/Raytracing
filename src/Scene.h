#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include <stdexcept>

#include <GL/glut.h>

const int MESHE = 0;
const int SPHERE = 1;
const int SQUARE = 2;

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

Vec3 getIntersection(RaySceneIntersection sceneInt){
    if (sceneInt.intersectionExists)
    {
        if (sceneInt.typeOfIntersectedObject == SPHERE) {
        return sceneInt.raySphereIntersection.intersection;
        }
        if (sceneInt.typeOfIntersectedObject == SQUARE) {
        return sceneInt.raySphereIntersection.intersection;
        }
    }
}

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
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        float distanceMinimum;
        float distanceObjet;
        result.intersectionExists=false;
        for (int i = 0; i < spheres.size(); i++)
        {
            RaySphereIntersection intersection = spheres[i].intersect(ray);
            if (intersection.intersectionExists)
            {
                distanceObjet = calculDistance(intersection.intersection, ray.origin());
                if (!result.intersectionExists)
                {
                    result.raySphereIntersection=intersection;
                    result.objectIndex=i;
                    result.typeOfIntersectedObject=SPHERE;
                    result.intersectionExists=true;
                    distanceMinimum=distanceObjet;
                }
                else if (distanceMinimum>distanceObjet)
                {
                    result.raySphereIntersection=intersection;
                    result.objectIndex=i;
                    result.typeOfIntersectedObject=SPHERE;
                    distanceMinimum=distanceObjet;
                }
            }
        }
        for (int i = 0; i < squares.size(); i++)
        {
            RaySquareIntersection intersection = squares[i].intersect(ray);
            if (intersection.intersectionExists)
            {
                distanceObjet = calculDistance(intersection.intersection, ray.origin());
                if (!result.intersectionExists)
                {
                    result.raySquareIntersection=intersection;
                    result.objectIndex=i;
                    result.typeOfIntersectedObject=SQUARE;
                    result.intersectionExists=true;
                    distanceMinimum=distanceObjet;
                }
                else if (distanceMinimum>distanceObjet)
                {
                    result.raySquareIntersection=intersection;
                    result.objectIndex=i;
                    result.typeOfIntersectedObject=SQUARE;
                    distanceMinimum=distanceObjet;
                }
            }
        }
        
        return result;
    }


    Vec3 phong (RaySceneIntersection intersectionObjet) {
        float Isa, Ka;
        float Isd, Kd;
        float Iss, Ks;
        Vec3 P, L, N, V, R;
        Vec3 color = Vec3(0., 0., 0.);
        Vec3 ambient = Vec3(0., 0., 0.);
        Vec3 sommeSpecular = Vec3(0,0,0);
        Vec3 sommeDifraction = Vec3(0,0,0);
        if (!intersectionObjet.intersectionExists) {
            return color;
        }
        for(int  iLight = 0; iLight < lights.size(); iLight++) {
            if (intersectionObjet.typeOfIntersectedObject == SPHERE) {
                RaySphereIntersection result_temp = intersectionObjet.raySphereIntersection;
                Sphere sph = spheres[intersectionObjet.objectIndex];

                P = result_temp.intersection;

                Vec3 L = lights[iLight].pos - P;
                L.normalize();

                Ray ray_lum = Ray(P, L);
                RaySceneIntersection RSI_lum = computeIntersection(ray_lum);

                result_temp.intersectionExists = RSI_lum.intersectionExists;

                /*
                Ray rayonOmbre= Ray(P, lights[iLight].pos-P );
                RaySceneIntersection intersectionOmbre= computeIntersection(rayonOmbre);
                
                if(intersectionOmbre.intersectionExists && (calculDistance(result_temp.intersection, lights[iLight].pos)>
                    calculDistance(getIntersection(intersectionOmbre), lights[iLight].pos))){result_temp.intersectionExists=false;}
                */
                if (result_temp.intersectionExists){
                    L = lights[iLight].pos - P;
                    L.normalize();

                    N = result_temp.normal;

                    V = -1. * (P);
                    V.normalize();

                    R = 2. * Vec3::dot(N, L) * N - L;
                    R.normalize();
                    
                    for(int i = 0; i < 3; i++) {

                        Isd = lights[iLight].material[i];
                        Kd = sph.material.diffuse_material[i];

                        Iss = lights[iLight].material[i];
                        Ks = sph.material.specular_material[i];
                        
                        sommeDifraction[i] += Isd * Kd * Vec3::dot(L, N);
                        sommeSpecular[i] += Iss * Ks * pow( fmax(0., Vec3::dot(R, V)), sph.material.shininess);
                    }
                }
            }

            else if (intersectionObjet.typeOfIntersectedObject == SQUARE) {
                RaySquareIntersection result_temp = intersectionObjet.raySquareIntersection;
                Square squ = squares[intersectionObjet.objectIndex];

                P = result_temp.intersection;
                
                Vec3 L = lights[iLight].pos - P;
                L.normalize();

                Ray rayonOmbre= Ray(P, L);
                RaySceneIntersection intersectionOmbre= computeIntersection(rayonOmbre);
                
                if(intersectionOmbre.intersectionExists/*&&
                    (calculDistance(result_temp.intersection, lights[iLight].pos)>
                    calculDistance(getIntersection(intersectionOmbre), lights[iLight].pos))*/)
                    {result_temp.intersectionExists=false;}

                if (result_temp.intersectionExists){
                    L = lights[iLight].pos - P;
                    L.normalize();

                    N = result_temp.normal;

                    V = -1. * (P);
                    V.normalize();

                    R = 2. * Vec3::dot(N, L) * N - L;
                    R.normalize();
                    
                    for(int i = 0; i < 3; i++) {

                        Isd = lights[iLight].material[i];
                        Kd = squ.material.diffuse_material[i];

                        Iss = lights[iLight].material[i];
                        Ks = squ.material.specular_material[i];
                        
                        sommeDifraction[i] += Isd * Kd * Vec3::dot(L, N);
                        sommeSpecular[i] += Iss * Ks * pow( fmax(0., Vec3::dot(R, V)), squ.material.shininess);
                    }
                }
            }
        }
        for(int i = 0; i < 3; i++){
            color[i]= ambient[i] + sommeDifraction[i] + sommeSpecular[i];
        }
        return color;
    }

    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces ) {
        //TODO 
        Vec3 color;
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        if (raySceneIntersection.intersectionExists)
        {
            color = phong(raySceneIntersection);
        }
        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
        //TODO appeler la fonction recursive
        
        Vec3 color = rayTraceRecursive(rayStart, 1);
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
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 0.01;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.0,0.0 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
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
            spheres.resize( spheres.size() + 2 );

            Sphere & s1 = spheres[spheres.size() - 1];
            s1.m_center = Vec3(1.0 , 0. , 0.);
            s1.m_radius = 1;
            s1.build_arrays();
            s1.material.type = Material_Mirror;
            s1.material.diffuse_material = Vec3( 1.,0.0,0.0 );
            s1.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s1.material.shininess = 20;

            Sphere & s2 = spheres[spheres.size() - 2];
            s2.m_center = Vec3(-1.0 , 0. , -20.);
            s2.m_radius = 1;
            s2.build_arrays();
            s2.material.type = Material_Mirror;
            s2.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s2.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s2.material.shininess = 20;
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
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_double_square() {
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
            light.material = Vec3(0,0,0);
            light.isInCamSpace = true;
        }
        /*
        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,1. );
            s.material.specular_material = Vec3( 1.,0.,1. );
            s.material.shininess = 16;
        }
        */
        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(1., 1., -0.5));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        {
            spheres.resize( spheres.size() + 1 );

            Sphere & s1 = spheres[spheres.size() - 1];
            s1.m_center = Vec3(1.0 , 0. , 0);
            s1.m_radius = 1;
            s1.build_arrays();
            s1.material.type = Material_Mirror;
            s1.material.diffuse_material = Vec3( 1.,0.0,0.0 );
            s1.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s1.material.shininess = 20;
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
        /*
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( -1.5, -1.5, 1 );
            light.radius = 1.f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(0.5,0,0);
            light.isInCamSpace = false;
        }
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0, 0, 2 );
            light.radius = 1.f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(0,0,1);
            light.isInCamSpace = false;
        }
        */
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
void setup_square_shadow(){
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
        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-0.25, -0.25, 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 0.5, 0.5);
            s.translate(Vec3(0., 0., -1.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }
    }

};



#endif
