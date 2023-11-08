#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>
#include "aide.h"

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }


    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        //TODO calcul l'intersection rayon sphere
        float rayon = this->m_radius;
        Vec3 centre = this->m_center;
        Vec3 direction = ray.direction();
        Vec3 origin =ray.origin();
        float a = Vec3::dot(direction, direction);
        float b = 2.0f * Vec3::dot(direction,origin-centre);
        float c = Vec3::dot(origin-centre, origin-centre)-pow(rayon,2);
        float discriminant =pow(b,2)-4*a*c;
        if (discriminant<0)
        {
            intersection.intersectionExists = false;
        }
        if (discriminant==0)
        {
            intersection.intersectionExists=true;
            intersection.intersection=direction*(-b/(2*a));
        }
        if (discriminant>0)
        {
            intersection.intersectionExists=true;
            Vec3 intersection1=ray.origin() + direction*((-b+sqrt(discriminant))/(2*a));
            Vec3 intersection2=ray.origin() +direction*((-b-sqrt(discriminant))/(2*a));
            if (calculDistance(ray.origin(), intersection1)<calculDistance(ray.origin(), intersection2))
            {
                intersection.intersection=intersection1;
                intersection.secondintersection=intersection2;
            }
            else
            {
                intersection.intersection=intersection2;
                intersection.secondintersection=intersection1;
            }
        }
        /*
        if (intersection.intersectionExists){
            std::cout << "intersection:" << std::endl;
            printVec(intersection.intersection);
            std::cout << "origine:" << std::endl;
            printVec(origin);
            std::cout << "direction:" << std::endl;
            printVec(direction);
            std::cout << "racine:" << std::endl;
            std::cout << ((-b+sqrt(discriminant))/(2*a)) << std::endl;
            std::cout << ((-b-sqrt(discriminant))/(2*a)) << std::endl;
            std::cout << "distance:" << std::endl;
            std::cout << calculDistance(intersection.intersection,origin) << std::endl;
        }
        */
        return intersection;
    }
};
#endif
