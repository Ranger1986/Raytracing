#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }

    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersectionResult;
        intersectionResult.intersectionExists = false;
        Vec3 bottom_left=vertices[0].position;
        Vec3 rvector =vertices[1].position-bottom_left;
        Vec3 uvector = vertices[3].position-bottom_left;
        Vec3 normal = normal.cross(rvector,uvector);
        float denominator = Vec3::dot(normal, ray.direction());
        if (std::abs(denominator) != 0) {
            Vec3 toPlane = bottom_left - ray.origin();
            float t = Vec3::dot(toPlane, normal) / denominator;
            if (t >= 1e-6) {
                Vec3 intersectionPoint = ray.origin() + ray.direction() * t;
                Vec3 intersectionToVertex = intersectionPoint - bottom_left;
                float u = Vec3::dot(intersectionToVertex, rvector) / rvector.squareLength();
                float v = Vec3::dot(intersectionToVertex, uvector) / uvector.squareLength();
                if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
                    intersectionResult.intersectionExists = true;
                    intersectionResult.t = t;
                    intersectionResult.u = u;
                    intersectionResult.v = v;
                    intersectionResult.intersection = intersectionPoint;
                    intersectionResult.normal = normal;
                    intersectionResult.normal.normalize();
                }
            }
        }

        return intersectionResult;
}
};
#endif // SQUARE_H
