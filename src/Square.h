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
        RaySquareIntersection intersection;
        //TODO calculer l'intersection rayon quand
        
        intersection.intersectionExists=false;
        float t;
        Vec3 bottomLeft=vertices[0].position;
        Vec3 rvect=vertices[1].position-bottomLeft;
        Vec3 uvect=vertices[3].position-bottomLeft;
        Vec3 normal = Vec3::cross(rvect , uvect);
        if (Vec3::dot(ray.direction(),normal)!=0){
            t=Vec3::dot((bottomLeft-ray.origin()),normal)/Vec3::dot(ray.direction(),normal);
        }
        if (t>0)
        {
            Vec3 supposed_intersection = ray.origin()+t*ray.direction();
            bool is_x_bound= ((vertices[0].position[0]<=supposed_intersection[0]&&supposed_intersection[0]<=vertices[1].position[0])||
            (vertices[0].position[0]<=supposed_intersection[0]&&supposed_intersection[0]<=vertices[2].position[0])||
            (vertices[0].position[0]<=supposed_intersection[0]&&supposed_intersection[0]<=vertices[3].position[0]))||
            (((vertices[0].position[0]>=supposed_intersection[0]&&supposed_intersection[0]>=vertices[1].position[0])||
            (vertices[0].position[0]>=supposed_intersection[0]&&supposed_intersection[0]>=vertices[2].position[0])||
            (vertices[0].position[0]>=supposed_intersection[0]&&supposed_intersection[0]>=vertices[3].position[0])));

            bool is_y_bound= (((vertices[0].position[1]<=supposed_intersection[1]&&supposed_intersection[1]<=vertices[1].position[1])||
            (vertices[0].position[1]<=supposed_intersection[1]&&supposed_intersection[1]<=vertices[2].position[1])||
            (vertices[0].position[1]<=supposed_intersection[1]&&supposed_intersection[1]<=vertices[3].position[1]))||
            ((vertices[0].position[1]>=supposed_intersection[1]&&supposed_intersection[1]>=vertices[1].position[1])||
            (vertices[0].position[1]>=supposed_intersection[1]&&supposed_intersection[1]>=vertices[2].position[1])||
            (vertices[0].position[1]>=supposed_intersection[1]&&supposed_intersection[1]>=vertices[3].position[1])));

            bool is_z_bound= (((vertices[0].position[2]<=supposed_intersection[2]&&supposed_intersection[2]<=vertices[1].position[2])||
            (vertices[0].position[2]<=supposed_intersection[2]&&supposed_intersection[2]<=vertices[2].position[2])||
            (vertices[0].position[2]<=supposed_intersection[2]&&supposed_intersection[2]<=vertices[3].position[2]))||
            ((vertices[0].position[2]>=supposed_intersection[2]&&supposed_intersection[2]>=vertices[1].position[2])||
            (vertices[0].position[2]>=supposed_intersection[2]&&supposed_intersection[2]>=vertices[2].position[2])||
            (vertices[0].position[2]>=supposed_intersection[2]&&supposed_intersection[2]>=vertices[3].position[2])));

            /*
            for (int i = 0; i < 3; i++)
            {
                float s_i=supposed_intersection[i];
                float v1=vertices[0].position[i];
                float v2=vertices[1].position[i];
                float v3=vertices[2].position[i];
                float v4=vertices[3].position[i];
                if ((supposed_intersection[i]<vertices[0].position[i])&&
                    (supposed_intersection[i]<vertices[1].position[i])&&
                    (supposed_intersection[i]<vertices[2].position[i])&&
                    (supposed_intersection[i]<vertices[3].position[i]))
                {
                    if (i==2)
                    {
                        std::cout<<"inf"<<std::endl;
                        std::cout<<supposed_intersection[i]<<std::endl;
                        std::cout<<vertices[0].position[i]<<std::endl;
                        std::cout<<vertices[1].position[i]<<std::endl;
                        std::cout<<vertices[2].position[i]<<std::endl;
                        std::cout<<vertices[3].position[i]<<std::endl;
                        std::cout<<(s_i<v1)<<std::endl;
                        std::cout<<(s_i<v2)<<std::endl;
                        std::cout<<(s_i<v3)<<std::endl;
                        std::cout<<(s_i<v4)<<std::endl;
                        std::cout<<false<<std::endl;

                    }
                    return intersection;
                }
                if ((supposed_intersection[i]>vertices[0].position[i])&&
                    (supposed_intersection[i]>vertices[1].position[i])&&
                    (supposed_intersection[i]>vertices[2].position[i])&&
                    (supposed_intersection[i]>vertices[3].position[i]))
                {
                    if (i==2)
                    {
                        std::cout<<"supp"<<std::endl;
                        std::cout<<supposed_intersection[i]<<std::endl;
                        std::cout<<vertices[0].position[i]<<std::endl;
                        std::cout<<vertices[1].position[i]<<std::endl;
                        std::cout<<vertices[2].position[i]<<std::endl;
                        std::cout<<vertices[3].position[i]<<std::endl;
                        std::cout<<((supposed_intersection[i]>vertices[0].position[i])&&
                        (supposed_intersection[i]>vertices[1].position[i])&&
                        (supposed_intersection[i]>vertices[2].position[i])&&
                        (supposed_intersection[i]>vertices[3].position[i]))<<std::endl;
                    }
                    return intersection;
                }
                */
                /*
                if ((vertices[0].position[i]<supposed_intersection[i]<vertices[2].position[i]))
                {
                    return intersection;
                }
                if ((vertices[0].position[i]<vertices[2].position[i])<supposed_intersection[i])
                {
                    return intersection;
                }
                */
                /*
                if (((vertices[2].position[i]<vertices[0].position[i]<supposed_intersection[i])||(vertices[0].position[i]<vertices[2].position[i]<supposed_intersection[i])))
                {
                    return intersection;
                }
                */
            if (is_x_bound&&is_y_bound&&is_z_bound)
            {
                intersection.intersectionExists=true;
                intersection.intersection=supposed_intersection;
            }
            
            
        }
        return intersection;
    }
};
#endif // SQUARE_H
