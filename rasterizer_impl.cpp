#include <cstdint>

#include "image.hpp"
#include "loader.hpp"
#include "rasterizer.hpp"


// helper functions

// returns the z value of the cross product (point - start) x (end - start)
inline double crossZ(const glm::vec2& point, const glm::vec2& start, const glm::vec2& end) {
    return (point.x - start.x)*(end.y - start.y) - (point.y - start.y)*(end.x - start.x);
}
// returns true if point is insdie the triangle given by v0, v1, and v2
inline bool insideTriangle(const glm::vec2& point, const glm::vec2& v0, const glm::vec2& v1, const glm::vec2& v2) {
    // calculate the z value of each cross product
    double z0 = crossZ(point, v0, v1);
    double z1 = crossZ(point, v1, v2);
    double z2 = crossZ(point, v2, v0);
    
   // compare their signs; if the point is on an edge, count it as inside
   return z0 >= 0 && z1 >= 0 && z2 >= 0 || z0 <= 0 && z1 <= 0 && z2 <= 0;
}



/** 
 * Given a single pixel in the screen space with the triangle in which the pixel is considered, determine the output color that should be rendered for the pixel.
 * This function will be called for every pixel in the bounding box of the triangle.
 * @param x: x coordinate of the pixel
 * @param y: y coordinate of the pixel
 * @param trig: the triangle in which the pixel is considered; see class `Triangle` in `entities.hpp`
 * @param config: the anti-aliasing configuration, which can be either `NONE` or `SSAA`
 * @param spp: the number of samples per pixel. Only useful if config is set to `SSAA`
 * @param image: the image to render the pixel on. See class `Image` in `image.hpp` for APIs of read/write operations
 * @param color: the color to render the pixel with, if the pixel is completely inside the triangle
*/
void Rasterizer::DrawPixel(uint32_t x, uint32_t y, Triangle trig, AntiAliasConfig config, uint32_t spp, Image& image, Color color)
{
    // define vertices of trig in 2D
        glm::vec2 v0(trig.pos[0].x, trig.pos[0].y);
        glm::vec2 v1(trig.pos[1].x, trig.pos[1].y);
        glm::vec2 v2(trig.pos[2].x, trig.pos[2].y);

    if (config == AntiAliasConfig::NONE)            // if anti-aliasing is off
    {
        // take one sample per pixel at (x + 0.5, y + 0.5)
        glm::vec2 check(x + 0.5, y + 0.5);

        // if point is inside the triangle, color the pixel
        if(insideTriangle(check, v0, v1, v2))
            image.Set(x, y, color); 
        
    }
    else if (config == AntiAliasConfig::SSAA)       // if anti-aliasing is on, use a grid algorithm with uniform distribution
    {
        // determine the smallest n such that we can take >=spp number of samples from an nxn grid
        uint32_t n = static_cast<int>(std::ceil(std::sqrt(spp)));
        double gridSize = 1.0 / n;

        double count = 0.0;
        for(int i = 0; i < n; i ++){
            for(int j = 0; j < n; j ++) {
                // get the x and y coordinates of the center of the (i,j)-th square
                double px = x + gridSize * (i + 0.5);
                double py = y + gridSize * (j + 0.5);
                if(insideTriangle(glm::vec2(px, py), v0, v1, v2))
                    count ++;
            }
        }

        // calculate the percentage of sampled points that were inside trig
        double inside = count / (n * n);

        // fill in the color scaled by the value of inside
        image.Set(x, y, color * inside);
    }
    return;
}
/**
 * Add the corresponding model transformation to the rasterizer. 
 * @param transform: the transformation to apply to the model. See class `MeshTransform` in `entities.hpp` for detailed definitions.
 * @param rotation: the rotation part of the transformation
*/
void Rasterizer::AddModel(MeshTransform transform, glm::mat4 rotation)
{
    // create translation matrix from transform.translation
    glm::mat4 Mtrans{
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        transform.translation.x, transform.translation.y, transform.translation.z, 1
    };

    // create scale matrix from transform.scale
    glm::mat4 Mscale{
        transform.scale.x, 0.0, 0.0, 0.0,
        0.0, transform.scale.y, 0.0, 0.0,
        0.0, 0.0, transform.scale.z, 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    // construct Mmodel by scaling first, then rotate, then translate
    glm::mat4 Mmodel = Mtrans * rotation * Mscale;

    model.push_back(Mmodel);
}

/**
 * Set the view transformation matrix. This function does not take any argument as `Rasterizer` contains the camera information in the `loader` member.
*/
void Rasterizer::SetView()
{
    // get camera from Loader
    const Camera& camera = this->loader.GetCamera();
    glm::vec3 cameraPos = camera.pos;
    glm::vec3 cameraLookAt = camera.lookAt;
    glm::vec3 cameraUp = camera.up;

    // construct orthogonal basis (f,h,-g) with the camera axes
    glm::vec3 g = glm::normalize(cameraLookAt - cameraPos); // g = lookAt vector
    glm::vec3 h = glm::normalize(cameraUp); // h = up vector
    glm::vec3 f = glm::normalize(glm::cross(g, h)); // f = binormal vector

    // create Mrot by orienting the camera towards the negative z-axis
    // this is the transpose of the transformation matrix from the canonical basis to (f,h,-g)
    glm::mat4 Mrot{
        f.x, h.x, -(g.x), 0.0,
        f.y, h.y, -(g.y), 0.0,
        f.z, h.z, (-g.z), 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    // create Mtrans by moving the camera to the origin
    glm::mat4 Mtrans{
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        -(cameraPos.x), -(cameraPos.y), -(cameraPos.z), 1
    };

    // multiply Mrot and Mtrans to get Mview
    // we translate the camera to the origin first, then rotate
    this->view = Mrot * Mtrans;
}

/**
* Set the projection transformation matrix. This function does not take any argument as `Rasterizer` contains the camera information in the `loader` member.
*/
void Rasterizer::SetProjection()
{
    // get camera parameters from Loader
    const Camera& camera = this->loader.GetCamera(); 

    float nearClip = camera.nearClip; // near clipping distance, strictly positive
    float farClip = camera.farClip; // far clipping distance, strictly positive
    
    float width = camera.width;
    float height = camera.height;
    
    // create Mpersp to squish the objects to conform to perspective
    // Mpersp: -n 0 0 0 
    //         0 -n 0 0
    //         a1 a2 a3 a4
    //         0 0 1 0       
    // where n = nearClip and a1 through a4 are constants

    // to find constants a1 and a2, we use the requirement that all points on the near clipping plane should be fixed.
    // so Mpersp * (x y -n 1)^T = (x y -n 1)^T for any x and y, meaning a1=a2=0 and a3 - 1/n(a4) = -n

    // to find constants a3 and a4, we use the requirement that the intersection on the z-axis and the far clipping plane should be fixed.
    // so Mpersp * (0 0 f 1)^T = (0 0 f 1) where f = farClip, and we have a linear system:
    // equation 1: a3 - 1/n(a4) = -n
    // equation 2: a3 - 1/f(a4) = -f
    // solving with n = nearClip and f = farClip:

    // a4 = (nearClip - farClip) / (1.0/nearClip - 1.0/farClip) = -(nearClip)*(farClip);
    // a3 = -(nearClip) + 1.0/nearClip * a4 = -(nearClip + farClip);

    float a4 = -(nearClip)*(farClip);
    float a3 = -(nearClip + farClip);

    glm::mat4 Mpersp{
        -(nearClip), 0.0, 0.0, 0.0,
        0.0, -(nearClip), 0.0, 0.0,
        0.0, 0.0, a3, 1.0,
        0.0, 0.0, a4, 0.0
    };

    // create Mortho to rescale and translate the objects to the canonical space
    
    glm::mat4 Mtrans{ // first translate the center to the origin
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 
        width/2.0, (height/2.0), (nearClip + farClip)/2.0, 1.0
    };

    glm::mat4 Mscale{ // then scale
        2.0/width, 0.0, 0.0, 0.0,
        0.0, 2.0/height, 0.0, 0.0, 
        0.0, 0.0, 2.0/(farClip-nearClip), 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    glm::mat4 Mortho = Mscale * Mtrans;

    // multiply Mortho and Mpersp to get Mprojection
    this->projection = Mortho * Mpersp;
}

/**
* Set the projection transformation matrix. This function does not take any argument as `Rasterizer` contains the camera information in the `loader` member.
*/
void Rasterizer::SetScreenSpace()
{
    // get screen parameters from Loader
    float width = this->loader.GetWidth();
    float height = this->loader.GetHeight();

    // create Mss based on the width and height of the image, and preserve the z value
    glm::mat4 Mss{
        width/2.0, 0.0, 0.0, 0.0,
        0.0, height/2.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    this->screenspace = Mss;
}


/**
 * Given a 2D position and a triangle, compute the barycentric coordinates of the position with respect to the triangle.
 * REQUIRE the `pos` fields of the triangle has `z` value 0 (i.e. it is guaranteed that the triangle is in same plane, and is inside the triangle)
 * @param pos: the 2D position to compute the barycentric coordinates for
 * @param trig: the triangle to compute the barycentric coordinates with respect to
 * @return: the barycentric coordinates of the position with respect to the triangle
 */

glm::vec3 Rasterizer::BarycentricCoordinate(glm::vec2 pos, Triangle trig)
{
    // solve for barycentric coordinates with the equations:
    // 1) (x,y) = alpha(A) + beta(B) + gamma(C)
    // 2) alpha + beta + gamma = 1
    glm::vec2 A(trig.pos[0].x, trig.pos[0].y);
    glm::vec2 B(trig.pos[1].x, trig.pos[1].y);
    glm::vec2 C(trig.pos[2].x, trig.pos[2].y);

    //float alpha = (-(pos.x - B.x)*(C.y - B.y) + (pos.y - B.y)*(C.x - B.x))/(-(A.x - B.x)*(C.y - B.y)+(A.y - B.y)*(C.x - B.x));
    //float beta = (-(pos.x - C.x)*(A.y - C.y) + (pos.y - C.y)*(A.x - C.x))/(-(B.x - C.x)*(A.y - C.y)+(B.y - C.y)*(A.x - C.x));
    //float gamma = 1 - alpha - beta;

    float area = crossZ(A, B, C);
    float alpha = crossZ(pos, B, C) / area;
    float beta = crossZ(pos, C, A) / area;
    float gamma = crossZ(pos, A, B) / area;

    return glm::vec3{alpha, beta, gamma};
}

// initialize zBuffer to the largest possible value
float Rasterizer::zBufferDefault = 1.0;

/**
 * Update the depth information at a single pixel in the ZBuffer. This function will be called for every pixel in the bounding box of the triangle.
 * @param x: x coordinate of the pixel
 * @param y: y coordinate of the pixel
 * @param original: the original triangle in the model space (before MVP transformation)
 * @param transformed: the transformed triangle in the screen space (after MVP transformation)
 * @param ZBuffer: the ZBuffer to update the depth information in. See class `Image` in `image.hpp` for APIs of read/write operations
 */
void Rasterizer::UpdateDepthAtPixel(uint32_t x, uint32_t y, Triangle original, Triangle transformed, ImageGrey& ZBuffer)
{
    // compute the barycentric coordinates of the pixel
    glm::vec3 bc = BarycentricCoordinate(glm::vec2{x, y},transformed);

    // compute the depth of the pixel using its barycentric coordinates
    float depth = bc[0]*transformed.pos[0].z + bc[1]*transformed.pos[1].z + bc[2]*transformed.pos[2].z;
     
    // if depth is smaller than ZBuffer then updte ZBuffer
    auto currentZBuffer = ZBuffer.Get(x,y);
    if(currentZBuffer.has_value()) {
        float ZB = currentZBuffer.value();
        if(depth < ZB) {
            ZBuffer.Set(x, y, depth);
        }
    }
}

/**
 * Shade the pixel at the given position, using Blinn-Phong shading model. This function will be called for every pixel in the bounding box of the triangle.
 * @param x: x coordinate of the pixel
 * @param y: y coordinate of the pixel
 * @param original: the original triangle in the model space (before MVP transformation)
 * @param transformed: the transformed triangle in the screen space (after MVP transformation)
 * @param image: the image to render the pixel on. See class `Image` in `image.hpp` for APIs of read/write operations
 */
void Rasterizer::ShadeAtPixel(uint32_t x, uint32_t y, Triangle original, Triangle transformed, Image& image)
{
    // compute the depth of the pixel using barycentric coordinates
    glm::vec3 bc = BarycentricCoordinate(glm::vec2{x, y},transformed);
    float depth = bc[0]*transformed.pos[0].z + bc[1]*transformed.pos[1].z + bc[2]*transformed.pos[2].z;
    
    auto currentZBuffer = ZBuffer.Get(x,y);
    if(currentZBuffer.has_value()) {
        float ZB = currentZBuffer.value();

        const float nearZero = 1e-5f;
        if (std::abs(depth - ZB) < nearZero) { // if depth is exactly the value in ZBuffer
            // compute the position corresponding to the pixel in the world space
            glm::vec3 pixel = bc[0] * glm::vec3(original.pos[0]) + bc[1] * glm::vec3(original.pos[1]) + bc[2] * glm::vec3(original.pos[2]);

            // retrieve its normal using barycentric coordinates
            glm::vec3 normal = bc[0] * glm::vec3(original.normal[0]) + bc[1] * glm::vec3(original.normal[1]) + bc[2] * glm::vec3(original.normal[2]);
            normal = glm::normalize(normal);

            // compute the color based on the normal and light directions, using the Blinn-Phong model
            // first load the camera, lights, exponent, and ambient from the Loader
            const Camera& camera = this->loader.GetCamera(); 
            const std::vector<Light>& lights = this->loader.GetLights(); 
            float exponent = this->loader.GetSpecularExponent();
            Color ambient = this->loader.GetAmbientColor();

            Color result;
            // 1. ambient component
            result = ambient;

            // iterate through all lights
            for(Light light : lights) {
                // compute squared between object and light
                float r2 = (light.pos.x - pixel.x)*(light.pos.x - pixel.x)
                         + (light.pos.y - pixel.y)*(light.pos.y - pixel.y)
                         + (light.pos.z - pixel.z)*(light.pos.z - pixel.z);

                // compute light direction
                glm::vec3 l = glm::normalize(light.pos - pixel);
                
                // 2. add the diffuse (lambertian) component
                float dot = glm::dot(l,normal);
                if(dot > 0) {
                    result.r += light.color.r * (light.intensity * 1.0 / r2 * dot);
                    result.g += light.color.g * (light.intensity * 1.0 / r2 * dot);
                    result.b += light.color.b * (light.intensity * 1.0 / r2 * dot);
                }
                    
                // compute direction from object to camera
                glm::vec3 v = glm::normalize(camera.pos - pixel);
                // compute h, the half vector of l and v
                glm::vec3 h{
                    (l.x + v.x) / 2.0,
                    (l.y + v.y) / 2.0,
                    (l.z + v.z) / 2.0
                };
                //glm::vec3 h = glm::normalize(l + v);

                // 3. add the specular component
                dot = glm::dot(normal, h);
                if(dot > 0) {
                    result.r += light.color.r * (light.intensity * 1.0/r2 * std::pow(dot, exponent));
                    result.g += light.color.g * (light.intensity * 1.0/r2 * std::pow(dot, exponent));
                    result.b += light.color.b * (light.intensity * 1.0/r2 * std::pow(dot, exponent));
                }
            }

            if(insideTriangle(glm::vec2{x, y}, glm::vec2(transformed.pos[0]), glm::vec2(transformed.pos[1]), glm::vec2(transformed.pos[2])))
                image.Set(x, y, result);
            
        }
    }
    
}

