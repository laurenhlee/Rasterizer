# Rasterizer
# Project Overview

This project is a **software rasterizer** implemented in C++ that renders 3D geometry using a configurable graphics pipeline.  
It supports **triangle rasterization**, **model/view/projection transformations**, **depth buffering**, **super-sampling anti-aliasing (SSAA)**, and **Blinn–Phong shading**.

The project was built on top of a provided graphics framework; my work focuses on the **core rasterization and shading logic**, implemented in `rasterizer_impl.cpp`.



## Features

- Triangle rasterization using edge tests
- Optional **Super Sampling Anti-Aliasing (SSAA)**
- Model, view, projection, and screen-space transformations
- Z-buffer–based hidden surface removal
- Per-pixel **Blinn–Phong shading**
  - Ambient, diffuse, and specular components
  - Multiple point lights
- YAML-based scene configuration



## Pipeline Overview

The rendering pipeline follows this order:

1. Load geometry from OBJ files
2. Apply **model transformations**
3. Apply **view (camera) transformation**
4. Apply **projection transformation**
5. Transform to **screen space**
6. Rasterize triangles
7. Perform **depth testing**
8. Shade visible fragments



## Input (YAML)

The renderer is driven by a YAML configuration file.  
The `task` field determines which stage(s) of the pipeline are active.

Supported tasks:
- `triangle`
- `transform`
- `shading`

### Common Fields

```yaml
task: triangle | transform | shading
antialias: none | SSAA
samples: <number of samples per pixel (ignored if antialias = none)>
resolution:
  width: <screen width in pixels>
  height: <screen height in pixels>
obj: <name of OBJ file, without .obj>
output: <name of output png file, wihtout .png>
```

#### Field descriptions:
* task: Selects which pipeline stages are enabled.
* antialias: Enables or disables super-sampling anti-aliasing.
* samples: Number of sub-pixel samples per pixel when SSAA is enabled.
* resolution: Output image resolution in pixels.
* obj: Name of the OBJ mesh file to render (located in objs/).
* output: Name of the output image file (written to outputs/).

### Task-Specific Fields
The following fields are only required for certain tasks. Rather than duplicating shared structures, they are grouped conceptually by purpose.
#### Transformation-Specific Fields (transform and shading tasks)
These fields control how geometry is positioned and viewed in the scene:
* Model transformation
    * Specifies translation, rotation, and scaling applied to the object in world space.
* Camera (view) parameters
    * Defines the camera position, look-at direction, and up vector.
* Projection parameters
    * Controls perspective projection settings such as field of view, near plane, and far plane.
Together, these fields define the full model–view–projection (MVP) transformation applied to each vertex.

#### Shading-Specific Fields (shading task only)
These fields control lighting and material appearance:
* Material properties
    * Ambient, diffuse, and specular reflection coefficients.
    * Shininess exponent for specular highlights.
* Light sources
    * One or more point lights defined by position and intensity.
* Ambient light
    * Global ambient illumination applied uniformly across the scene.
These parameters are used in a per-pixel Blinn–Phong shading model after rasterization and depth testing.

## Usage
```
./rasterizer <yaml file>
```
Example:
```
./rasterizer task-triangle.yaml
```

## Requirements
* C++17-compatible compiler (e.g., g++, clang++)
* GLM (OpenGL Mathematics) library
* Input yaml files and obj files

## Future Improvements
Potential extensions to the transformation and rendering pipeline include:
* Implement Multi-Sample Anti-Aliasing in the shading model
* Implement deferred shading
* Implement texture mapping
* Implement shadows
