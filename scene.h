#pragma once

#include "lajolla.h"
#include "camera.h"
#include "shape.h"

#include <memory>
#include <vector>

struct Scene {
    Scene(const RTCDevice &embree_device,
          const Camera &camera,
          const std::vector<std::shared_ptr<const Shape>> &shapes);
    ~Scene();

    RTCDevice embree_device;
    RTCScene embree_scene;
    // We decide to maintain a copy of the scene here.
    // This allows us to manage the memory of the scene ourselves and decouple
    // from the scene parser, but it's obviously less efficient.
    Camera camera;
    // For now we use a stl vector to store a list of pointers to our shapes & materials.
    // This wouldn't work if we want to extend this to run on GPUs.
    // If we want to port this to GPUs later, we need to flatten the array.
    std::vector<std::shared_ptr<const Shape>> shapes;
};