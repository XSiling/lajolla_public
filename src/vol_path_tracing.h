#pragma once

#include "scene.h"
#include "pcg.h"

int update_medium(const PathVertex& vertex, const Ray& ray, int medium_id) {
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        if (dot(ray.dir, vertex.geometric_normal)>0) {
            medium_id = vertex.exterior_medium_id;
        }
        else {
            medium_id = vertex.interior_medium_id;
        }
    }
    return medium_id;
}



// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    PathVertex vertex = *vertex_;

    if (!vertex_) {
        return make_zero_spectrum();
    }
    else {
        Vector3 p = vertex.position;
        Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], ray.org);
        Real t = distance(ray.org, p);
        Spectrum transmittance = exp(-sigma_a * t);
        Spectrum Le = Vector3(0, 0, 0);
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return transmittance * Le;
    }

}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    PathVertex vertex = *vertex_;

    // if no hit, set to infinite far
    Real t_hit = 0;
    if (!vertex_) {
        t_hit = infinity<Real>();
    }
    else {
        t_hit = distance(ray.org, vertex.position);
    }

    Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], ray.org);
    Spectrum sigma_s = get_sigma_s(scene.media[scene.camera.medium_id], ray.org);
    Spectrum sigma_t = sigma_a + sigma_s;
    PhaseFunction phase_function = get_phase_function(scene.media[scene.camera.medium_id]);

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t[0]; //since mono, from vector to real.

    if (t < t_hit){
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Spectrum p = ray.org + t * ray.dir;
        
        int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
        Real light_pdf = light_pmf(scene, light_id);

        Vector2 rnd_param_uv = Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
        Real rnd_param_w = next_pcg32_real<Real>(rng);
        PointAndNormal light_point = sample_point_on_light(scene.lights[light_id], p, rnd_param_uv, rnd_param_w, scene);
        Real light_point_pdf = pdf_point_on_light(scene.lights[light_id], light_point, p, scene);

        Real L_s1_pdf = light_pdf * light_point_pdf;

        Spectrum L_s1_estimate;
        Vector3 w1 = normalize(light_point.position - p);
         
        Real distance1 = distance(light_point.position, p);

        L_s1_estimate = eval(phase_function, -ray.dir, w1) * \
            emission(scene.lights[light_id], -w1, Real(0), light_point, scene) * \
            exp(-sigma_t * distance1) * \
            std::abs(dot(w1, light_point.normal)) / \
            pow(distance1 ,2);

        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);
        // return make_zero_spectrum();

    }
    else {
        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum Le = Vector3(0, 0, 0);
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return (transmittance / trans_pdf) * Le;
    }

    return make_zero_spectrum();
}

// The third volumetric renderer (not so simple anymore): 
// ...? simple above?
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    int current_medium_id = scene.camera.medium_id;

    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    int max_depth = scene.options.max_depth, rr_depth = scene.options.rr_depth, bounces = 0;
    Spectrum radiance = Vector3(0, 0, 0), current_path_throughput = Vector3(1, 1, 1), transmittance, trans_pdf, sigma_a, sigma_s, sigma_t;
    bool scatter;
    Real t_hit, t, rr_prob, pdf_phase;
    PhaseFunction phase_function;
    std::optional<Vector3> next_dir_;
    Vector3 next_dir, eval_phase;
    Real eta_scale = 1.0;

    while (true) {
        scatter = false;
        vertex_ = intersect(scene, ray, ray_diff);
        vertex = *vertex_;
        transmittance = Vector3(1, 1, 1);
        trans_pdf = Vector3(1, 1, 1);
        
        // printf("1\n");
        if (current_medium_id >= 0) {
            sigma_a = get_sigma_a(scene.media[current_medium_id], ray.org);
            sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            sigma_t = sigma_a + sigma_s;
            if (!vertex_) {
                t_hit = infinity<Real>();
            }
            else {
                t_hit = distance(ray.org, vertex.position);
            }
            t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[0];


            if (t < t_hit) {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                ray.org = ray.org + t * ray.dir;
            }
            else {
                scatter = false;
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                //
            }

        }
        current_path_throughput *= (transmittance / trans_pdf);
        
        // printf("2\n");
        if (!scatter) {
            // reach a surface, include emission
            // assertion failed: light_id >= 0
            if (is_light(scene.shapes[vertex.shape_id])) {
                radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
            }
        }


        if ((bounces == (max_depth - 1)) && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        //printf("3\n");
        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                ray.org = ray.org + t_hit * ray.dir; // corresponding to line 181
                continue;
            }
        }

        // printf("4\n");
        // sample next direct & update path throughput
        if (scatter) {
            phase_function = get_phase_function(scene.media[current_medium_id]);
            next_dir_ = sample_phase_function(phase_function, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            next_dir = *next_dir_;
            eval_phase = eval(phase_function, -ray.dir, next_dir);
            pdf_phase = pdf_sample_phase(phase_function, -ray.dir, next_dir);
            current_path_throughput *= (eval_phase / pdf_phase * sigma_s);
            // update ray.dir
            ray.dir = next_dir;
        }
        else {
            // hit a surface -- don't need to deal with this yet
            break;
        }

        // Russian Roulette
        // printf("5\n");
        rr_prob = 1;
        if (bounces >= rr_depth) {
            rr_prob = min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}


// help function of vol_path_tracing_4
Spectrum next_event_estimation(const Scene &scene, const Vector3 w, Vector3 p, const int &current_medium_id, pcg32_state &rng, const int &bounces) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Real light_pdf = light_pmf(scene, light_id);
    Vector2 rnd_param_uv = Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    Real rnd_param_w = next_pcg32_real<Real>(rng);
    PointAndNormal light_point = sample_point_on_light(scene.lights[light_id], p, rnd_param_uv, rnd_param_w, scene);
    Real light_point_pdf = pdf_point_on_light(scene.lights[light_id], light_point, p, scene);

    Vector3 p_prime = light_point.position;
    Vector3 p0 = p;

    // Vector3 T_light = Vector3(1,1,1);
    Real T_light = 1.0;
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    int max_depth = scene.options.max_depth;
    //Vector3 p_trans_dir = Vector3(1, 1,1);
    Real p_trans_dir = 1.0;

    Ray shadow_ray;
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    Real next_t;
    Spectrum sigma_a, sigma_s, sigma_t;
    PhaseFunction phase_function;

    while (true) {
        shadow_ray = Ray{p, normalize(p_prime - p), get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) * distance(p_prime, p)};
        vertex_ = intersect(scene, shadow_ray, ray_diff);
        vertex = *vertex_;
        next_t = distance(p, p_prime);
        if (vertex_) {
            next_t = distance(p, vertex.position);
        }

        if (shadow_medium_id != -1) {
            sigma_a = get_sigma_a(scene.media[shadow_medium_id], shadow_ray.org);
            sigma_s = get_sigma_s(scene.media[shadow_medium_id], shadow_ray.org);
            sigma_t = sigma_a + sigma_s; // another typo here

            T_light *= exp(-sigma_t.x * next_t);
            p_trans_dir *= exp(-sigma_t.x * next_t);
        }

        if (!vertex_) {
            break;
        }
        else {
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
            shadow_bounces += 1;
            if (max_depth != -1 && (bounces + shadow_bounces + 1) >= max_depth) {
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(vertex, shadow_ray, shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }

    }
    //https://computergraphics.stackexchange.com/questions/13320/confusion-about-pdf-defined-in-solid-angle-area-measure
    if (T_light > 0) {
        Vector3 w1 = normalize(p_prime - p0);
        Real G = std::abs(dot(light_point.normal, w1)) / pow(distance(p0, p_prime), 2);
        phase_function = get_phase_function(scene.media[shadow_medium_id]);
        Spectrum rho = eval(phase_function, w, w1);
        Real pdf_nee = light_pdf * light_point_pdf;
        Vector3 contrib = T_light * G * rho * emission(scene.lights[light_id], -w1, Real(0), light_point, scene) /pdf_nee;
        Real pdf_phase = pdf_sample_phase(phase_function, w, w1) * G * p_trans_dir;
        Real omega = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return omega * contrib;
    }

    return make_zero_spectrum();
}


// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
// the only difference: integrate next event estimation into the sampling. so others are the same
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    int current_medium_id = scene.camera.medium_id;

    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    int max_depth = scene.options.max_depth, rr_depth = scene.options.rr_depth, bounces = 0;
    Spectrum radiance = Vector3(0, 0, 0), current_path_throughput = Vector3(1, 1, 1), sigma_a, sigma_s, sigma_t;
    Real trans_pdf=1, transmittance=1;
    bool scatter;
    Real t_hit, t, rr_prob, pdf_phase;
    PhaseFunction phase_function;
    std::optional<Vector3> next_dir_;
    Vector3 next_dir, eval_phase, eval_bsdf;
    Spectrum nee;

    // the cache variables for the multiple importance sampling.
    Vector3 nee_p_cache, w1;
    //Vector3 multi_trans_pdf = Vector3(1, 1, 1);
    Real multi_trans_pdf = 1.0;
    bool never_scatter = true;
    int light_id;
    PointAndNormal light_point;
    Real pdf_nee, pdf_bsdf;
    Real dir_pdf_, G, dir_pdf = 0;
    Real omega;

    while (true) {
        scatter = false;
        vertex_ = intersect(scene, ray, ray_diff);
        vertex = *vertex_;
        transmittance = 1;
        trans_pdf = 1;

        //printf("1\n");
        if (current_medium_id >= 0) {
            sigma_a = get_sigma_a(scene.media[current_medium_id], ray.org);
            sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            sigma_t = sigma_a + sigma_s;

            if (!vertex_) {
                t_hit = infinity<Real>();
            }
            else {
                t_hit = distance(ray.org, vertex.position);
            }
            t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[0];


            if (t < t_hit) {
                scatter = true;
                never_scatter = false;
                trans_pdf = exp(-sigma_t[0] * t) * sigma_t[0];
                transmittance = exp(-sigma_t[0] * t);
                ray.org += t * ray.dir;
            }
            else {
                scatter = false;
                trans_pdf = exp(-sigma_t[0] * t_hit);
                transmittance = exp(-sigma_t[0] * t_hit);
            }

        }


        if (current_medium_id == -1) {
            if (vertex_) {
                ray.org = vertex.position;
            }
            else {
                break;
            }
        }

        current_path_throughput *= (transmittance / trans_pdf);

        //printf("2\n");
        // need modification here as well, but for multiple importance sampling weight
        if (!scatter) {
            // reach a surface, include emission
            // assertion failed: light_id >= 0
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
                }
            }
            else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    light_point.position = vertex.position;
                    light_point.normal = vertex.geometric_normal;
                    light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    pdf_nee = pdf_point_on_light(scene.lights[light_id], light_point, nee_p_cache, scene) * light_pmf(scene, light_id);
                    w1 = normalize(vertex.position - nee_p_cache);
                    G = std::abs(dot(vertex.geometric_normal, w1)) / pow(distance(nee_p_cache, vertex.position), 2);
                    dir_pdf_ = pdf_phase * multi_trans_pdf * G;
                    omega = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * omega;
                   
                }
            }

        }


        if ((bounces == (max_depth - 1)) && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        //printf("3\n");
        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                //ray.org = ray.org + t_hit * ray.dir; // corresponding to line 181
                ray = Ray{ vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>() };
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        //printf("4\n");
        // sample next direct & update path throughput
        // modification here!
        if (scatter) {
            phase_function = get_phase_function(scene.media[current_medium_id]);
            next_dir_ = sample_phase_function(phase_function, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            next_dir = *next_dir_;
            sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            eval_phase = eval(phase_function, -ray.dir, next_dir);
            pdf_phase = pdf_sample_phase(phase_function, -ray.dir, next_dir);
            multi_trans_pdf = 1;

            // add a next event estimation
            nee = next_event_estimation(scene, -ray.dir, ray.org, current_medium_id, rng, bounces);
            nee_p_cache = ray.org;
            radiance += current_path_throughput * nee * sigma_s;

            current_path_throughput *= (eval_phase / pdf_phase * sigma_s);
            // update ray.dir
            ray.dir = next_dir;
        }
        else {
            // hit a surface -- don't need to deal with this yet
            break;
        }

        // Russian Roulette
        // printf("5\n");
        rr_prob = 1;
        if (bounces >= rr_depth) {
            rr_prob = min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}


// change rho and bsdf to the bsdf version rather than phase
Spectrum next_event_estimation_bsdf(const Scene& scene, const Vector3 w, Vector3 p, const int& current_medium_id, const int &current_material_id, pcg32_state& rng, const int& bounces, const PathVertex &vertex0) {
    Material mat;
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Real light_pdf = light_pmf(scene, light_id);
    Vector2 rnd_param_uv = Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    Real rnd_param_w = next_pcg32_real<Real>(rng);
    PointAndNormal light_point = sample_point_on_light(scene.lights[light_id], p, rnd_param_uv, rnd_param_w, scene);
    Real light_point_pdf = pdf_point_on_light(scene.lights[light_id], light_point, p, scene);

    Vector3 p_prime = light_point.position;
    Vector3 p0 = p;

    // Vector3 T_light = Vector3(1,1,1);
    Real T_light = 1.0;
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    int max_depth = scene.options.max_depth;
    //Vector3 p_trans_dir = Vector3(1, 1,1);
    Real p_trans_dir = 1.0;

    Ray shadow_ray;
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    Real next_t;
    Spectrum sigma_a, sigma_s, sigma_t;
    PhaseFunction phase_function;

    while (true) {
        shadow_ray = Ray{ p, normalize(p_prime - p), get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) * distance(p_prime, p) };
        vertex_ = intersect(scene, shadow_ray, ray_diff);
        vertex = *vertex_;
        next_t = distance(p, p_prime);
        mat = scene.materials[current_material_id];
        if (vertex_) {
            next_t = distance(p, vertex.position);
        }

        if (shadow_medium_id != -1) {
            sigma_a = get_sigma_a(scene.media[shadow_medium_id], shadow_ray.org);
            sigma_s = get_sigma_s(scene.media[shadow_medium_id], shadow_ray.org);
            sigma_t = sigma_a + sigma_s; // another typo here

            T_light *= exp(-sigma_t.x * next_t);
            p_trans_dir *= exp(-sigma_t.x * next_t);
        }

        if (!vertex_) {
            break;
        }
        else {
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
            shadow_bounces += 1;
            if (max_depth != -1 && (bounces + shadow_bounces + 1) >= max_depth) {
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(vertex, shadow_ray, shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }

    }
    //https://computergraphics.stackexchange.com/questions/13320/confusion-about-pdf-defined-in-solid-angle-area-measure
    if (T_light > 0) {
        Vector3 w1 = normalize(p_prime - p0);
        Real G = std::abs(dot(light_point.normal, w1)) / pow(distance(p0, p_prime), 2);
        phase_function = get_phase_function(scene.media[shadow_medium_id]);
        // change the rho here
        //Spectrum rho = eval(phase_function, w, w1);
        //mat = scene.materials[current_material_id];
        mat = scene.materials[vertex0.material_id];
        Spectrum rho = eval(mat, w, w1, vertex0, scene.texture_pool);
        Real pdf_nee = light_pdf * light_point_pdf;
        Vector3 contrib = T_light * G * rho * emission(scene.lights[light_id], -w1, Real(0), light_point, scene) / pdf_nee;
        Real pdf_bsdf = pdf_sample_bsdf(mat, w, w1, vertex0, scene.texture_pool) * G * p_trans_dir;
        Real omega = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
        return omega * contrib;
    }

    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
// the difference: 1) hit a surface; 2) change the phase to include the bsdf
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    int current_medium_id = scene.camera.medium_id;
    int current_material_id;
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    int max_depth = scene.options.max_depth, rr_depth = scene.options.rr_depth, bounces = 0;
    Spectrum radiance = Vector3(0, 0, 0), current_path_throughput = Vector3(1, 1, 1), sigma_a, sigma_s, sigma_t;
    Real trans_pdf = 1, transmittance = 1;
    bool scatter;
    Real t_hit, t, rr_prob, pdf_phase, pdf_bsdf, dir_pdf;
    PhaseFunction phase_function;
    std::optional<Vector3> next_dir_;
    std::optional<BSDFSampleRecord> bsdf_;
    BSDFSampleRecord bsdf;
    Vector3 next_dir, eval_phase, eval_bsdf;
    Spectrum nee, nee_bsdf;

    Vector3 nee_p_cache, w1;
    Real multi_trans_pdf = 1.0;
    bool never_scatter = true;
    int light_id;
    PointAndNormal light_point;
    Real pdf_nee, pdf_nee_bsdf;
    Real dir_pdf_, G;
    Real omega;

    while (true) {
        scatter = false;
        vertex_ = intersect(scene, ray, ray_diff);
        vertex = *vertex_;
        transmittance = 1;
        trans_pdf = 1;

        // printf("1\n");
        if (current_medium_id >= 0) {
            sigma_a = get_sigma_a(scene.media[current_medium_id], ray.org);
            sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            sigma_t = sigma_a + sigma_s;

            
            if (!vertex_) {
                t_hit = infinity<Real>();
            }
            else {
                t_hit = distance(ray.org, vertex.position);
            }
            t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[0];


            if (t < t_hit) {
                scatter = true;
                never_scatter = false;
                trans_pdf = exp(-sigma_t[0] * t) * sigma_t[0];
                transmittance = exp(-sigma_t[0] * t);
                ray.org += t * ray.dir;
            }
            else {
                scatter = false;
                trans_pdf = exp(-sigma_t[0] * t_hit);
                transmittance = exp(-sigma_t[0] * t_hit);
                ray.org += t_hit * ray.dir;
            }

        }


        if (current_medium_id == -1) {
            if (vertex_) {
                ray.org = vertex.position;
            }
            else {
                break;
            }
        }

        current_path_throughput *= (transmittance / trans_pdf);

        // printf("2\n");
        if (!scatter) {
            // reach a surface, include emission
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
                }
            }
            else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    light_point.position = vertex.position;
                    light_point.normal = vertex.geometric_normal;
                    light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    pdf_nee = pdf_point_on_light(scene.lights[light_id], light_point, nee_p_cache, scene) * light_pmf(scene, light_id);
                    w1 = normalize(vertex.position - nee_p_cache);
                    G = std::abs(dot(vertex.geometric_normal, w1)) / pow(distance(nee_p_cache, vertex.position), 2);
                    dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    omega = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * omega;

                }
            }

        }


        if ((bounces == (max_depth - 1)) && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        // printf("3\n");
        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                //ray.org = ray.org + t_hit * ray.dir; // corresponding to line 181
                ray = Ray{ vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>() };
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        // printf("4\n");
        // sample next direct & update path throughput
        
        if (scatter) {
            never_scatter = false;
            phase_function = get_phase_function(scene.media[current_medium_id]);
            next_dir_ = sample_phase_function(phase_function, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            next_dir = *next_dir_;
            
            sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            eval_phase = eval(phase_function, -ray.dir, next_dir);
            pdf_phase = pdf_sample_phase(phase_function, -ray.dir, next_dir);
            multi_trans_pdf = 1;

            // add a next event estimation
            nee = next_event_estimation(scene, -ray.dir, ray.org, current_medium_id, rng, bounces);
            nee_p_cache = ray.org;
            radiance += current_path_throughput * nee * sigma_s;

            current_path_throughput *= (eval_phase / pdf_phase * sigma_s);
            dir_pdf = pdf_phase;
            // update ray.dir
            ray.dir = next_dir;
        }
        else {
            never_scatter = false;
            int current_material_id = vertex.material_id;
            Spectrum nee_bsdf = next_event_estimation_bsdf(scene, -ray.dir, ray.org, current_medium_id, current_material_id, rng, bounces, vertex);

            bsdf_ = sample_bsdf(scene.materials[current_material_id], -ray.dir, vertex, scene.texture_pool, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng));
            if (!bsdf_) {
                break;
            }
            else {
                bsdf = *bsdf_;
            }

            Spectrum eval_bsdf = eval(scene.materials[current_material_id], -ray.dir, bsdf.dir_out, vertex, scene.texture_pool);
            Real pdf_bsdf = pdf_sample_bsdf(scene.materials[current_material_id], -ray.dir, bsdf.dir_out, vertex, scene.texture_pool);

            ray = Ray{ ray.org, bsdf.dir_out, get_intersection_epsilon(scene), infinity<Real>() };
            current_medium_id = update_medium(vertex, ray, current_medium_id);

            dir_pdf = pdf_bsdf;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1.0;

            radiance += current_path_throughput * nee_bsdf;
            current_path_throughput *= eval_bsdf / pdf_bsdf;

        }
        

        // Russian Roulette
        // printf("5\n");
        rr_prob = 1;
        if (bounces >= rr_depth) {
            rr_prob = min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}


Spectrum get_nee(const Scene& scene, const Vector3 w, Vector3 p, const int& current_medium_id, pcg32_state& rng, const int& bounces) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Real light_pdf = light_pmf(scene, light_id);
    Vector2 rnd_param_uv = Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    Real rnd_param_w = next_pcg32_real<Real>(rng);
    PointAndNormal light_point = sample_point_on_light(scene.lights[light_id], p, rnd_param_uv, rnd_param_w, scene);
    Real light_point_pdf = pdf_point_on_light(scene.lights[light_id], light_point, p, scene);

    Vector3 p_prime = light_point.position;
    Vector3 p0 = p;

    Spectrum T_light = make_const_spectrum(1);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    int max_depth = scene.options.max_depth;
    //Vector3 p_trans_dir = Vector3(1, 1,1);
    Spectrum p_trans_nee = make_const_spectrum(1);
    Spectrum p_trans_dir = make_const_spectrum(1);
    Spectrum real_prob;

    Real u, t, dt, accum_t;
    int channel, iteration;
    int max_null_collisions = scene.options.max_null_collisions;
    Spectrum majorant;

    Ray shadow_ray;
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    Real next_t;
    Spectrum sigma_a, sigma_s, sigma_t;
    PhaseFunction phase_function;

    while (true) {
        shadow_ray = Ray{ p, normalize(p_prime - p), get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) * distance(p_prime, p) };
        vertex_ = intersect(scene, shadow_ray, ray_diff);
        
        next_t = distance(p, p_prime);
        if (vertex_) {
            vertex = *vertex_;
            next_t = distance(p, vertex.position);
        }

        if (shadow_medium_id != -1) {
            u = next_pcg32_real<Real>(rng);
            channel = std::clamp(int(u * 3), 0, 2);
            iteration = 0;
            dt = 0;
            accum_t = 0;
            majorant = get_majorant(scene.media[shadow_medium_id], shadow_ray);
            while (true) {
                // to avoid strange white line
                //if (majorant[channel] <= 0) {
                //    break;
                //}
                if (iteration >= max_null_collisions) {
                    break;
                }
                t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                dt = next_t - accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    // did not hit the surface, so this is a null-scattering event
                    sigma_a = get_sigma_a(scene.media[shadow_medium_id], shadow_ray.org + accum_t * shadow_ray.dir);
                    sigma_s = get_sigma_s(scene.media[shadow_medium_id], shadow_ray.org + accum_t * shadow_ray.dir);
                    sigma_t = sigma_a + sigma_s;
                    T_light *= (exp(-majorant * t) * (majorant - sigma_t) / max(majorant));
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                    if (max(T_light) <= 0) {
                        //break;
                        break;
                    }
                }
                else {
                    // hit the surface
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    break;
                }
                iteration += 1;
            }
        }

        if (!vertex_) {
            break;
        }
        else {
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
            shadow_bounces += 1;
            if (max_depth != -1 && (bounces + shadow_bounces + 1) >= max_depth) {
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(vertex, shadow_ray, shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }

    }
    //https://computergraphics.stackexchange.com/questions/13320/confusion-about-pdf-defined-in-solid-angle-area-measure
    if (max(T_light) > 0) {
        Vector3 w1 = normalize(p_prime - p0);
        Real G = std::abs(dot(light_point.normal, w1)) / pow(distance(p0, p_prime), 2);
        phase_function = get_phase_function(scene.media[shadow_medium_id]);
        Spectrum rho = eval(phase_function, w, w1);
        Spectrum pdf_nee = light_pdf * light_point_pdf * p_trans_nee;
        Vector3 contrib = T_light * G * rho * emission(scene.lights[light_id], -w1, Real(0), light_point, scene) / avg(pdf_nee);
        Spectrum pdf_phase = pdf_sample_phase(phase_function, w, w1) * G * p_trans_dir;
        Spectrum omega = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return omega * contrib;
    }

    return make_zero_spectrum();
}

Spectrum get_nee_bsdf(const Scene& scene, const Vector3 w, Vector3 p, const int& current_medium_id, const int& current_material_id, pcg32_state& rng, const int& bounces, const PathVertex& vertex0) {
    Material mat;
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Real light_pdf = light_pmf(scene, light_id);
    Vector2 rnd_param_uv = Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    Real rnd_param_w = next_pcg32_real<Real>(rng);
    PointAndNormal light_point = sample_point_on_light(scene.lights[light_id], p, rnd_param_uv, rnd_param_w, scene);
    Real light_point_pdf = pdf_point_on_light(scene.lights[light_id], light_point, p, scene);

    Vector3 p_prime = light_point.position;
    Vector3 p0 = p;

    // Vector3 T_light = Vector3(1,1,1);
    Spectrum T_light = make_const_spectrum(1);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    int max_depth = scene.options.max_depth;
    //Vector3 p_trans_dir = Vector3(1, 1,1);
    Spectrum p_trans_nee = make_const_spectrum(1);
    Spectrum p_trans_dir = make_const_spectrum(1);
    Spectrum real_prob;

    Real u, t, dt, accum_t;
    int channel, iteration;
    int max_null_collisions = scene.options.max_null_collisions;
    Spectrum majorant;

    Ray shadow_ray;
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    Real next_t;
    Spectrum sigma_a, sigma_s, sigma_t;
    PhaseFunction phase_function;

    while (true) {
        shadow_ray = Ray{ p, normalize(p_prime - p), get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) * distance(p_prime, p) };
        vertex_ = intersect(scene, shadow_ray, ray_diff);
        vertex = *vertex_;
        next_t = distance(p, p_prime);
        mat = scene.materials[current_material_id];
        if (vertex_) {
            next_t = distance(p, vertex.position);
        }

        if (shadow_medium_id != -1) {
            u = next_pcg32_real<Real>(rng);
            channel = std::clamp(int(u * 3), 0, 2);
            iteration = 0;
            dt = 0;
            accum_t = 0;
            majorant = get_majorant(scene.media[shadow_medium_id], shadow_ray);
            while (true) {
                // to avoid strange white line
                //if (majorant[channel] <= 0) {
                //    break;
                //}
                if (iteration >= max_null_collisions) {
                    break;
                }
                t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                dt = next_t - accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    sigma_a = get_sigma_a(scene.media[shadow_medium_id], shadow_ray.org + accum_t * shadow_ray.dir);
                    sigma_s = get_sigma_s(scene.media[shadow_medium_id], shadow_ray.org + accum_t * shadow_ray.dir);
                    sigma_t = sigma_a + sigma_s;
                    T_light *= (exp(-majorant * t) * (majorant - sigma_t) / max(majorant));
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                    if (max(T_light) <= 0) {
                        break;
                    }
                }
                else {
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    break;
                }
                iteration += 1;
            }
        }

        if (!vertex_) {
            break;
        }
        else {
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
            shadow_bounces += 1;
            if (max_depth != -1 && (bounces + shadow_bounces + 1) >= max_depth) {
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(vertex, shadow_ray, shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }

    }
    //https://computergraphics.stackexchange.com/questions/13320/confusion-about-pdf-defined-in-solid-angle-area-measure
    if (max(T_light) > 0) {
        Vector3 w1 = normalize(p_prime - p0);
        Real G = std::abs(dot(light_point.normal, w1)) / pow(distance(p0, p_prime), 2);
        mat = scene.materials[vertex0.material_id];
        Spectrum rho = eval(mat, w, w1, vertex0, scene.texture_pool);
        Spectrum pdf_nee = light_pdf * light_point_pdf * p_trans_nee;
        Vector3 contrib = T_light * G * rho * emission(scene.lights[light_id], -w1, Real(0), light_point, scene) / avg(pdf_nee);
        Spectrum pdf_bsdf = pdf_sample_bsdf(mat, w, w1, vertex0, scene.texture_pool) * G * p_trans_dir;
        Spectrum omega = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
        return omega * contrib;
    }

    return make_zero_spectrum();
}


// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
// ......
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    int current_medium_id = scene.camera.medium_id;
    int current_material_id;
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    int max_depth = scene.options.max_depth, rr_depth = scene.options.rr_depth, bounces = 0;
    Spectrum radiance = Vector3(0, 0, 0), current_path_throughput = Vector3(1, 1, 1), sigma_a, sigma_s, sigma_t;
    Spectrum trans_dir_pdf, trans_nee_pdf, transmittance; // from real to spectrum
    Real dir_pdf = 0;
    Spectrum majorant, real_prob;
    Real u;
    int channel, iteration, max_null_collisions = scene.options.max_null_collisions;
    Real accum_t, dt;
    bool scatter;
    Real t_hit, t, rr_prob, pdf_phase, pdf_bsdf;
    PhaseFunction phase_function;
    std::optional<Vector3> next_dir_;
    std::optional<BSDFSampleRecord> bsdf_;
    BSDFSampleRecord bsdf;
    Vector3 next_dir, eval_phase, eval_bsdf;
    Spectrum nee, nee_bsdf;

    Vector3 nee_p_cache, w1;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    bool never_scatter = true;
    int light_id;
    PointAndNormal light_point;
    Real pdf_nee_bsdf;
    Real G;
    Spectrum dir_pdf_, omega, pdf_nee;

    Real eta_scale = 1.0;

    
    while (true) {
        scatter = false;
        vertex_ = intersect(scene, ray, ray_diff);
        if (vertex_) {
            vertex = *vertex_;
            t_hit = distance(ray.org, vertex.position);
        }
        else {
            t_hit = infinity<Real>();
        }
        
        transmittance = make_const_spectrum(1);
        trans_dir_pdf = make_const_spectrum(1);
        trans_nee_pdf = make_const_spectrum(1);

        //printf("1\n");
        // the free flight sampling position - here - in mudium
        
        if (current_medium_id >= 0) {
            majorant = get_majorant(scene.media[current_medium_id], ray);

            u = next_pcg32_real<Real>(rng);
            channel = std::clamp(int(u * 3), 0, 2);

            accum_t = 0;
            iteration = 0;

           
            while (true) {
                // it brings strange white line.
                //if (majorant[channel] <= 0) {
                //    break;
                //}
                if (iteration >= max_null_collisions) {
                    break;
                }
                t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                dt = t_hit - accum_t;
                // update accumulated distance
                accum_t = min(accum_t + t, t_hit);
                sigma_a = get_sigma_a(scene.media[current_medium_id], ray.org + accum_t * ray.dir);
                sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org + accum_t * ray.dir);
                sigma_t = sigma_a + sigma_s;
                if (t < dt) { // have not reached the surface
                    // sample from real/fake particle events
                    real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        // hit a real particle
                        scatter = true;
                        never_scatter = false;
                        transmittance *= (exp(-majorant * t) / max(majorant));
                        trans_dir_pdf *= (exp(-majorant * t) * majorant * real_prob / max(majorant));
                        // don't need to account for trans_nee_pdf since we scatter
                        ray.org += accum_t * ray.dir;
                        break;

                    }
                    else {
                        // hit a fake particle
                        transmittance *= (exp(-majorant * t) * (majorant - sigma_t)/ max(majorant));
                        trans_dir_pdf *= (exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant));
                        trans_nee_pdf *= (exp(-majorant * t) * majorant / max(majorant));
                    }
                }
                else { // reach the surface
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    ray.org = vertex.position;
                    break;
                }
                iteration += 1;
            }
            multi_trans_pdf *= trans_dir_pdf;
        }

        //yes
        if (current_medium_id == -1) {
            if (vertex_) {
                ray.org = vertex.position;
            }
            else {
                break; 
            }
        }
        
        


        current_path_throughput *= (transmittance / avg(trans_dir_pdf));

        //printf("2\n");
        if (!scatter && vertex_) {
            // if (!vertex_) {
            //     printf("here1!\n");
            //     printf("%d", vertex.shape_id);
            // }
            // reach a surface, include emission
            if (bounces == 0) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
                }
            }
            else {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    light_point.position = vertex.position;
                    light_point.normal = vertex.geometric_normal;
                    light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    pdf_nee = pdf_point_on_light(scene.lights[light_id], light_point, nee_p_cache, scene) * light_pmf(scene, light_id) * trans_nee_pdf;
                    w1 = normalize(vertex.position - nee_p_cache);
                    
                    // G problem
                    G = std::abs(dot(vertex.geometric_normal, w1)) / pow(distance(nee_p_cache, vertex.position), 2);

                    dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    omega = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * omega;

                }
            }

        }


        if ((bounces >= (max_depth - 1)) && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        // printf("3\n");
        if (!scatter && vertex_) {
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                ray = Ray{ vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>() };
                multi_trans_pdf *= trans_dir_pdf;
                continue;
            }
        }

        //printf("4\n");
        // sample next direct & update path throughput
        if (scatter && current_medium_id != -1) {
            never_scatter = false;
            sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            nee = get_nee(scene, -ray.dir, ray.org, current_medium_id, rng, bounces);

            radiance += current_path_throughput * nee * sigma_s;

            phase_function = get_phase_function(scene.media[current_medium_id]);
            next_dir_ = sample_phase_function(phase_function, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            next_dir = *next_dir_;
            
            if (!next_dir_) {
                break;
            }
            
            // add a next event estimation



            eval_phase = eval(phase_function, -ray.dir, next_dir);
            pdf_phase = pdf_sample_phase(phase_function, -ray.dir, next_dir);
            multi_trans_pdf = make_const_spectrum(1);



            current_path_throughput *= (eval_phase / pdf_phase * sigma_s);

            // update ray.dir
            dir_pdf = pdf_phase;
            ray.dir = next_dir;
            nee_p_cache = ray.org;
        }
        else {
            // if (!vertex_) {
            //   //printf("here\n");
            //     continue;
            // }
            //if (vertex_) {
            if (vertex_) {
                never_scatter = false;
                int current_material_id = vertex_->material_id;
                Spectrum nee_bsdf = get_nee_bsdf(scene, -ray.dir, ray.org, current_medium_id, current_material_id, rng, bounces, *vertex_);
                radiance += current_path_throughput * nee_bsdf;
                
                bsdf_ = sample_bsdf(scene.materials[current_material_id], -ray.dir, *vertex_, scene.texture_pool, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng));
                if (!bsdf_) {
                    break;
                }
                else {
                    bsdf = *bsdf_;
                }

                
                Spectrum eval_bsdf = eval(scene.materials[current_material_id], -ray.dir, bsdf.dir_out, *vertex_, scene.texture_pool);
                Real pdf_bsdf = pdf_sample_bsdf(scene.materials[current_material_id], -ray.dir, bsdf.dir_out, *vertex_, scene.texture_pool);

                if (bsdf.eta == 0) {
                    ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf.roughness);
                }
                else {
                    ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf.eta, bsdf.roughness);
                    eta_scale /= (bsdf.eta * bsdf.eta);
                    current_medium_id = update_medium(vertex, ray, current_medium_id);
                }

                ray = Ray{ ray.org, bsdf.dir_out, get_intersection_epsilon(scene), infinity<Real>() };
                nee_p_cache = ray.org;
                multi_trans_pdf = make_const_spectrum(1);
                dir_pdf = pdf_bsdf;

                current_path_throughput *= eval_bsdf / pdf_bsdf;



            }

        }



        // Russian Roulette
        // printf("5\n");
        rr_prob = 1;
        if (bounces >= rr_depth) {
            rr_prob = min(max(current_path_throughput/eta_scale), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }

    return radiance;
}
