#include "../microfacet.h"

// reference source: https://github.com/knightcrawler25/GLSL-PathTracer/blob/master/src/shaders/common/disney.glsl
// line 90

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    roughness = std::clamp(roughness, Real(0.01), Real(1));
    //anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real HDotL = dot(half_vector, dir_out);

    // get Fm
    Spectrum Fm = schlick_fresnel(base_color, std::abs(HDotL));
    //base_color + (1 - base_color) * pow(1 - std::abs(HDotL), 5);
    
    // get Dm
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);

    // github version x hw1 version YES quite different
    Real NDotH = dot(frame.n, half_vector);
    Vector3 Hl = to_local(frame, half_vector);
    Real Dm = GGX_metal(Hl, alphax, alphay);

    // get Gm
    Real NDotV = dot(frame.n, dir_in);
    Real GV = Smith_metal(to_local(frame,dir_in), alphax, alphay);
    Real GL = Smith_metal(to_local(frame,dir_out), alphax, alphay);
    Real Gm = GV * GL;
    Spectrum R = Fm * Dm * Gm;
    return R / (4 * std::abs(NDotV));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // same as roughplastic
    // currently exactly copy
    
    Vector3 half_vector = normalize(dir_in + dir_out);
    Vector3 Hl = to_local(frame, half_vector);
    Real NDotV = dot(frame.n, dir_in);
    Real NDotL = dot(frame.n, dir_out);
    Real NDotH = dot(frame.n, half_vector);

    // no specular_reflectance and diffuse reflectance here
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);

    Real G = Smith_metal(to_local(frame, dir_in), alphax, alphay);
    Real D = GGX_metal(Hl, alphax, alphay);
    Real spec_prob = (G * D) / (4 * NDotV);
    return spec_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // use the part of specular of roughplastic
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    //Real alpha = roughness * roughness;
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);
    Vector3 local_micro_normal = sample_visible_normals_metal(local_dir_in, alphax, alphay, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0), roughness
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
