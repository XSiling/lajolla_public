#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // the disney bsdf is the combination of different materials
    // restrucutre tbd
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    //Real eta = bsdf.eta;

    // the clamp of some parameters to avoid some issues
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // common parameters
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real HDotL = dot(half_vector, dir_out);
    Real HDotV = dot(half_vector, dir_in);
    Real NDotL = dot(frame.n, dir_out);
    Real NDotV = dot(frame.n, dir_in);
    Vector3 Hl = to_local(frame, half_vector);
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);

    // calculate the separate materials BRDF
    Spectrum f_diffuse = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });

    Spectrum f_sheen = operator()(DisneySheen{ bsdf.base_color, bsdf.sheen_tint });

    Spectrum f_metal; // to be edit Fm
    Real Dm = GGX_metal(Hl, alphax, alphay);
    Real Gm = Smith_metal(to_local(frame, dir_in), alphax, alphay) * Smith_metal(to_local(frame, dir_out), alphax, alphay);
    Vector3 Ctint = luminance(base_color) > 0 ? (base_color / luminance(base_color)) : Vector3(1.0, 1.0, 1.0);
    Spectrum Ks = (1 - specular_tint) + specular_tint * Ctint;
    Spectrum C0 = specular * R0(eta) * (1 - metallic) * Ks + metallic * base_color;
    Spectrum Fm = C0 + (1 - C0) * pow(1 - HDotL, 5);
    f_metal = Fm * Dm * Gm / (4 * std::abs(NDotV));

    Spectrum f_clearcoat = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
    Spectrum f_glass = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });

    // in-material judgement
    if (dot(dir_in, frame.n) <= 0) {
        f_diffuse = make_zero_spectrum();
        f_metal = make_zero_spectrum();
        f_clearcoat = make_zero_spectrum();
        f_sheen = make_zero_spectrum();
    }
    else {
        if (dot(vertex.geometric_normal, dir_out) <= 0) {
            f_metal = make_zero_spectrum();
        }
    }


    return (1 - specular_transmission) * (1 - metallic) * f_diffuse \
        + (1 - metallic) * sheen * f_sheen \
        + (1 - specular_transmission * (1 - metallic)) * f_metal \
        + 0.25 * clearcoat * f_clearcoat \
        + (1 - metallic) * specular_transmission * f_glass;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);


    // the weight
    Real diffuse_weight = (1 - metallic) * (1 - specular_transmission);
    Real metal_weight = (1 - specular_transmission * (1 - metallic));
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real clearcoat_weight = 0.25 * clearcoat;
    Real sum = diffuse_weight + metal_weight + glass_weight + clearcoat_weight;
    diffuse_weight /= sum;
    metal_weight /= sum;
    glass_weight /= sum;
    clearcoat_weight /= sum;

    // importance sampling.
    Real pdf_diffuse = operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    Real pdf_clearcoat = operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
    Real pdf_metal = operator()(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    Real pdf_glass = operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });

    // dot(dir_in, frame.n<=0)
    if (dot(dir_in, vertex.geometric_normal) <= 0) {
        return pdf_glass;
    }

    return pdf_diffuse * diffuse_weight + pdf_glass * glass_weight + pdf_clearcoat * clearcoat_weight + pdf_metal * metal_weight;
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyBSDF& bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    // the weight
    Real diffuse_weight = (1 - metallic) * (1 - specular_transmission);
    Real metal_weight = (1 - specular_transmission * (1 - metallic));
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real clearcoat_weight = 0.25 * clearcoat;
    Real sum = diffuse_weight + metal_weight + glass_weight + clearcoat_weight;
    diffuse_weight /= sum;
    metal_weight /= sum;
    glass_weight /= sum;
    clearcoat_weight /= sum;

    metal_weight += diffuse_weight;
    glass_weight += metal_weight;
    clearcoat_weight += glass_weight;

    assert(std::abs(clearcoat_weight - 1) < 0.001);

    if (dot(dir_in, vertex.geometric_normal) <= 0) {
        return operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }

    // control the random number seed;
    if (rnd_param_w < diffuse_weight) {
        return operator()(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    }
    if (rnd_param_w < metal_weight) {
        return operator()(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    }
    if (rnd_param_w < glass_weight) {
        return operator()(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }
    if (rnd_param_w < clearcoat_weight) {
        return operator()(DisneyClearcoat{ bsdf.clearcoat_gloss });
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF& bsdf) const {
    return bsdf.base_color;
}