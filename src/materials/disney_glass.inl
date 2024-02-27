#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // refer to homework1 and roughdielectric.inl
    // src: https://github.com/3AceShowHand/LuisaCompute/blob/d39996c05e7e1029259f2917750e71f6f7550310/src/py/disney.py#L309
    // src: https://github.com/NJUCG/Moer/blob/03f7ae2931da63beca710a1acb17df8c5dedf5fe/src/FunctionLayer/Material/DisneyBSDF.cpp#L17

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    roughness = std::clamp(roughness, Real(0.01), Real(1));
    //anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    else {
        half_vector = normalize(dir_in + dir_out * eta);
    }

    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real HDotV = dot(half_vector, dir_in);
    Real HDotL = dot(half_vector, dir_out);
    Real NDotV = dot(frame.n, dir_in);
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);
   
    //Real Rs = (HDotV - eta * HDotL) / (HDotV + eta * HDotL);
    //Real Rp = (eta * HDotV - HDotL) / (eta * HDotV + HDotL);
    // this version's result is weird
    //Real Fg = fresnel_dielectric(HDotV, HDotL, eta);
    Real Fg = fresnel_dielectric(HDotV ,eta);
    // Dg and Gg are the same as the metal situation;
    Vector3 Hl = to_local(frame, half_vector);
    Real Dg = GGX_metal(Hl, alphax, alphay);
    Real Gg = Smith_metal(to_local(frame, dir_in),alphax, alphay) * Smith_metal(to_local(frame, dir_out), alphax, alphay);



    //if (dot(frame.n, dir_in) * dot(frame.n, dir_out) > 0) {
    if (reflect) {
        return base_color * Fg * Dg * Gg / (4 * std::abs(NDotV));
    }
    else {
        return sqrt(base_color) * (1 - Fg) * Dg * Gg * std::abs(HDotL) * std::abs(HDotV) / (std::abs(NDotV) * sqr(HDotV + eta * HDotL));
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    roughness = std::clamp(roughness, Real(0.01), Real(1));
    //anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    assert(eta > 0);
    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    else {
        half_vector = normalize(dir_in + dir_out * eta);
    }

    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);
    Vector3 Hl = to_local(frame, half_vector);
    Real HDotV = dot(half_vector, dir_in);
    Real HDotL = dot(half_vector, dir_out);
    Real F = fresnel_dielectric(HDotV, eta);
    Real D = GGX_metal(Hl, alphax, alphay);
    Real G_in = Smith_metal(to_local(frame, dir_in), alphax, alphay);
    Real G_out = Smith_metal(to_local(frame, dir_out), alphax, alphay);


    if (reflect) {
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    }
    else {
        Real HDotL = dot(half_vector, dir_out);
        Real sqrt_denom = HDotV + eta * HDotL;
        Real dh_dout = eta * eta * HDotL / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G_in * fabs(dh_dout * HDotV / dot(frame.n, dir_in));
    }


    return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    roughness = std::clamp(roughness, Real(0.01), Real(1));
    //anisotropic = std::clamp(anisotropic, Real(0.01), Real(1));

    // edit here to change the alpha
    Real aspect = std::sqrt(1 - 0.9 * anisotropic);
    Real alphax = std::max(0.0001, roughness * roughness / aspect);
    Real alphay = std::max(0.0001, roughness * roughness * aspect);
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal = sample_visible_normals_metal(local_dir_in, alphax, alphay, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    Real HDotV = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(HDotV, eta);

    if (rnd_param_w <= F) {
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{ reflected, 0, roughness };
    }
    else {
        Real h_dot_out_sq = 1 - (1 - HDotV * HDotV) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            return {};
        }

        if (HDotV < 0) {
            half_vector = -half_vector;
        }
        Real HDotL = sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(HDotV) / eta - HDotL) * half_vector;
        return BSDFSampleRecord{ refracted, eta, roughness };

    }
    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
