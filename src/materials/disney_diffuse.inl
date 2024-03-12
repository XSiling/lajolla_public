Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Vector3 half_vector = normalize(dir_in + dir_out);

    Real HDotL = std::abs(dot(half_vector, dir_out));
    Real Fd90 = 0.5 + 2 * roughness * HDotL * HDotL;

    Real NDotL = std::abs(dot(frame.n, dir_out));
    Real NDotV = std::abs(dot(frame.n, dir_in));
    Real Fd = schlick_fresnel_approximation(Fd90, NDotV) * schlick_fresnel_approximation(Fd90, NDotL);

    Spectrum fbD = base_color / c_PI * Fd * NDotL;

    Real Fss90 = roughness * HDotL * HDotL;
    Real Fss = schlick_fresnel_approximation(Fss90, NDotV) * schlick_fresnel_approximation(Fss90, NDotL);
    Spectrum fss = 1.25 * base_color / c_PI * (Fss * (Real(1)/(NDotL + NDotV) - Real(0.5)) + Real(0.5)) * NDotL;
    return (1 - subsurface) * fbD + subsurface * fss;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    
    // cosine hemisphere
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: copy lambertian
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0), roughness
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
