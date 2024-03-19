#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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

    // the pdf solution is totally a mess.
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    //clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));
    Real alpha = (Real(1) - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real HDotL = dot(half_vector, dir_out);
    Real r0 = R0(Real(1.5));
    Vector3 Hl = to_local(frame, half_vector);

    Real Dr = GTR1(Hl.z, alpha);
    Real Fr = r0 + (1-r0) * pow(1 - std::abs(HDotL), 5);
    Real Gr = Smith_metal(to_local(frame, dir_in), 0.25, 0.25) * Smith_metal(to_local(frame, dir_out), 0.25, 0.25);

    Real NDotV = dot(frame.n, dir_in);
    
    Real f = Fr * Gr * Dr / (4 * std::abs(NDotV));
    // solve divisionbyzero problem
    return Vector3(f, f, f);

}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat& bsdf) const {
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

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    //clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));
    Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
    Vector3 Hl = to_local(frame, half_vector);
    Real Dr = GTR1(Hl.z, alpha);

    Real HDotL = dot(half_vector, dir_out);
    Real NDotH = dot(frame.n, half_vector);
    return Dr * abs(NDotH) / (Real(4) * abs(HDotL));

    // Homework 1: implement this!
    // source: https://github.com/mmp/pbrt-v3/blob/master/src/materials/disney.cpp
    // source: https://schuttejoe.github.io/post/disneybsdf/

}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    //clearcoat_gloss = std::clamp(clearcoat_gloss, Real(0.01), Real(1));

    Real alphag = (Real(1) - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
    Real alpha2 = alphag * alphag; //cannot believe the typo error here...

    Real cosTheta = std::sqrt((Real(1) - std::pow(alpha2, Real(1) - rnd_param_uv.x)) / (Real(1) - alpha2));
    Real sinTheta = std::sqrt(Real(1) - cosTheta * cosTheta);
    Real phi = 2 * c_PI * rnd_param_uv.y;
    Vector3 wm = Vector3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
    //wm = normalize(wm);


    Vector3 half_vector = to_world(frame, wm);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0), alphag
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}

