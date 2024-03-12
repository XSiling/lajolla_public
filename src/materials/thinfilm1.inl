#include "../microfacet.h"

// reference source: https://github.com/knightcrawler25/GLSL-PathTracer/blob/master/src/shaders/common/disney.glsl
// line 90

Spectrum eval_op::operator()(const ThinFilm& bsdf) const {
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



    float eta_2 = mix(1.0, eta2, smoothstep(0.0, 0.03, Dinc));

    // Compute dot products
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if (NdotL < 0 || NdotV < 0) return vec3(0);
    vec3 H = normalize(L + V);
    float NdotH = dot(N, H);
    float cosTheta1 = dot(H, L);
    float cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1 - sqr(cosTheta1)));

    // First interface
    vec2 R12, phi12;
    fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
    vec2 R21 = R12;
    vec2 T121 = vec2(1.0) - R12;
    vec2 phi21 = vec2(PI) - phi12;

    // Second interface
    vec2 R23, phi23;
    fresnelConductor(cosTheta2, eta_2, eta3, kappa3, R23, phi23);

    // Phase shift
    float OPD = Dinc * cosTheta2;
    vec2 phi2 = phi21 + phi23;

    // Compound terms
    vec3 I = vec3(0);
    vec2 R123 = R12 * R23;
    vec2 r123 = sqrt(R123);
    vec2 Rs = sqr(T121) * R23 / (1 - R123);

    // Reflectance term for m=0 (DC term amplitude)
    vec2 C0 = R12 + Rs;
    vec3 S0 = evalSensitivity(0.0, 0.0);
    I += depol(C0) * S0;

    // Reflectance term for m>0 (pairs of diracs)
    vec2 Cm = Rs - T121;
    for (int m = 1; m <= 3; ++m) {
        Cm *= r123;
        vec3 SmS = 2.0 * evalSensitivity(m * OPD, m * phi2.x);
        vec3 SmP = 2.0 * evalSensitivity(m * OPD, m * phi2.y);
        I += depolColor(Cm.x * SmS, Cm.y * SmP);
    }

    // Convert back to RGB reflectance
    I = clamp(XYZ_TO_RGB * I, vec3(0.0), vec3(1.0));

    // Microfacet BRDF formula
    float D = GGX(NdotH, alpha);
    float G = smithG_GGX(NdotL, NdotV, alpha);
    return D * G * I / (4 * NdotL * NdotV);
}

Real pdf_sample_bsdf_op::operator()(const ThinFilm& bsdf) const {
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
    Real roughness = 0.5;
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = 0.8;
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
sample_bsdf_op::operator()(const ThinFilm& bsdf) const {
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
    Real roughness = 0.5;
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = 0.8;
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

TextureSpectrum get_texture_op::operator()(const ThinFilm& bsdf) const {
    return {};
}
