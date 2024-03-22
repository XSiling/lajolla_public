#include "../microfacet.h"

Spectrum eval_op::operator()(const ThinFilm& bsdf) const {

	if (dot(vertex.geometric_normal, dir_in) < 0 ||
		dot(vertex.geometric_normal, dir_out) < 0) {
		return make_zero_spectrum();
	}

	Frame frame = vertex.shading_frame;
	if (dot(frame.n, dir_in) < 0) {
		frame = -frame;
	}
	Vector3 wavelengths = { 580.0, 550.0, 450.0 };

	Real alpha = eval(bsdf.alpha, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real eta = eval(bsdf.eta, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real k = eval(bsdf.k, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real filmEta = eval(bsdf.filmEta, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real height = eval(bsdf.height, vertex.uv, vertex.uv_screen_size, texture_pool);

	//Real eta_2 = mix(1.0, filmEta, smoothstep(0.0, 0.03, height));

	Real NDotL = dot(frame.n, dir_out);
	Real NDotV = dot(frame.n, dir_in);
	Vector3 half_vector = normalize(dir_in + dir_out);
	Real NDotH = dot(frame.n, half_vector);

	Real eta_2 = filmEta;
	Real eta_3 = eta;
	Real kappa3 = k;
	Real Dinc = 0.5; // TODO: Different unit? 
	// Real Dinc = height;

	// Force eta_2 -> 1.0 when Dinc -> 0.0
	// Equivalent to: float eta_2 = mix(1.0, eta2, smoothstep(0.0, 0.03, Dinc));
	// https://thebookofshaders.com/glossary/?search=smoothstep
	Real smoothstep;
	smoothstep = std::clamp((Dinc - 0.0) / (0.03 - 0.0), 0.0, 1.0);
    smoothstep = smoothstep * smoothstep * (3.0 - 2.0 * smoothstep);
	eta_2 = 1.0 * (1 - smoothstep) + eta_2 * smoothstep;

	// V=dir_in, L=vir_out
	Real cosTheta1 = dot(half_vector, dir_out);
	Real cosTheta2 = sqrt(1.0 - sqr(1.0 / eta_2) * (1 - sqr(cosTheta1)));

	// First interface
	Vector2 R12, phi12;
	fresnelDielectric(cosTheta1, 1.0, eta_2, R12, phi12);
	// std::cout << R12 << std::endl;
	// Vector2 R21 = R12;
	Vector2 T121 = Vector2(1.0, 1.0) - R12;
	Vector2 phi21 = Vector2(c_PI, c_PI) - phi12;

	// Second interface
	Vector2 R23, phi23;
	fresnelConductor(cosTheta2, eta_2, eta_3, kappa3, R23, phi23);
	// std::cout << R23 << std::endl;

	// Phase shift
	float OPD = Dinc * cosTheta2;
	Vector2 phi2 = phi21 + phi23;

	// Compound terms
	Vector3 I = Vector3(0, 0, 0);
	Vector2 R123 = R12 * R23;
	// std::cout << "R123" << R123.x << R123.y << std::endl;
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Vector2 r123 = Vector2(sqrt(R123.x), sqrt(R123.y));
	r123 = Vector2(0.95, 0.95);
	// std::cout << r123.x << r123.y << std::endl;
	Vector2 Rs = Vector2(sqr(T121).x * R23.x, sqr(T121).y * R23.y) / (1 - R123);

	// Reflectance term for m=0 (DC term amplitude)
	Vector2 C0 = R12 + Rs;
	Vector3 S0 = evalSensitivity_bsdf(0.0, 0.0);
	I += (Real)depol(C0) * S0;

	// Reflectance term for m>0 (pairs of diracs)
	Vector2 Cm = Rs - T121;
	for (int m = 1; m <= 3; ++m) {
		Cm = Cm * r123;
		Vector3 SmS = 2.0 * evalSensitivity_bsdf(m * OPD, m * phi2.x);
		Vector3 SmP = 2.0 * evalSensitivity_bsdf(m * OPD, m * phi2.y);
		I += depolColor(Cm.x * SmS, Cm.y * SmP);
	}

	Real XYZ_TO_RGB[9] = { 2.3706743, -0.5138850, 0.0052982, -0.9000405, 1.4253036, -0.0146949, -0.4706338, 0.0885814, 1.0093968 };

	for (int i = 0; i < 3; ++i) {
		Real x = 0.0;
		for (int j = 0; j < 3; ++j) {
			int index = i + 3 * j;
			x += I[j] * XYZ_TO_RGB[index];
		}
		x = std::clamp(x, 0.0, 1.0);
		I[i] = x;

	}

	// TEST: try use IridescenceTerm()
	// std::cout << "before: " << I << std::endl;
	// I = IridescenceTerm(make_const_spectrum(Dinc), cosTheta1, make_const_spectrum(1.0), make_const_spectrum(eta_2), make_const_spectrum(eta_3), make_const_spectrum(kappa3), wavelengths);
	// std::cout << "after: " << I << std::endl;

	Real D = GTR2(NDotH, alpha);
	if (D == 0) {
		return make_zero_spectrum();
	}
	Real G = smith_masking_gtr2(to_local(frame, dir_in), alpha) * smith_masking_gtr2(to_local(frame, dir_out), alpha);

	return D * G * I / (4.0 * dot(frame.n, dir_in));



}

Real pdf_sample_bsdf_op::operator()(const ThinFilm& bsdf) const {
	if (dot(vertex.geometric_normal, dir_in) < 0 ||
		dot(vertex.geometric_normal, dir_out) < 0) {
		return 0.0;
	}

	Frame frame = vertex.shading_frame;
	if (dot(frame.n, dir_in) < 0) {
		frame = -frame;
	}

	Real alpha = eval(bsdf.alpha, vertex.uv, vertex.uv_screen_size, texture_pool);
	Vector3 half_vector = normalize(dir_in + dir_out);
	Real NDotH = dot(frame.n, half_vector);
	Real NDotV = dot(frame.n, dir_in);

	Real D = GTR2(NDotH, alpha);
	Real G1 = smith_masking_gtr2(to_local(frame, dir_in), alpha);

	return D * G1 / (4.0 * NDotV);
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const ThinFilm& bsdf) const {
	if (dot(vertex.geometric_normal, dir_in) < 0) {
		return {};
	}

	Frame frame = vertex.shading_frame;
	if (dot(frame.n, dir_in) < 0) {
		frame = -frame;
	}


	Real alpha = eval(bsdf.alpha, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real eta = eval(bsdf.eta, vertex.uv, vertex.uv_screen_size, texture_pool);

	Vector3 local_dir_in = to_local(frame, dir_in);
	Vector3 local_micro_normal = sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
	Vector3 half_vector = to_world(frame, local_micro_normal);

	Spectrum reflected = normalize(dir_in + 2 * dot(dir_in, half_vector) * half_vector);
	reflected = 2 * dot(dir_in, half_vector) * half_vector - dir_in;
	return BSDFSampleRecord{ reflected, eta, alpha };
}


TextureSpectrum get_texture_op::operator()(const ThinFilm& bsdf) const {
	return make_constant_spectrum_texture(make_zero_spectrum());
}