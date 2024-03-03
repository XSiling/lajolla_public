#include "../microfacet.h"

Spectrum eval_op::operator()(const ThinFilm &bsdf) const {
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
	Spectrum eta = eval(bsdf.eta, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum k = eval(bsdf.k, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum filmEta = eval(bsdf.filmEta, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum height = eval(bsdf.height, vertex.uv, vertex.uv_screen_size, texture_pool);

	//Real eta_2 = mix(1.0, filmEta, smoothstep(0.0, 0.03, height));

	Real NDotL = dot(frame.n, dir_out);
	Real NDotV = dot(frame.n, dir_in);
	Vector3 half_vector = normalize(dir_in + dir_out);
	Real NDotH = dot(frame.n, half_vector);
	
	Real D = GTR2(NDotH, alpha);
	if (D == 0) {
		return make_zero_spectrum();
	}

	Spectrum m_eta1 = make_const_spectrum(1.0);
	Spectrum m_eta = eta;
	Spectrum m_eta2 = filmEta;

	Spectrum I = IridescenceTerm(height, dot(dir_in, half_vector), m_eta1, m_eta2, m_eta, k, wavelengths, true, false);// * (1.0, 1.0, 1.0) specular

	Real G = smith_masking_gtr2(to_local(frame, dir_in), alpha) * smith_masking_gtr2(to_local(frame, dir_out), alpha);
	
	return D * I * G / (4.0 * NDotV);
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

	Real pdf;


	Vector3 wavelengths = { 580.0, 550.0, 450.0 };

	Real alpha = eval(bsdf.alpha, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum eta = eval(bsdf.eta, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum k = eval(bsdf.k, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum filmEta = eval(bsdf.filmEta, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum height = eval(bsdf.height, vertex.uv, vertex.uv_screen_size, texture_pool);

	Spectrum m_eta1 = make_const_spectrum(1.0);
	Spectrum m_eta = eta;
	Spectrum m_eta2 = filmEta;

	Vector3 local_dir_in = to_local(frame, dir_in);
	Vector3 local_micro_normal = sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
	Vector3 half_vector = to_world(frame, local_micro_normal);

	Spectrum reflected = normalize(dir_in + 2 * dot(dir_in, half_vector) * half_vector);
	reflected = 2 * dot(dir_in, half_vector) * half_vector - dir_in;
	return BSDFSampleRecord{ reflected, 1.0, alpha};

}


TextureSpectrum get_texture_op::operator()(const ThinFilm& bsdf) const {
	return make_constant_spectrum_texture(make_zero_spectrum());
}