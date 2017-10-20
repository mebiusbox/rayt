/*
* Copyright (c) 2017 mebiusbox software. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
*/
#include "rayt.h"

namespace rayt {

	class Shape;
	class Material;
	typedef std::shared_ptr<Shape> ShapePtr;
	typedef std::shared_ptr<Material> MaterialPtr;

	//-------------------------------------------------------------------------

	class HitRec {
	public:
		float t;
		float u;
		float v;
		vec3 p;
		vec3 n;
		MaterialPtr mat;
	};

	//-------------------------------------------------------------------------

	class Pdf {
	public:
		virtual float value(const HitRec& hrec, const vec3& direction) const = 0;
		virtual vec3 generate(const HitRec& hrec) const = 0;
	};

	class CosinePdf : public Pdf {
	public:
		CosinePdf() { }

		virtual float value(const HitRec& hrec, const vec3& direction) const override {
			float cosine = dot(normalize(direction), hrec.n);
			if (cosine > 0) {
				return cosine / PI;
			}
			else {
				return 0;
			}
		}

		virtual vec3 generate(const HitRec& hrec) const override {
			ONB uvw; uvw.build_from_w(hrec.n);
			vec3 v = uvw.local(random_cosine_direction());
			return v;
		}
	};

	class MixturePdf : public Pdf {
	public:
		MixturePdf(const Pdf* p0, const Pdf* p1) { m_pdfs[0] = p0; m_pdfs[1] = p1; }

		virtual float value(const HitRec& hrec, const vec3& direction) const override {
			float pdf0 = m_pdfs[0]->value(hrec, direction);
			float pdf1 = m_pdfs[1]->value(hrec, direction);
			return 0.5f*pdf0 + 0.5f*pdf1;
		}

		virtual vec3 generate(const HitRec& hrec) const override {
			if (drand48()<0.5f) {
				return m_pdfs[0]->generate(hrec);
			}
			else {
				return m_pdfs[1]->generate(hrec);
			}
		}

	private:
		const Pdf* m_pdfs[2];
	};

	//-------------------------------------------------------------------------

	class ScatterRec {
	public:
		Ray	ray;
		vec3 albedo;
		const Pdf* pdf;
		bool is_specular;
	};

	class Material {
	public:
		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const = 0;
		virtual float scattering_pdf(const Ray& r, const HitRec& hrec) const { return 0; }
		virtual vec3 emitted(const Ray& r, const HitRec& hrec) const { return vec3(0); }
	};

	class Lambertian : public Material {
	public:
		Lambertian(const TexturePtr& a)
			: m_albedo(a) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			srec.albedo = m_albedo->value(hrec.u, hrec.v, hrec.p);
			srec.pdf = &m_pdf;
			srec.is_specular = false;
			return true;
		};

		virtual float scattering_pdf(const Ray& r, const HitRec& hrec) const override {
			return std::max(dot(hrec.n, normalize(r.direction())), 0.0f) / PI;
		}

	private:
		TexturePtr m_albedo;
		CosinePdf  m_pdf;
	};

	class Metal : public Material {
	public:
		Metal(const TexturePtr& a, float fuzz)
			: m_albedo(a)
			, m_fuzz(fuzz) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			vec3 reflected = reflect(normalize(r.direction()), hrec.n);
			reflected += m_fuzz*random_in_unit_sphere();
			srec.ray = Ray(hrec.p, reflected);
			srec.albedo = m_albedo->value(hrec.u, hrec.v, hrec.p);
			srec.pdf = nullptr;
			srec.is_specular = true;
			return dot(srec.ray.direction(), hrec.n) > 0;
		}

	private:
		TexturePtr m_albedo;
		float m_fuzz;
	};

	class Dielectric : public Material {
	public:
		Dielectric(float ri)
			: m_ri(ri) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			srec.is_specular = true;
			srec.pdf = nullptr;
			vec3 outward_normal;
			vec3 reflected = reflect(r.direction(), hrec.n);
			float ni_over_nt;
			float reflect_prob;
			float cosine;
			if (dot(r.direction(), hrec.n) > 0) {
				outward_normal = -hrec.n;
				ni_over_nt = m_ri;
				cosine = m_ri * dot(r.direction(), hrec.n) / length(r.direction());
			}
			else {
				outward_normal = hrec.n;
				ni_over_nt = recip(m_ri);
				cosine = -dot(r.direction(), hrec.n) / length(r.direction());
			}

			srec.albedo = vec3(1);

			vec3 refracted;
			if (refract(-r.direction(), outward_normal, ni_over_nt, refracted)) {
				reflect_prob = schlick(cosine, m_ri);
			}
			else {
				reflect_prob = 1;
			}

			if (drand48() < reflect_prob) {
				srec.ray = Ray(hrec.p, reflected);
			}
			else {
				srec.ray = Ray(hrec.p, refracted);
			}

			return true;
		}

	private:
		float m_ri;
	};

	class DiffuseLight : public Material {
	public:
		DiffuseLight(const TexturePtr& emit)
			: m_emit(emit) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			return false;
		}

		virtual vec3 emitted(const Ray& r, const HitRec& hrec) const override {
			if (dot(hrec.n, r.direction()) < 0) {
				return m_emit->value(hrec.u, hrec.v, hrec.p);
			}
			else {
				return vec3(0);
			}
		}

	private:
		TexturePtr m_emit;
	};

	//-------------------------------------------------------------------------

	class Shape {
	public:
		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const = 0;
		virtual float pdf_value(const vec3& o, const vec3& v) const { return 0; }
		virtual vec3 random(const vec3& o) const { return vec3(1, 0, 0); }
	};

	class ShapeList : public Shape {
	public:
		ShapeList() {}

		void add(const ShapePtr& shape) {
			m_list.push_back(shape);
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			HitRec temp_rec;
			bool hit_anything = false;
			float closest_so_far = t1;
			for (auto& p : m_list) {
				if (p->hit(r, t0, closest_so_far, temp_rec)) {
					hit_anything = true;
					closest_so_far = temp_rec.t;
					hrec = temp_rec;
				}
			}
			return hit_anything;
		}

	private:
		std::list<ShapePtr> m_list;
	};

	class Sphere : public Shape {
	public:
		Sphere() {}
		Sphere(const vec3& c, float r, const MaterialPtr& mat)
			: m_center(c)
			, m_radius(r)
			, m_material(mat) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			vec3 oc = r.origin() - m_center;
			float a = dot(r.direction(), r.direction());
			float b = 2.0f*dot(oc, r.direction());
			float c = dot(oc, oc) - pow2(m_radius);
			float D = b*b - 4 * a*c;
			if (D > 0) {
				float root = sqrtf(D);
				float temp = (-b - root) / (2.0f*a);
				if (temp < t1 && temp > t0) {
					hrec.t = temp;
					hrec.p = r.at(hrec.t);
					hrec.n = (hrec.p - m_center) / m_radius;
					hrec.mat = m_material;
					get_sphere_uv(hrec.n, hrec.u, hrec.v);
					return true;
				}
				temp = (-b + root) / (2.0f*a);
				if (temp < t1 && temp > t0) {
					hrec.t = temp;
					hrec.p = r.at(hrec.t);
					hrec.n = (hrec.p - m_center) / m_radius;
					hrec.mat = m_material;
					get_sphere_uv(hrec.n, hrec.u, hrec.v);
					return true;
				}
			}

			return false;
		}

	private:
		vec3 m_center;
		float m_radius;
		MaterialPtr m_material;
	};

	class Rect : public Shape {
	public:
		enum AxisType {
			kXY = 0,
			kXZ,
			kYZ
		};
		Rect() {}
		Rect(float x0, float x1, float y0, float y1, float k, AxisType axis, const MaterialPtr& m)
			: m_x0(x0)
			, m_x1(x1)
			, m_y0(y0)
			, m_y1(y1)
			, m_k(k)
			, m_axis(axis)
			, m_material(m) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {

			int xi, yi, zi;
			vec3 axis;
			switch (m_axis) {
			case kXY: xi = 0; yi = 1; zi = 2; axis = vec3::zAxis(); break;
			case kXZ: xi = 0; yi = 2; zi = 1; axis = vec3::yAxis(); break;
			case kYZ: xi = 1; yi = 2; zi = 0; axis = vec3::xAxis(); break;
			}

			float t = (m_k - r.origin()[zi]) / r.direction()[zi];
			if (t < t0 || t > t1) {
				return false;
			}

			float x = r.origin()[xi] + t*r.direction()[xi];
			float y = r.origin()[yi] + t*r.direction()[yi];
			if (x < m_x0 || x > m_x1 || y < m_y0 || y > m_y1) {
				return false;
			}

			hrec.u = (x - m_x0) / (m_x1 - m_x0);
			hrec.v = (y - m_y0) / (m_y1 - m_y0);
			hrec.t = t;
			hrec.mat = m_material;
			hrec.p = r.at(t);
			hrec.n = axis;
			return true;
		}

		virtual float pdf_value(const vec3& o, const vec3& v) const override {
			if (m_axis != kXZ) return 0;
			HitRec hrec;
			if (this->hit(Ray(o, v), 0.001f, FLT_MAX, hrec)) {
				float area = (m_x1 - m_x0) * (m_y1 - m_y0);
				float distance_squared = hrec.t * hrec.t * lengthSqr(v);
				float cosine = fabs(dot(v, hrec.n)) / length(v);
				return distance_squared / (cosine * area);
			}
			else {
				return 0;
			}
		}

		virtual vec3 random(const vec3& o) const override {
			if (m_axis != kXZ) return vec3(1, 0, 0);
			float x = m_x0 + drand48()*(m_x1 - m_x0);
			float y = m_y0 + drand48()*(m_y1 - m_y0);
			vec3 random_point;
			switch (m_axis) {
			case kXY:
				random_point = vec3(x, y, m_k);
				break;
			case kXZ:
				random_point = vec3(x, m_k, y);
				break;
			case kYZ:
				random_point = vec3(m_k, x, y);
				break;
			}
			vec3 v = random_point - o;
			return v;
		}

	private:
		float m_x0, m_x1, m_y0, m_y1, m_k;
		AxisType m_axis;
		MaterialPtr m_material;
	};

	class FlipNormals : public Shape {
	public:
		FlipNormals(const ShapePtr& shape)
			: m_shape(shape) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			if (m_shape->hit(r, t0, t1, hrec)) {
				hrec.n = -hrec.n;
				return true;
			}
			else {
				return false;
			}
		}

	private:
		ShapePtr m_shape;
	};

	class Translate : public Shape {
	public:
		Translate(const ShapePtr& sp, const vec3& displacement)
			: m_shape(sp)
			, m_offset(displacement) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			Ray move_r(r.origin() - m_offset, r.direction());
			if (m_shape->hit(move_r, t0, t1, hrec)) {
				hrec.p += m_offset;
				return true;
			}
			else {
				return false;
			}
		}

	private:
		ShapePtr m_shape;
		vec3 m_offset;
	};

	class Rotate : public Shape {
	public:
		Rotate(const ShapePtr& sp, const vec3& axis, float angle)
			: m_shape(sp)
			, m_quat(Quat::rotation(radians(angle), axis)) {
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			Quat revq = conj(m_quat);
			vec3 origin = rotate(revq, r.origin());
			vec3 direction = rotate(revq, r.direction());
			Ray rot_r(origin, direction);
			if (m_shape->hit(rot_r, t0, t1, hrec)) {
				hrec.p = rotate(m_quat, hrec.p);
				hrec.n = rotate(m_quat, hrec.n);
				return true;
			}
			else {
				return false;
			}
		}

	private:
		ShapePtr m_shape;
		Quat m_quat;
	};

	class Box : public Shape {
	public:
		Box() {}
		Box(const vec3& p0, const vec3& p1, const MaterialPtr& m)
			: m_p0(p0)
			, m_p1(p1)
			, m_list(std::make_unique<ShapeList>()) {

			ShapeList* l = new ShapeList();
			l->add(std::make_shared<Rect>(
				p0.getX(), p1.getX(), p0.getY(), p1.getY(), p1.getZ(), Rect::kXY, m));
			l->add(std::make_shared<FlipNormals>(std::make_shared<Rect>(
				p0.getX(), p1.getX(), p0.getY(), p1.getY(), p0.getZ(), Rect::kXY, m)));
			l->add(std::make_shared<Rect>(
				p0.getX(), p1.getX(), p0.getZ(), p1.getZ(), p1.getY(), Rect::kXZ, m));
			l->add(std::make_shared<FlipNormals>(std::make_shared<Rect>(
				p0.getX(), p1.getX(), p0.getZ(), p1.getZ(), p0.getY(), Rect::kXZ, m)));
			l->add(std::make_shared<Rect>(
				p0.getY(), p1.getY(), p0.getZ(), p1.getZ(), p1.getX(), Rect::kYZ, m));
			l->add(std::make_shared<FlipNormals>(std::make_shared<Rect>(
				p0.getY(), p1.getY(), p0.getZ(), p1.getZ(), p0.getX(), Rect::kYZ, m)));
			m_list.reset(l);
		}

		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const override {
			return m_list->hit(r, t0, t1, hrec);
		}

	private:
		vec3 m_p0, m_p1;
		std::unique_ptr<ShapeList> m_list;
	};

	//-------------------------------------------------------------------------

	class ShapePdf : public Pdf {
	public:
		ShapePdf(const ShapePtr& p, const vec3& origin) : m_ptr(p), m_origin(origin) {
		}

		virtual float value(const HitRec& hrec, const vec3& direction) const override {
			return m_ptr->pdf_value(m_origin, direction);
		}

		virtual vec3 generate(const HitRec& hrec) const override {
			return m_ptr->random(m_origin);
		}

	private:
		ShapePtr m_ptr;
		vec3 m_origin;
	};

	class ShapeBuilder {
	public:
		ShapeBuilder() {}
		ShapeBuilder(const ShapePtr& sp)
			: m_ptr(sp) {
		}

		ShapeBuilder& reset(const ShapePtr& sp) {
			m_ptr = sp;
			return *this;
		}

		ShapeBuilder& sphere(const vec3& c, float r, const MaterialPtr& m) {
			m_ptr = std::make_shared<Sphere>(c, r, m);
			return *this;
		}

		ShapeBuilder& rect(float x0, float x1, float y0, float y1, float k, Rect::AxisType axis, const MaterialPtr& m) {
			m_ptr = std::make_shared<Rect>(x0, x1, y0, y1, k, axis, m);
			return *this;
		}
		ShapeBuilder& rectXY(float x0, float x1, float y0, float y1, float k, const MaterialPtr& m) {
			m_ptr = std::make_shared<Rect>(x0, x1, y0, y1, k, Rect::kXY, m);
			return *this;
		}
		ShapeBuilder& rectXZ(float x0, float x1, float y0, float y1, float k, const MaterialPtr& m) {
			m_ptr = std::make_shared<Rect>(x0, x1, y0, y1, k, Rect::kXZ, m);
			return *this;
		}
		ShapeBuilder& rectYZ(float x0, float x1, float y0, float y1, float k, const MaterialPtr& m) {
			m_ptr = std::make_shared<Rect>(x0, x1, y0, y1, k, Rect::kYZ, m);
			return *this;
		}

		ShapeBuilder& rect(const vec3& p0, const vec3& p1, float k, Rect::AxisType axis, const MaterialPtr& m) {
			switch (axis) {
			case Rect::kXY:
				m_ptr = std::make_shared<Rect>(
					p0.getX(), p1.getX(), p0.getY(), p1.getY(), k, axis, m);
				break;
			case Rect::kXZ:
				m_ptr = std::make_shared<Rect>(
					p0.getX(), p1.getX(), p0.getZ(), p1.getZ(), k, axis, m);
				break;
			case Rect::kYZ:
				m_ptr = std::make_shared<Rect>(
					p0.getY(), p1.getY(), p0.getZ(), p1.getZ(), k, axis, m);
				break;
			}
			return *this;
		}
		ShapeBuilder& rectXY(const vec3& p0, const vec3& p1, float k, const MaterialPtr& m) {
			return rect(p0, p1, k, Rect::kXY, m);
		}
		ShapeBuilder& rectXZ(const vec3& p0, const vec3& p1, float k, const MaterialPtr& m) {
			return rect(p0, p1, k, Rect::kXZ, m);
		}
		ShapeBuilder& rectYZ(const vec3& p0, const vec3& p1, float k, const MaterialPtr& m) {
			return rect(p0, p1, k, Rect::kYZ, m);
		}

		ShapeBuilder& box(const vec3& p0, const vec3& p1, const MaterialPtr& m) {
			m_ptr = std::make_shared<Box>(p0, p1, m);
			return *this;
		}

		ShapeBuilder& flip() {
			m_ptr = std::make_shared<FlipNormals>(m_ptr);
			return *this;
		}

		ShapeBuilder& translate(const vec3& t) {
			m_ptr = std::make_shared<Translate>(m_ptr, t);
			return *this;
		}

		ShapeBuilder& rotate(const vec3& axis, float angle) {
			m_ptr = std::make_shared<Rotate>(m_ptr, axis, angle);
			return *this;
		}

		const ShapePtr& get() const { return m_ptr; }

	private:
		ShapePtr m_ptr;
	};

	//-------------------------------------------------------------------------
	class Scene {
	public:
		Scene(int width, int height, int samples)
			: m_image(std::make_unique<Image>(width, height))
			, m_backColor(0.1f)
			, m_samples(samples) {
		}

		void build() {

			m_backColor = vec3(0);

			// Camera

			vec3 lookfrom(278, 278, -800);
			vec3 lookat(278, 278, 0);
			vec3 vup(0, 1, 0);
			float aspect = float(m_image->width()) / float(m_image->height());
			m_camera = std::make_unique<Camera>(lookfrom, lookat, vup, 40, aspect);

			// Shapes

			MaterialPtr red = std::make_shared<Lambertian>(
				std::make_shared<ColorTexture>(vec3(0.65f, 0.05f, 0.05f)));
			MaterialPtr white = std::make_shared<Lambertian>(
				std::make_shared<ColorTexture>(vec3(0.73f)));
			MaterialPtr green = std::make_shared<Lambertian>(
				std::make_shared<ColorTexture>(vec3(0.12f, 0.45f, 0.15f)));
			MaterialPtr light = std::make_shared<DiffuseLight>(
				std::make_shared<ColorTexture>(vec3(15.0f)));
			MaterialPtr aluminum = std::make_shared<Metal>(
				std::make_shared<ColorTexture>(vec3(0.8f, 0.85f, 0.88f)), 0.0f);

			ShapeList* world = new ShapeList();
			ShapeBuilder builder;
			world->add(builder.rectYZ(0, 555, 0, 555, 555, green).flip().get());
			world->add(builder.rectYZ(0, 555, 0, 555, 0, red).get());
			world->add(builder.rectXZ(213, 343, 227, 332, 554, light).flip().get());
			world->add(builder.rectXZ(0, 555, 0, 555, 555, white).flip().get());
			world->add(builder.rectXZ(0, 555, 0, 555, 0, white).get());
			world->add(builder.rectXY(0, 555, 0, 555, 555, white).flip().get());
			world->add(builder.box(vec3(0), vec3(165), white)
				.rotate(vec3::yAxis(), -18)
				.translate(vec3(130, 0, 65))
				.get());
			world->add(builder.box(vec3(0), vec3(165, 330, 165), aluminum)
				.rotate(vec3::yAxis(), 15)
				.translate(vec3(265, 0, 295))
				.get());
			m_world.reset(world);
		}

		vec3 color(const rayt::Ray& r, const Shape* world, int depth) {
			HitRec hrec;
			if (world->hit(r, 0.001f, FLT_MAX, hrec)) {
				vec3 emitted = hrec.mat->emitted(r, hrec);
				ScatterRec srec;
				if (depth < MAX_DEPTH && hrec.mat->scatter(r, hrec, srec)) {
					if (srec.is_specular) {
						return emitted + mulPerElem(srec.albedo, color(srec.ray, world, depth + 1));
					}
					else {
						ShapePtr p(std::make_shared<Rect>(
							213, 343, 227, 332, 554, Rect::kXZ, MaterialPtr()));
						ShapePdf shapePdf(p, hrec.p);
						MixturePdf mixPdf(&shapePdf, srec.pdf);
						srec.ray = Ray(hrec.p, mixPdf.generate(hrec));
						float pdf_value = mixPdf.value(hrec, srec.ray.direction());
						if (pdf_value > 0) {
							float spdf_value = hrec.mat->scattering_pdf(srec.ray, hrec);
							vec3 albedo = srec.albedo * spdf_value;
							return emitted + mulPerElem(albedo, color(srec.ray, world, depth + 1)) / pdf_value;
						}
						else {
							return emitted;
						}
					}
				}
				else {
					return emitted;
				}
			}
			return background(r.direction());
		}

		vec3 background(const vec3& d) const {
			return m_backColor;
		}

		vec3 backgroundSky(const vec3& d) const {
			vec3 v = normalize(d);
			float t = 0.5f * (v.getY() + 1.0f);
			return lerp(t, vec3(1), vec3(0.5f, 0.7f, 1.0f));
		}

		void render() {

			build();

			int nx = m_image->width();
			int ny = m_image->height();
#pragma omp parallel for schedule(dynamic,1) num_threads(NUM_THREAD)
			for (int j = 0; j<ny; j++) {
				std::cerr << "Rendering (y = " << j << ") " << (100.0 * j / (ny - 1)) << "%" << std::endl;
				for (int i = 0; i<nx; ++i) {
					vec3 c(0);
					for (int s = 0; s<m_samples; ++s) {
						float u = (float(i) + drand48()) / float(nx);
						float v = (float(j) + drand48()) / float(ny);
						Ray r = m_camera->getRay(u, v);
						c += color(r, m_world.get(), 0);
					}
					c /= m_samples;
					m_image->write(i, (ny - j - 1), c.getX(), c.getY(), c.getZ());
				}
			}

			stbi_write_bmp("render.bmp", nx, ny, sizeof(Image::rgb), m_image->pixels());
		}

	private:
		std::unique_ptr<Camera> m_camera;
		std::unique_ptr<Image> m_image;
		std::unique_ptr<Shape> m_world;
		vec3 m_backColor;
		int m_samples;
	};

} // namespace rayt

int main()
{
	int nx = 200;
	int ny = 200;
	int ns = 500;
	std::unique_ptr<rayt::Scene> scene(std::make_unique<rayt::Scene>(nx, ny, ns));
	scene->render();
	return 0;
}

