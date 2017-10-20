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
		vec3 p;
		vec3 n;
		MaterialPtr mat;
	};

	//-------------------------------------------------------------------------

	class ScatterRec {
	public:
		Ray	ray;
		vec3 albedo;
	};

	class Material {
	public:
		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const = 0;
	};

	class Lambertian : public Material {
	public:
		Lambertian(const vec3& c)
			: m_albedo(c) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			vec3 target = hrec.p + hrec.n + random_in_unit_sphere();
			srec.ray = Ray(hrec.p, target - hrec.p);
			srec.albedo = m_albedo;
			return true;
		};

	private:
		vec3 m_albedo;
	};

	class Metal : public Material {
	public:
		Metal(const vec3& c, float fuzz)
			: m_albedo(c)
			, m_fuzz(fuzz) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			vec3 reflected = reflect(normalize(r.direction()), hrec.n);
			reflected += m_fuzz*random_in_unit_sphere();
			srec.ray = Ray(hrec.p, reflected);
			srec.albedo = m_albedo;
			return dot(srec.ray.direction(), hrec.n) > 0;
		}

	private:
		vec3 m_albedo;
		float m_fuzz;
	};

	class Dielectric : public Material {
	public:
		Dielectric(float ri)
			: m_ri(ri) {
		}

		virtual bool scatter(const Ray& r, const HitRec& hrec, ScatterRec& srec) const override {
			vec3 outward_normal;
			vec3 reflected = reflect(r.direction(), hrec.n);
			float ni_over_nt;
			if (dot(r.direction(), hrec.n) > 0) {
				outward_normal = -hrec.n;
				ni_over_nt = m_ri;
			}
			else {
				outward_normal = hrec.n;
				ni_over_nt = recip(m_ri);
			}

			srec.albedo = vec3(1);

			vec3 refracted;
			if (refract(-r.direction(), outward_normal, ni_over_nt, refracted)) {
				srec.ray = Ray(hrec.p, refracted);
			}
			else {
				srec.ray = Ray(hrec.p, reflected);
				return false;
			}

			return true;
		}

	private:
		float m_ri;
	};

	//-------------------------------------------------------------------------

	class Shape {
	public:
		virtual bool hit(const Ray& r, float t0, float t1, HitRec& hrec) const = 0;
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
		std::vector<ShapePtr> m_list;
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
					return true;
				}
				temp = (-b + root) / (2.0f*a);
				if (temp < t1 && temp > t0) {
					hrec.t = temp;
					hrec.p = r.at(hrec.t);
					hrec.n = (hrec.p - m_center) / m_radius;
					hrec.mat = m_material;
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

	//-------------------------------------------------------------------------
	class Scene {
	public:
		Scene(int width, int height, int samples)
			: m_image(std::make_unique<Image>(width, height))
			, m_backColor(0.2f)
			, m_samples(samples) {
		}

		void build() {

			// Camera

			vec3 w(-2.0f, -1.0f, -1.0f);
			vec3 u(4.0f, 0.0f, 0.0f);
			vec3 v(0.0f, 2.0f, 0.0f);
			m_camera = std::make_unique<Camera>(u, v, w);

			ShapeList* world = new ShapeList();
			world->add(std::make_shared<Sphere>(
				vec3(0.6, 0, -1), 0.5f,
				std::make_shared<Lambertian>(vec3(0.1f, 0.2f, 0.5f))));
			world->add(std::make_shared<Sphere>(
				vec3(-0.6, 0, -1), 0.5f,
				std::make_shared<Dielectric>(1.5f)));
			world->add(std::make_shared<Sphere>(
				vec3(0, -0.35, -0.8f), 0.15f,
				std::make_shared<Metal>(vec3(0.8f, 0.8f, 0.8f), 0.2f)));
			world->add(std::make_shared<Sphere>(
				vec3(0, -100.5, -1), 100,
				std::make_shared<Lambertian>(vec3(0.8f, 0.8f, 0.0f))));
			m_world.reset(world);
		}

		vec3 color(const rayt::Ray& r, const Shape* world, int depth) {
			HitRec hrec;
			if (depth < MAX_DEPTH && world->hit(r, 0.001f, FLT_MAX, hrec)) {
				ScatterRec srec;
				if (hrec.mat->scatter(r, hrec, srec)) {
					return mulPerElem(srec.albedo, color(srec.ray, world, depth + 1));
				}
				else {
					return vec3(0);
				}
			}
			return backgroundSky(r.direction());
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
#pragma omp parallel for schedule(dynamic, 1) num_threads(NUM_THREAD)
			for (int j = 0; j<ny; ++j) {
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
	int ny = 100;
	int ns = 100;
	std::unique_ptr<rayt::Scene> scene(std::make_unique<rayt::Scene>(nx, ny, ns));
	scene->render();
	return 0;
}

