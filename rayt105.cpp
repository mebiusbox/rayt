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

	//-------------------------------------------------------------------------
	class Scene {
	public:
		Scene(int width, int height)
			: m_image(std::make_unique<Image>(width, height))
			, m_backColor(0.2f) {
		}

		void build() {

			// Camera

			vec3 w(-2.0f, -1.0f, -1.0f);
			vec3 u(4.0f, 0.0f, 0.0f);
			vec3 v(0.0f, 2.0f, 0.0f);
			m_camera = std::make_unique<Camera>(u, v, w);
		}

		float hit_sphere(const vec3& center, float radius, const rayt::Ray& r) const {
			vec3 oc = r.origin() - center;
			float a = dot(r.direction(), r.direction());
			float b = 2.0f * dot(r.direction(), oc);
			float c = dot(oc, oc) - pow2(radius);
			float D = b*b - 4 * a*c;
			if (D < 0) {
				return -1.0f;
			}
			else {
				return (-b - sqrtf(D)) / (2.0f*a);
			}
		}

		vec3 color(const rayt::Ray& r) {
			vec3 c(0, 0, -1);
			float t = hit_sphere(c, 0.5f, r);
			if (t > 0.0f) {
				vec3 N = normalize(r.at(t) - c);
				return 0.5f*(N + vec3(1.0f));
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
					float u = float(i) / float(nx);
					float v = float(j) / float(ny);
					Ray r = m_camera->getRay(u, v);
					vec3 c = color(r);
					m_image->write(i, (ny - j - 1), c.getX(), c.getY(), c.getZ());
				}
			}

			stbi_write_bmp("render.bmp", nx, ny, sizeof(Image::rgb), m_image->pixels());
		}

	private:
		std::unique_ptr<Camera> m_camera;
		std::unique_ptr<Image> m_image;
		vec3 m_backColor;
	};

} // namespace rayt

int main()
{
	int nx = 200;
	int ny = 100;
	std::unique_ptr<rayt::Scene> scene(std::make_unique<rayt::Scene>(nx, ny));
	scene->render();
	return 0;
}

