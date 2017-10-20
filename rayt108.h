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
#include <memory>
#include <iostream>
#include <vector>

#define NUM_THREAD 8

#include <float.h>	// FLT_MIN, FLT_MAX
#define PI 3.14159265359f
#define PI2 6.28318530718f
#define RECIP_PI 0.31830988618f
#define RECIP_PI2 0.15915494f
#define LOG2 1.442695f
#define EPSILON 1e-6f
#define GAMMA_FACTOR 2.2f

// https://github.com/kikikikina/drand48/blob/master/main.cpp
#include <random>
inline float drand48() {
	return float(((double)(rand()) / (RAND_MAX))); /* RAND_MAX = 32767 */
}

inline float pow2(float x) { return x*x; }
inline float pow3(float x) { return x*x*x; }
inline float pow4(float x) { return x*x*x*x; }
inline float pow5(float x) { return x*x*x*x*x; }
inline float clamp(float x, float a, float b) { return x < a ? a : x > b ? b : x; }
inline float saturate(float x) { return x < 0.f ? 0.f : x > 1.f ? 1.f : x; }
inline float recip(float x) { return 1.f / x; }
inline float mix(float a, float b, float t) { return a*(1.f - t) + b*t; /* return a + (b-a) * t; */ }
inline float step(float edge, float x) { return (x < edge) ? 0.f : 1.f; }
inline float smoothstep(float a, float b, float t) { if (a >= b) return 0.f; float x = saturate((t - a) / (b - a)); return x*x*(3.f - 2.f * t); }
inline float radians(float deg) { return (deg / 180.f)*PI; }
inline float degrees(float rad) { return (rad / PI) * 180.f; }
inline float nearyeq(float a, float b, float eps = EPSILON) { return fabsf(a - b) <= eps; }
inline float iszero(float a) { return fabsf(a) <= EPSILON; };
inline float safe_recip(float x) { return 1.f / (x + EPSILON); }
//inline float safe_recip(float x) { return 1.0f / (iszero(x)?EPSILON:x); }

// https://github.com/erwincoumans/sce_vectormath
//#define _VECTORMATH_DEBUG
#include <vectormath/scalar/cpp/vectormath_aos.h>
using namespace Vectormath::Aos;
typedef Vector3 vec3;
typedef Vector3 col3;

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

namespace rayt {

	//-------------------------------------------------------------------------

	inline vec3 random_vector() {
		return vec3(drand48(), drand48(), drand48());
	}

	inline vec3 random_in_unit_sphere() {
		vec3 p;
		do {
			p = 2.f * random_vector() - vec3(1.f);
		} while (lengthSqr(p) >= 1.f);
		return p;
	}

	//-------------------------------------------------------------------------

	class Image {
	public:
		struct rgb {
			unsigned char r;
			unsigned char g;
			unsigned char b;
		};

		Image() : m_pixels(nullptr) { }
		Image(int w, int h) {
			m_width = w;
			m_height = h;
			m_pixels.reset(new rgb[m_width*m_height]);
		}

		int width() const { return m_width; }
		int height() const { return m_height; }
		void* pixels() const { return m_pixels.get(); }

		void write(int x, int y, float r, float g, float b) {
			int index = m_width*y + x;
			m_pixels[index].r = static_cast<unsigned char>(r*255.99f);
			m_pixels[index].g = static_cast<unsigned char>(g*255.99f);
			m_pixels[index].b = static_cast<unsigned char>(b*255.99f);
		}

	private:
		int m_width;
		int m_height;
		std::unique_ptr<rgb[]> m_pixels;
	};

	//-------------------------------------------------------------------------

	class Ray {
	public:
		Ray() {}
		Ray(const vec3& o, const vec3& dir)
			: m_origin(o)
			, m_direction(dir) {
		}

		const vec3& origin() const { return m_origin; }
		const vec3& direction() const { return m_direction; }
		vec3 at(float t) const { return m_origin + t*m_direction; }

	private:
		vec3 m_origin;
		vec3 m_direction;
	};

	//-------------------------------------------------------------------------

	class Camera {
	public:
		Camera() {}
		Camera(const vec3& u, const vec3& v, const vec3& w) {
			m_origin = vec3(0);
			m_uvw[0] = u;
			m_uvw[1] = v;
			m_uvw[2] = w;
		}
		Camera(const vec3& lookfrom, const vec3& lookat, const vec3& vup, float vfov, float aspect) {
			vec3 u, v, w;
			float halfH = tanf(radians(vfov) / 2.0f);
			float halfW = aspect * halfH;
			m_origin = lookfrom;
			w = normalize(lookfrom - lookat);
			u = normalize(cross(vup, w));
			v = cross(w, u);
			m_uvw[2] = m_origin - halfW*u - halfH*v - w;
			m_uvw[0] = 2.0f * halfW * u;
			m_uvw[1] = 2.0f * halfH * v;
		}

		Ray getRay(float u, float v) const {
			return Ray(m_origin, m_uvw[2] + m_uvw[0] * u + m_uvw[1] * v - m_origin);
		}

	private:
		vec3 m_origin;
		vec3 m_uvw[3];
	};
} // namespace rayt
