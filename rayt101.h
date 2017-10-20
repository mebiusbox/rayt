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

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

namespace rayt {

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
} // namespace rayt
