#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cmath>
#include "CImg.h"
using namespace std;
using namespace cimg_library;
const CImg<unsigned char> readBMP(const char* filename)
{
	FILE* f = fopen(filename, "rb");
	unsigned char info[54];
	fread(info, sizeof(unsigned char), 54, f); // 读入文件头
	int data_offset = *(int*)(&info[0x0A]);
	if (data_offset > 54) {
		fseek(f, (long int)(data_offset - 54), SEEK_CUR);
	}
	// 从文件头读出长和宽
	int width = *(int*)&info[18];
	width = width + width % 2;
	int height = *(int*)&info[22];
	int size = 3 * width * height;
	unsigned char* data = new unsigned char[size]; // 每个像素3个字节
	fread(data, sizeof(unsigned char), size, f); // 读完其他所有字节
	fclose(f);

	/*for (int i = 0; i < size; i += 3)
	{
		unsigned char tmp = data[i];
		data[i] = data[i + 2];
		data[i + 2] = tmp;
	}*/
	CImg<unsigned char> image(width, height, 1, 3, 0);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++) 
		{
			int current = 3 * ((height - i) * width + j);
			image(j, i, 0) = data[current + 2];
			image(j, i, 1) = data[current + 1];
			image(j, i, 2) = data[current];
		}
	}
	return image;
}
int main(int argc, char * argv[])
{
	CImg<unsigned char> image = readBMP(argv[1]);
	int height = image.height();
	int width = image.width();
	double rgb2yuv[3][3] = { {0.299, 0.587, 0.114}, {-0.147, -0.289, 0.436}, {0.615, -0.515, -0.1} };
	double rgb2yiq[3][3] = { {0.299, 0.587, 0.114}, {0.596, -0.274, -0.322}, {0.211, -0.523, 0.312} };
	double rgb2ycbcr[3][3] = { {0.299, 0.587, 0.114}, {0.500, -0.4187, -0.0813}, {-0.1687, -0.3313, 0.5} };
	double rgb2xyz[3][3] = { {0.49, 0.31, 0.2}, {0.177, 0.813, 0.011}, {0, 0.01, 0.99} };

	CImg<unsigned char> image_rgb_r(width, height, 1, 1, 0);
	CImg<unsigned char> image_rgb_g(width, height, 1, 1, 0);
	CImg<unsigned char> image_rgb_b(width, height, 1, 1, 0);

	CImg<double> image_yiq_y(width, height, 1, 1, 0);
	CImg<double> image_yiq_i(width, height, 1, 1, 0);
	CImg<double> image_yiq_q(width, height, 1, 1, 0);

	CImg<double> image_hsi_h(width, height, 1, 1, 0);
	CImg<double> image_hsi_s(width, height, 1, 1, 0);
	CImg<double> image_hsi_i(width, height, 1, 1, 0);

	CImg<double> image_ycbcr_y(width, height, 1, 1, 0);
	CImg<double> image_ycbcr_cb(width, height, 1, 1, 0);
	CImg<double> image_ycbcr_cr(width, height, 1, 1, 0);

	CImg<double> image_xyz_x(width, height, 1, 1, 0);
	CImg<double> image_xyz_y(width, height, 1, 1, 0);
	CImg<double> image_xyz_z(width, height, 1, 1, 0);

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int r = image(j, i, 0);
			int g = image(j, i, 1);
			int b = image(j, i, 2);
			image_rgb_r(j, i) = r;
			image_rgb_g(j, i) = g;
			image_rgb_b(j, i) = b;
		}
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int r = image(j, i, 0);
			int g = image(j, i, 1);
			int b = image(j, i, 2);
			image_yiq_y(j, i) = rgb2yiq[0][0] * r + rgb2yiq[0][1] * g + rgb2yiq[0][2] * b;
			image_yiq_i(j, i) = rgb2yiq[1][0] * r + rgb2yiq[1][1] * g + rgb2yiq[1][2] * b;
			image_yiq_q(j, i) = rgb2yiq[2][0] * r + rgb2yiq[2][1] * g + rgb2yiq[2][2] * b;
		}
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int r = image(j, i, 0);
			int g = image(j, i, 1);
			int b = image(j, i, 2);
			image_hsi_i(j, i) = (r + g + b)/3;
			image_hsi_s(j, i) = 1 - min(min(r, g),b) / image_hsi_i(j, i);
			image_hsi_h(j, i) = acos((r-g+r-b)/2/sqrt(pow(r-g,2)+(r-b)*(g-b)));
		}
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int r = image(j, i, 0);
			int g = image(j, i, 1);
			int b = image(j, i, 2);
			image_ycbcr_y(j, i) = rgb2ycbcr[0][0] * r + rgb2ycbcr[0][1] * g + rgb2ycbcr[0][2] * b;
			image_ycbcr_cb(j, i) = rgb2ycbcr[1][0] * r + rgb2ycbcr[1][1] * g + rgb2ycbcr[1][2] * b;
			image_ycbcr_cr(j, i) = rgb2ycbcr[2][0] * r + rgb2ycbcr[2][1] * g + rgb2ycbcr[2][2] * b;
		}
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int r = image(j, i, 0);
			int g = image(j, i, 1);
			int b = image(j, i, 2);
			image_xyz_x(j, i) = rgb2xyz[0][0] * r + rgb2xyz[0][1] * g + rgb2xyz[0][2] * b;
			image_xyz_y(j, i) = rgb2xyz[1][0] * r + rgb2xyz[1][1] * g + rgb2xyz[1][2] * b;
			image_xyz_z(j, i) = rgb2xyz[2][0] * r + rgb2xyz[2][1] * g + rgb2xyz[2][2] * b;
		}
	}
	image.display("bmp");

	image_rgb_r.display("rgb_r");
	image_rgb_g.display("rgb_g");
	image_rgb_b.display("rgb_b");

	image_yiq_y.display("yiq_y");
	image_yiq_i.display("yiq_i");
	image_yiq_q.display("yiq_q");

	CImg<unsigned char>(image_hsi_h * 255).display("hsi_h");
	CImg<unsigned char>(image_hsi_s * 255).display("hsi_s");
	image_hsi_i.display("hsi_i");

	image_ycbcr_y.display("ycbcr_y");
	image_ycbcr_cb.display("ycbcr_cb");
	image_ycbcr_cr.display("ycbcr_cr");

	image_xyz_x.display("xyz: x");
	image_xyz_y.display("xyz: y");
	image_xyz_z.display("xyz: z");

	cout <<"pixel in(50,50)"<<"r:"<< (int)image(50, 50, 0) <<"  g:"<< (int)image(50, 50, 1) <<"  b:"<< (int)image(50, 50, 2) << endl;
}