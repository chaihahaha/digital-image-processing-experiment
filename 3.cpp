#include <iostream>
//#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iterator>
#include <random>
#include <algorithm>
#include "CImg.h"
#define sizexsize 90
using namespace std;
using namespace cimg_library;
CImg<unsigned char> readBMP(char* filename)
{
	FILE* f;
	errno_t err;
	err = fopen_s(&f, filename, "rb");
	/*assert(err == 0);*/
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
CImg<int> median_filter(CImg<int> image_salt_pepper_noise, int size)
{
	const int width = image_salt_pepper_noise.width();
	const int height = image_salt_pepper_noise.height();
	CImg<int> image_median_filted(width, height, 1, 1, 0);
	image_median_filted = image_salt_pepper_noise;
	int subimage_arr[sizexsize];
	for (int i = size / 2; i < height - size / 2; i++)
	{
		for (int j = size / 2; j < width - size / 2; j++)
		{
			for (int i1 = -size / 2; i1 <= size / 2; i1++)
			{
				for (int j1 = -size / 2; j1 <= size / 2; j1++)
				{
					subimage_arr[size * (i1 + size / 2) + j1 + size / 2] = image_salt_pepper_noise(j + j1, i + i1);
				}
			}
			sort(subimage_arr, subimage_arr + size * size);
			image_median_filted(j, i) = subimage_arr[size * size / 2 + 1];
		}
	}
	return image_median_filted;
}
CImg<int> mean_filter(CImg<int> image_salt_pepper_noise, int size)
{
	const int width = image_salt_pepper_noise.width();
	const int height = image_salt_pepper_noise.height();
	CImg<int> image_mean_filted(width, height, 1, 1, 0);
	image_mean_filted = image_salt_pepper_noise;
	CImg<int> subimage(size, size, 1, 1, 0);
	for (int i = size / 2; i < height - size / 2; i++)
	{
		for (int j = size / 2; j < width - size / 2; j++)
		{
			for (int i1 = -size / 2; i1 <= size / 2; i1++)
			{
				for (int j1 = -size / 2; j1 <= size / 2; j1++)
				{
					subimage(j1 + size / 2, i1 + size / 2) = image_salt_pepper_noise(j + j1, i + i1);
				}
			}
			image_mean_filted(j, i) = (int)subimage.mean();
		}
	}
	return image_mean_filted;
}
CImg<int> operator_filter(CImg<int> image_salt_pepper_noise, int op[3][3])
{
	const int width = image_salt_pepper_noise.width();
	const int height = image_salt_pepper_noise.height();
	CImg<int> image_operator_filted(width, height, 1, 1, 0);
	image_operator_filted = image_salt_pepper_noise;
	for (int i = 1; i < height - 1; i++)
	{
		for (int j = 1; j < width - 1; j++)
		{
			int summation = 0;
			for (int i1 = -1; i1 <= 1; i1++)
			{
				for (int j1 = -1; j1 <= 1; j1++)
				{
					summation += image_salt_pepper_noise(j + j1, i + i1) * op[i1 + 1][j1 + 1];
				}
			}
			image_operator_filted(j, i) = summation;
		}
	}
	return image_operator_filted;
}
CImg<double> DWT_1d(CImg<double> input)
{
	int width = input.width();
	if (input.height() != 1)
		cout << "fuck height>1" << endl;
	CImg<double> result(width, 1, 1, 1, 0);
	for (int i = 0; i < width; i += 2)
	{
		result(i / 2, 0) = (input(i, 0) + input(i + 1, 0)) / 2.0;
		result(i / 2 + width / 2, 0) = input(i, 0) - result(i / 2, 0);
	}
	return result;
}

CImg<double> DWT_2d_horizontal(CImg<double> input)
{
	int width = input.width();
	int height = input.height();
	CImg<double> tmp_row(width, 1);
	CImg<double> result(width, height, 1, 1, 0);
	for (int i = 0; i < height; i++)
	{
		tmp_row = DWT_1d(input.get_crop(0, i, width - 1, i));
		for (int k = 0; k < width; k++)
		{
			result(k, i) = tmp_row(k, 0);
		}
		//row
	}
	return result;
}

CImg<double> DWT_2d(CImg<double> input)
{
	int width = input.width();
	int height = input.height();
	CImg<double> result(width, height, 1, 1, 0);
	result = DWT_2d_horizontal(input);
	result = DWT_2d_horizontal(result.get_transpose()).get_transpose();
	return result;
}
CImg<double> DWT(CImg<int> input, int level)
{
	int height = input.height();
	int width = input.width();
	CImg<double> result(width, height, 1, 1, 0);
	if (level == 0)
		return DWT_2d(input);
	else
	{
		result = DWT_2d(input);
		CImg<double> tmp(width / 2, height / 2, 1, 1, 0);
		tmp = DWT(result.get_crop(0, 0, width / 2 - 1, height / 2 - 1), level - 1);
		for (int i = 0; i < height / 2; i++)
		{
			for (int j = 0; j < width / 2; j++)
			{
				result(j, i) = tmp(j, i);
			}
		}
		return result;
	}
}

CImg<double> IDWT_1d(CImg<double> input)
{
	int width = input.width();
	if (input.height() != 1)
		cout << "fuck height>1" << endl;
	CImg<double> result(width, 1, 1, 1, 0);
	for (int i = 0; i < width; i += 2)
	{
		result(i) = input(i / 2, 0) + input(i / 2 + width / 2, 0);
		result(i + 1) = input(i / 2, 0) - input(i / 2 + width / 2, 0);
	}
	return result;
}

CImg<double> IDWT_2d_horizontal(CImg<double> input)
{
	int width = input.width();
	int height = input.height();
	CImg<double> tmp_row(width, 1);
	CImg<double> result(width, height, 1, 1, 0);
	for (int i = 0; i < height; i++)
	{
		tmp_row = IDWT_1d(input.get_crop(0, i, width - 1, i));
		for (int k = 0; k < width; k++)
		{
			result(k, i) = tmp_row(k, 0);
		}
		//row
	}
	return result;
}

CImg<double> IDWT_2d(CImg<double> input)
{
	int width = input.width();
	int height = input.height();
	CImg<double> result(width, height, 1, 1, 0);
	result = IDWT_2d_horizontal(input);
	result = IDWT_2d_horizontal(result.get_transpose()).get_transpose();
	return result;
}
CImg<double> IDWT(CImg<int> input, int level)
{
	int height = input.height();
	int width = input.width();
	CImg<double> result(width, height, 1, 1, 0);
	if (level == 0)
		return IDWT_2d(input);
	else
	{
		result = input;
		CImg<double> tmp(width / 2, height / 2, 1, 1, 0);
		tmp = IDWT(result.get_crop(0, 0, width / 2 - 1, height / 2 - 1), level - 1);
		for (int i = 0; i < height / 2; i++)
		{
			for (int j = 0; j < width / 2; j++)
			{
				result(j, i) = tmp(j, i);
			}
		}
		return IDWT_2d(result);
	}
}
CImg<double> wavelet_filter(CImg<double> image_salt_pepper_noise)
{
	int width = image_salt_pepper_noise.width();
	int height = image_salt_pepper_noise.height();
	CImg<double> image_dwt(width, height, 1, 1, 0);
	image_dwt = DWT(image_salt_pepper_noise, 3);
	image_dwt.display("离散哈尔小波变换后");
	image_dwt.save("离散哈尔小波变换后.bmp");
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if (abs(image_dwt(j, i)) < 8)
				image_dwt(j, i) = 0;
		}
	}
	image_dwt.display("离散哈尔小波变换滤波后");
	image_dwt.save("离散哈尔小波变换滤波后.bmp");
	CImg<double> image_idwt(width, height, 1, 1, 0);
	image_idwt = IDWT(image_dwt, 3);
	image_idwt.display("离散哈尔小波逆变换后");
	image_idwt.save("离散哈尔小波逆变换后.bmp");
	return image_idwt;
}
double SSIM(CImg<double> ima, CImg<double> imb, float block_size)
{
	if (ima.width() != imb.width() || ima.height() != imb.height())
		return -1;
	int width = ima.width();
	int height = ima.height();
	CImg<double> ua(width, height, 1, 1, 0), ub(width, height, 1, 1, 0), ua_sq(width, height, 1, 1, 0), ub_sq(width, height, 1, 1, 0), ua_ub(width, height, 1, 1, 0), siga_sq(width, height, 1, 1, 0), sigb_sq, sigab(width, height, 1, 1, 0), ssim_map(width, height, 1, 1, 0);
	double K1 = 0.01;
	double K2 = 0.03;
	double L = 255;
	double C1 = pow((K1*L), 2);
	double C2 = pow((K2*L), 2);
	ua = ima.blur_box(block_size).normalize(0, 255);
	ub = imb.blur_box(block_size).normalize(0, 255);
	ua_sq = ua.mul(ua);
	ub_sq = ub.mul(ub);
	ua_ub = ua.mul(ub);
	siga_sq = ima.pow(2).blur_box(block_size).normalize(0, 255) - ua_sq;
	sigb_sq = imb.pow(2).blur_box(block_size).normalize(0, 255) - ub_sq;
	sigab = (ima.mul(imb)).blur_box(block_size).normalize(0, 255) - ua_ub;
	ssim_map = ((2 * ua_ub + C1).mul(2 * sigab + C2)).div((ua_sq + ub_sq + C1).mul(siga_sq + sigb_sq + C2));
	return ssim_map.mean();
}
int main(int argc, char * argv[])
{
	// 读取文件，获得原图像
	CImg<unsigned char> image = readBMP(argv[1]);
	int height = image.height();
	int width = image.width();
	double rgb2ycbcr[3][3] = { {0.299, 0.587, 0.114}, {0.500, -0.4187, -0.0813}, {-0.1687, -0.3313, 0.5} };

	CImg<int> image_ycbcr_y(width, height, 1, 1, 0);
	CImg<int> image_ycbcr_cb(width, height, 1, 1, 0);
	CImg<int> image_ycbcr_cr(width, height, 1, 1, 0);

	// 转换到ycbcr颜色空间
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int r = image(j, i, 0);
			int g = image(j, i, 1);
			int b = image(j, i, 2);
			image_ycbcr_y(j, i) = (int)(rgb2ycbcr[0][0] * r + rgb2ycbcr[0][1] * g + rgb2ycbcr[0][2] * b);
			image_ycbcr_cb(j, i) = (int)(rgb2ycbcr[1][0] * r + rgb2ycbcr[1][1] * g + rgb2ycbcr[1][2] * b);
			image_ycbcr_cr(j, i) = (int)(rgb2ycbcr[2][0] * r + rgb2ycbcr[2][1] * g + rgb2ycbcr[2][2] * b);
		}
	}
	// 显示ycbcr图
	image_ycbcr_y.display("ycbcr_y");
	image_ycbcr_y.save("ycbcr_y.bmp");

	//对灰度图像加上高斯噪声
	CImg<int> image_gaussian_noise(width, height, 1, 1, 0);
	const double mean = 0.0;//均值
	const double stddev = 30;//标准差
	default_random_engine generator;
	normal_distribution<double> dist(mean, stddev);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			image_gaussian_noise(j, i) = image_ycbcr_y(j, i) + (int)dist(generator);
		}
	}
	image_gaussian_noise.display("高斯噪声");
	image_gaussian_noise.save("高斯噪声.bmp");

	//对灰度图像加上椒盐噪声
	CImg<int> image_salt_pepper_noise(width, height, 1, 1, 0);
	image_salt_pepper_noise = image_ycbcr_y;
	int n = 10000;
	for (int k = 0; k < n; k++)
	{
		int i = rand() % height;
		int j = rand() % width;
		image_salt_pepper_noise(j, i) = rand() % 2 ? 255 : 0;
	}
	image_salt_pepper_noise.display("椒盐噪声");
	image_salt_pepper_noise.save("椒盐噪声.bmp");

	//中值滤波 size x size
	CImg<int> image_salt_median_filted(width, height, 1, 1, 0);
	image_salt_median_filted = median_filter(image_salt_pepper_noise, 3);
	image_salt_median_filted.display("椒盐噪声中值滤波后");
	image_salt_median_filted.save("椒盐噪声中值滤波后.bmp");

	CImg<int> image_gaussian_median_filted(width, height, 1, 1, 0);
	image_gaussian_median_filted = median_filter(image_gaussian_noise, 3);
	image_gaussian_median_filted.display("高斯噪声中值滤波后");
	image_gaussian_median_filted.save("高斯噪声中值滤波后.bmp");
	/*(image_ycbcr_y - image_salt_median_filted).display();//残差，测试用*/

	//均值滤波 size x size
	CImg<int> image_salt_mean_filted(width, height, 1, 1, 0);
	image_salt_mean_filted = mean_filter(image_salt_pepper_noise, 5);
	image_salt_mean_filted.display("椒盐噪声均值滤波后");
	image_salt_mean_filted.save("椒盐噪声均值滤波后.bmp");

	CImg<int> image_gaussian_mean_filted(width, height, 1, 1, 0);
	image_gaussian_mean_filted = mean_filter(image_gaussian_noise, 5);
	image_gaussian_mean_filted.display("高斯噪声均值滤波后");
	image_gaussian_mean_filted.save("高斯噪声均值滤波后.bmp");
	/*(image_ycbcr_y - image_salt_median_filted).display();//残差，测试用*/

	//Sobel算子边缘检测
	int sobel[3][3] = { {-1, 0, 1},
						{-2, 0, 2},
						{-1, 0, 1} };
	CImg<int> sobel_filted(width, height, 1, 1, 0);
	sobel_filted = operator_filter(image_ycbcr_y, sobel);
	sobel_filted.display("Sobel边缘检测");
	sobel_filted.save("Sobel边缘检测.bmp");

	//离散哈尔小波变换滤波
	CImg<double> image_salt_wavelet_filted(width, height, 1, 1, 0);
	image_salt_wavelet_filted = wavelet_filter(image_salt_pepper_noise);

	CImg<double> image_gaussian_wavelet_filted(width, height, 1, 1, 0);
	image_gaussian_wavelet_filted = wavelet_filter(image_gaussian_noise);

	//PSNR SSIM
	double max_ssim = 1.02 * SSIM(image_ycbcr_y, image_ycbcr_y, 3);
	cout << "高斯噪声图像中值滤波PSNR: " << image_ycbcr_y.PSNR(image_gaussian_median_filted) << endl;
	cout << "高斯噪声图像中值滤波SSIM: " << SSIM(image_ycbcr_y, image_gaussian_median_filted, 3) / max_ssim << endl;
	cout << "椒盐噪声图像中值滤波PSNR: " << image_ycbcr_y.PSNR(image_salt_median_filted) << endl;
	cout << "椒盐噪声图像中值滤波SSIM: " << SSIM(image_ycbcr_y, image_salt_median_filted, 3) / max_ssim << endl;

	cout << "高斯噪声图像均值滤波PSNR: " << image_ycbcr_y.PSNR(image_gaussian_mean_filted) << endl;
	cout << "高斯噪声图像均值滤波SSIM: " << SSIM(image_ycbcr_y, image_gaussian_mean_filted, 3) / max_ssim << endl;
	cout << "椒盐噪声图像均值滤波PSNR: " << image_ycbcr_y.PSNR(image_salt_mean_filted) << endl;
	cout << "椒盐噪声图像均值滤波SSIM: " << SSIM(image_ycbcr_y, image_salt_mean_filted, 3) / max_ssim << endl;

	cout << "高斯噪声图像小波滤波PSNR: " << image_ycbcr_y.PSNR(image_gaussian_wavelet_filted) << endl;
	cout << "高斯噪声图像小波滤波SSIM: " << SSIM(image_ycbcr_y, image_gaussian_wavelet_filted, 3) / max_ssim << endl;
	cout << "椒盐噪声图像小波滤波PSNR: " << image_ycbcr_y.PSNR(image_salt_wavelet_filted) << endl;
	cout << "椒盐噪声图像小波滤波SSIM: " << SSIM(image_ycbcr_y, image_salt_wavelet_filted, 3) / max_ssim << endl;

}