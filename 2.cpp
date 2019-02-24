#include <iostream>
#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <string>
//#include <cassert>
#include "CImg.h"
using namespace std;
using namespace cimg_library;
complex<double> complex_i(0, 1);
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
CImg<complex<double>> DFT(CImg<complex<double>> image)
{
	int width = image.width();
	int height = image.height();
	CImg<complex<double>> image_dft(width, height, 1, 1); // 这里加了第5个参数（初始化参数）会有问题
	// 每个8x8块以它左上的像素下标表示
	for (int block_i = 0; block_i < height; block_i +=8)
	{
		for (int block_j = 0; block_j < width; block_j +=8)
		{
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++)
				{
					for (int i1 = 0; i1 < 8; i1++)
					{
						for (int j1 = 0; j1 < 8; j1++)
						{
							image_dft(block_j + j, block_i + i) += exp(-2 * M_PI * complex_i * (double)(i * i1 + j * j1) / 8.0) * image(block_j + j1, block_i + i1);
						}
					}
				}
			}
		}
	}
	return image_dft;
}
 CImg<complex<double>> IDFT(CImg<complex<double>> image)
{
	int width = image.width();
	int height = image.height();
	CImg<complex<double>> image_idft(width, height, 1, 1);
	
	for (int block_i = 0; block_i < height; block_i += 8)
	{
		for (int block_j = 0; block_j < width; block_j += 8)
		{
			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++)
				{
					for (int i1 = 0; i1 < 8; i1++)
					{
						for (int j1 = 0; j1 < 8; j1++)
						{
							image_idft(block_j + j, block_i + i) += exp(2 * M_PI * complex_i * (double)(i * i1 + j * j1) / 8.0) * image(block_j + j1, block_i + i1) / 64.0;
						}
					}
				}
			}
		}
	}
	return image_idft;
}
 CImg<double> amplitude(CImg<complex<double>> image)
 {
	 int width = image.width();
	 int height = image.height();
	 CImg<double> image_amplitude(width, height, 1, 1, 0);
	 for (int i = 0; i < height; i++)
	 {
		 for (int j = 0; j < width; j++)
		 {
			 image_amplitude(j, i) = abs(image(j, i));
		 }
	 }
	 return image_amplitude;
 }
 CImg<double> phase(CImg<complex<double>> image)
 {
	 int width = image.width();
	 int height = image.height();
	 CImg<double> image_phase(width, height, 1, 1, 0);
	 for (int i = 0; i < height; i++)
	 {
		 for (int j = 0; j < width; j++)
		 {
			 image_phase(j, i) = arg(image(j, i));
		 }
	 }
	 return image_phase;
 }
 CImg<complex<double>> exp_i_phase(CImg<complex<double>> image)
 {
	 int width = image.width();
	 int height = image.height();
	 CImg<complex<double>> image_phase(width, height, 1, 1);
	 for (int i = 0; i < height; i++)
	 {
		 for (int j = 0; j < width; j++)
		 {
			 image_phase(j, i) = exp(complex_i * arg(image(j, i)));
		 }
	 }
	 return image_phase;
 }

 CImg<double> DCT(CImg<double> image, int reserve_matrix[8][8])
 {
	 int width = image.width();
	 int height = image.height();
	 CImg<double> image_dct(width, height, 1, 1, 0); // 这里加了第5个参数（初始化参数）会有问题
	 // 每个8x8块以它左上的像素下标表示
	 for (int block_i = 0; block_i < height; block_i += 8)
	 {
		 for (int block_j = 0; block_j < width; block_j += 8)
		 {
			 for (int i = 0; i < 8; i++)
			 {
				 for (int j = 0; j < 8; j++)
				 {
					 for (int i1 = 0; i1 < 8; i1++)
					 {
						 for (int j1 = 0; j1 < 8; j1++)
						 {
							 image_dct(block_j + j, block_i + i) += cos(M_PI * i * (i1 + 0.5) / 8.0) * cos(M_PI * j * (j1 + 0.5) / 8.0) * image(block_j + j1, block_i + i1);
						 }
					 }
					 image_dct(block_j + j, block_i + i) *= reserve_matrix[i][j];
				 }
			 }
		 }
	 }
	 return image_dct;
 }
 CImg<double> IDCT(CImg<double> image)
 {
	 int width = image.width();
	 int height = image.height();
	 CImg<double> image_idct(width, height, 1, 1, 0); // 这里加了第5个参数（初始化参数）会有问题
	 // 每个8x8块以它左上的像素下标表示
	 for (int block_i = 0; block_i < height; block_i += 8)
	 {
		 for (int block_j = 0; block_j < width; block_j += 8)
		 {
			 for (int i = 0; i < 8; i++)
			 {
				 for (int j = 0; j < 8; j++)
				 {
					 image_idct(block_j + j, block_i + i) += image(block_j, block_i) * 2.0 / 128.0;
					 for (int i1 = 0; i1 < 8; i1++)
					 {
						 for (int j1 = 0; j1 < 8; j1++)
						 {
							 if (!(i1 == 0 && j1 == 0))
							 {
								 image_idct(block_j + j, block_i + i) += cos(M_PI * i1 * (i + 0.5) / 8.0) * cos(M_PI * j1 * (j + 0.5) / 8.0) * image(block_j + j1, block_i + i1) * 2.0 / 64.0;
							 }
						 }
					 }
				 }
			 }
		 }
	 }
	 return image_idct;
 }
void DCT_display_PSNR(CImg<double> image_ycbcr_y, int reserve_matrx[8][8], int number)
{
	//DCT
	CImg<double> image_8x8DCT = DCT(image_ycbcr_y, reserve_matrx);
	//IDCT
	CImg<double> image_8x8IDCT = IDCT(image_8x8DCT);
	//显示还原后的图像
	string p = "保留" + to_string(number) + "个系数.bmp";
	image_8x8IDCT.display(p.data());
	image_8x8IDCT.save(p.data());
	//计算PSNR
	cout << "\nPSNR(" << "保留" << number << "个系数):" << image_ycbcr_y.PSNR(image_8x8IDCT) << endl << endl;
}
int main(int argc, char * argv[])
{
	// 读取文件，获得原图像
	CImg<unsigned char> image = readBMP(argv[1]);
	int height = image.height();
	int width = image.width();
	double rgb2ycbcr[3][3] = { {0.299, 0.587, 0.114}, {0.500, -0.4187, -0.0813}, {-0.1687, -0.3313, 0.5} };

	CImg<double> image_ycbcr_y(width, height, 1, 1, 0);
	CImg<double> image_ycbcr_cb(width, height, 1, 1, 0);
	CImg<double> image_ycbcr_cr(width, height, 1, 1, 0);

	// 转换到ycbcr颜色空间
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
	// 显示原图和ycbcr图
	image.display("bmp");
	image.save("raw.bmp");
	image_ycbcr_y.display("ycbcr_y");
	image_ycbcr_y.save("ycbcr_y.bmp");



	///////////////////////////////////////DFT//////////////////////////////
	// 计时开始
	clock_t start = clock();
	// DFT
	CImg<complex<double>> image_8x8DFT = DFT(image_ycbcr_y);
	// 计时结束
	clock_t end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "\n傅里叶变换时间：" << time << "s" << endl << endl;

	// 频域取幅度
	CImg<double> image_8x8DFT_amplitude = amplitude(image_8x8DFT);
	// 频域取相位
	CImg<double> image_8x8DFT_phase = phase(image_8x8DFT);
	// 显示幅度和相位
	image_8x8DFT_amplitude.display("傅里叶变换后的频域幅度");
	image_8x8DFT_phase.display("傅里叶变换后的频域相位");
	image_8x8DFT_amplitude.save("傅里叶变换后的频域幅度.bmp");
	image_8x8DFT_phase.save("傅里叶变换后的频域相位.bmp");

	// 幅度IDFT
	CImg<double> image_8x8IDFT_amplitude = amplitude(IDFT(image_8x8DFT_amplitude));
	// 频域取相位的exp()，IDFT
	CImg<complex<double>> image_8x8DFT_exp_i_phase = exp_i_phase(image_8x8DFT);
	CImg<double> image_8x8IDFT_phase = amplitude(IDFT(image_8x8DFT_exp_i_phase));
	// 显示幅度和相位的IDFT
	image_8x8IDFT_amplitude.display("由幅度还原的图像");
	image_8x8IDFT_phase.display("由相位还原的图像");
	image_8x8IDFT_amplitude.save("由幅度还原的图像.bmp");
	image_8x8IDFT_phase.normalize(0, 255).save("由相位还原的图像.bmp");
	// 显示原图像傅里叶变换后的图像，分块做逆变换的图像，测试算法正确性
	/*CImg<double> image_IDFT_test = amplitude(IDFT(image_8x8DFT));
	image_IDFT_test.display();*/




	/////////////////////////////////DCT//////////////////////////////
	// 保留频域像素位置的矩阵
	int reserve_matrx1[8][8] = { {1,0,0,0,0,0,0,0},
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} };
	int reserve_matrx2[8][8] = { {1,1,0,0,0,0,0,0},
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} };
	int reserve_matrx4[8][8] = { {1,1,0,0,0,0,0,0},
								{1,1,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} };
	int reserve_matrx6[8][8] = { {1,1,1,0,0,0,0,0},
								{1,1,0,0,0,0,0,0} ,
								{1,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} };
	int reserve_matrx8[8][8] = { {1,1,1,1,0,0,0,0},
								{1,1,0,0,0,0,0,0} ,
								{1,0,0,0,0,0,0,0} ,
								{1,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} };
	int reserve_matrx10[8][8] ={{1,1,1,1,0,0,0,0},
								{1,1,1,0,0,0,0,0} ,
								{1,1,0,0,0,0,0,0} ,
								{1,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} ,
								{0,0,0,0,0,0,0,0} };
	DCT_display_PSNR(image_ycbcr_y, reserve_matrx1, 1);
	DCT_display_PSNR(image_ycbcr_y, reserve_matrx2, 2);
	DCT_display_PSNR(image_ycbcr_y, reserve_matrx4, 4);
	DCT_display_PSNR(image_ycbcr_y, reserve_matrx6, 6);
	DCT_display_PSNR(image_ycbcr_y, reserve_matrx8, 8);
	DCT_display_PSNR(image_ycbcr_y, reserve_matrx10, 10);

	//显示8x8DCT基图像
	CImg<double> DCT_base_function(8*8, 8*8, 1, 1, 0);
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			for (int i1 = 0; i1 < 8; i1++)
			{
				for (int j1 = 0; j1 < 8; j1++)
				{
					DCT_base_function(i * 8 + i1, j * 8 + j1) = cos(M_PI * i * (i1 + 0.5) / 8.0) * cos(M_PI * j * (j1 + 0.5) / 8.0);
				}
			}
		}
	}
	DCT_base_function.display();
	DCT_base_function.normalize(0, 255).save("base_function.bmp");


}