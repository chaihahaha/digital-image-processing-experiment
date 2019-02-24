import cv2
import numpy as np
import os
import matplotlib.pyplot as plt


def error_diffusion(input_arr):
    a = input_arr.copy()
    height, width = a.shape
    for i in range(height):
        for j in range(width):
            error = a[i, j] - 127
            a[i, j] = 255 if error > 0 else 0
            if i + 1 < height and j - 1 >= 0:
                a[i + 1, j - 1] += error * 3/16
            if i + 1 < height:
                a[i + 1, j] += error * 5/16
            if j + 1 < width:
                a[i, j + 1] += error * 7/16
            if i + 1 < height and j + 1 < width:
                a[i + 1, j + 1] += error * 1/16
    return a


def gaussian_noise(a):
    return (a + np.random.normal(0, 30, a.shape)).astype(np.int)


def salt_pepper_noise(input_arr, k):
    a = input_arr.copy()
    rand_arr = np.random.uniform(0, 1, a.shape)
    a[rand_arr < k] = 0
    a[rand_arr < 2*k] = 255
    return a


def compute_thickness(file_name):
    # 输入：细胞图像文件名
    # 细胞图像厚度计算使用以圆心为中心的12个方向的线扫描，对可能的膜（边缘）位置做聚类，得到5个膜的位置，再计算厚度
    # 输出：细胞各层厚度
    cell = cv2.imread(file_name, 0)  # 读入图像
    print('\n正在处理的文件名:', file_name)

    mean_center = get_center(cell)

    # cell_tmp = cell.copy()
    # for i in circles[0]:
    #     # 画圆
    #     cv2.circle(cell_tmp, (i[0], i[1]), i[2], (255, 255, 255), 2)
    #     # 画圆心
    #     cv2.circle(cell_tmp, (i[0], i[1]), 2, (255, 255, 255), 3)
    # cell1 = cv2.resize(src=cell_tmp, dsize=(800, 600), interpolation=cv2.INTER_AREA)
    # plt.imshow(cell1, cmap='gray')
    # plt.show()

    cell_cannyed = cv2.Canny(cell, 50, 80)
    # cell1 = cv2.resize(src=cell_cannyed, dsize=(800, 600), interpolation=cv2.INTER_AREA)
    # plt.imshow(cell1, cmap='gray')
    # plt.show()

    # scan_and_plot(cell_cannyed, mean_center[0], mean_center[1], 900)
    radius_points = scaned_radius(cell_cannyed, mean_center[0], mean_center[1], 900)
    radius_points = [i for i in radius_points if i > 350]
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0)
    ret, label, clustered_radius_points = cv2.kmeans(np.array(radius_points).astype(np.float32), 5, None, criteria, 10,
                                            cv2.KMEANS_RANDOM_CENTERS)
    sorted_radius_points = np.sort(clustered_radius_points.reshape(-1))
    print('聚类结果，膜与中心距离：', sorted_radius_points)
    thickness = []
    for i in range(4):
        thickness.append(sorted_radius_points[i + 1] - sorted_radius_points[i])
    return thickness


def get_center(cell):
    # 输入：细胞图像
    # 使用霍夫圆变换，通过圆心位置计算中心位置
    # 输出：圆心坐标（x,y）
    circles = None
    p1, p2 = 200, 200  # 80 90正常
    while type(circles) != np.ndarray:
        p1 -= 40
        p2 -= 40
        circles = cv2.HoughCircles(cell, cv2.HOUGH_GRADIENT, 1, 0.001, param1=p1, param2=p2, minRadius=500, maxRadius=1000)
    circle_center = circles[0, :, :2]
    mean_center = np.mean(circle_center, 0).tolist()
    mean_center = [int(i) for i in mean_center]
    print('圆心位置: ', mean_center)
    return mean_center


# def scan_and_plot(cell_cannyed, x, y, radius):
#     # 以一定半径进行线扫描，并绘制此线上的图像值
#     for i in range(1, 13):
#         scan_result = scan_in_a_direction(cell_cannyed, x, y, radius, 2 * np.pi / i, radius * 2)
#         plt.plot(np.arange(-radius, radius), scan_result)
#         plt.show()


def scaned_radius(cell_cannyed, x, y, radius):
    # 输入：边缘检测后的细胞图像cell_cannyed，圆心坐标(x,y)，半径radius（以像素为单位）
    # 输出：扫描线上的膜可能所在位置（像素值为255）的点，即膜与圆心距离（每个像素值为255的像素的索引构成的数组）
    radius_points = []
    for i in range(1, 13):
        scan_result = scan_in_a_direction(cell_cannyed, x, y, radius, 2 * np.pi / i, radius * 2)
        radius_points.extend(combine_halves(scan_result, radius))
    return radius_points


def scan_in_a_direction(matrix, x, y, radius, angle, resolution):
    # 输入: 图像matrix，中心点坐标(x,y)，扫描半径radius，扫描角度angle，分辨率resolution（输出的图像像素值数组长度）
    # 输出：图像任一中心点（x,y）以一倾斜角度向两侧扫描的像素值数组（直径）
    height, width = matrix.shape
    y0, x0 = max(0, x - radius * np.cos(angle)), max(0, y - radius * np.sin(angle))
    y1, x1 = min(width - 1, x + radius * np.cos(angle)), min(height - 1, y + radius * np.sin(angle))
    x_index, y_index = np.linspace(x0, x1, resolution).astype(np.int), np.linspace(y0, y1, resolution).astype(np.int)
    return [matrix[x_index[i], y_index[i]] for i in range(resolution)]


def combine_halves(array, radius):
    # 输入：图像任一中心点（x,y）以一倾斜角度向两侧扫描的像素值数组array（直径），半径radius（以像素为单位）
    # 输出：将直径数组array对折，输出长度为array一半的数组（半径）
    l = len(array)
    half = int(l/2)
    array1, array2 = np.array(array[: half]), np.array(array[half: l])
    return (np.hstack([half - np.where(array1 != 0)[0], np.where(array2 != 0)[0]]) * 2 * radius / l).tolist()


img = cv2.imread('lena.bmp', 0)   # 读入图像

th, dst = cv2.threshold(img, 0, 255, cv2.THRESH_OTSU)   # OTSU二值化
print("OTSU阈值：", th)
plt.imshow(dst, cmap='gray')
plt.show()

kernel = np.ones((5, 5), np.uint8)
erosion = cv2.erode(dst, kernel, iterations=1)      # 腐蚀
plt.imshow(erosion, cmap='gray')
plt.show()

binaryzation = error_diffusion(img)     # 基于误差扩散的半色调二值化
plt.imshow(binaryzation, cmap='gray')
plt.show()

gaussian = gaussian_noise(img)          # 产生高斯噪声
salt = salt_pepper_noise(img, 0.05)     # 产生椒盐噪声
plt.imshow(gaussian, cmap='gray')
plt.show()
plt.imshow(salt, cmap='gray')
plt.show()

median = cv2.medianBlur(salt, 5)
plt.imshow(median, cmap='gray')
plt.show()

#########################
# 扫描picture文件下的所有细胞图片，计算每层平均厚度
layer_thickness = []
for root, dirs, files in os.walk('pictures'):
    for file_name in files:
        layer_thickness.append(compute_thickness(root + '\\' + file_name))
print('平均厚度: ', np.mean(layer_thickness, 0))   # [  76.07469177   79.96934509  101.90732574  122.3502655 ]
