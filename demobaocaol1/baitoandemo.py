import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# CẤU HÌNH THUẬT TOÁN
# ==========================================
class FireflyAlgorithm:
    def __init__(self, func, dim, pop_size=20, max_gen=50, alpha=0.5, beta_0=1.0, gamma=1.0):
        self.func = func          # Hàm mục tiêu cần tối ưu
        self.dim = dim            # Số chiều của bài toán 
        self.pop_size = pop_size  # Số lượng đom đóm
        self.max_gen = max_gen    # Số vòng lặp tối đa
        self.alpha = alpha        # Hệ số ngẫu nhiên (Randomness)
        self.beta_0 = beta_0      # Độ hấp dẫn tại khoảng cách r=0
        self.gamma = gamma        # Hệ số hấp thụ ánh sáng (Light Absorption)
        
        # Khởi tạo quần thể đom đóm ngẫu nhiên trong khoảng [-5, 5]
        self.fireflies = np.random.uniform(-5, 5, (self.pop_size, self.dim))
        self.intensity = np.zeros(self.pop_size)
        
    def calculate_intensity(self):
        """Tính độ sáng (nghịch đảo của giá trị hàm mục tiêu vì ta đang tìm Min)"""
        for i in range(self.pop_size):
            # Với bài toán tìm Min, giá trị hàm càng nhỏ -> độ sáng càng lớn
            # Ta dùng trực tiếp giá trị hàm để so sánh: bé hơn là "sáng hơn" (tốt hơn)
            self.intensity[i] = self.func(self.fireflies[i])

    def run(self):
        # Lưu lại lịch sử để vẽ đồ thị
        history = []
        
        # Đánh giá ban đầu
        self.calculate_intensity()
        best_val = np.min(self.intensity)
        print(f"Khởi tạo: Best Value = {best_val:.5f}")

        # Vòng lặp chính (Generations)
        for t in range(self.max_gen):
            
            # Loop lồng nhau để so sánh từng cặp đom đóm
            for i in range(self.pop_size):
                for j in range(self.pop_size):
                    
                    # Nếu đom đóm j sáng hơn đom đóm i (j tốt hơn i)
                    # Lưu ý: Bài toán Min, giá trị nhỏ hơn nghĩa là sáng hơn
                    if self.intensity[j] < self.intensity[i]:
                        
                        # 1. Tính khoảng cách Euclid giữa i và j
                        r = np.linalg.norm(self.fireflies[i] - self.fireflies[j])
                        
                        # 2. Tính độ hấp dẫn (Beta)
                        # Beta = Beta0 * exp(-gamma * r^2)
                        beta = self.beta_0 * np.exp(-self.gamma * (r ** 2))
                        
                        # 3. Cập nhật vị trí của đom đóm i (di chuyển về phía j)
                        # Công thức: x_new = x_old + Attraction + Randomness
                        tmp_pos = self.fireflies[i] + \
                                  beta * (self.fireflies[j] - self.fireflies[i]) + \
                                  self.alpha * (np.random.rand(self.dim) - 0.5)
                        
                        # Đánh giá vị trí mới
                        tmp_intensity = self.func(tmp_pos)
                        
                        # Nếu vị trí mới tốt hơn, cập nhật
                        if tmp_intensity < self.intensity[i]:
                            self.fireflies[i] = tmp_pos
                            self.intensity[i] = tmp_intensity

            # Giảm hệ số ngẫu nhiên theo thời gian (để hội tụ dần)
            self.alpha *= 0.97
            
            # Ghi log kết quả tốt nhất hiện tại
            current_best = np.min(self.intensity)
            history.append(current_best)
            if t % 10 == 0:
                print(f"Gen {t}: Best Value = {current_best:.5f}")

        # Tìm kết quả cuối cùng
        best_index = np.argmin(self.intensity)
        return self.fireflies[best_index], self.intensity[best_index], history

# ==========================================
# ĐỊNH NGHĨA BÀI TOÁN & CHẠY
# ==========================================

# Hàm mục tiêu: Sphere Function f(x) = x1^2 + x2^2 + ...
def sphere_function(x):
    return np.sum(x**2)

if __name__ == "__main__":
    # Cấu hình: 2 chiều (để dễ hình dung x, y), 20 con đom đóm, chạy 50 vòng
    fa = FireflyAlgorithm(func=sphere_function, dim=2, pop_size=30, max_gen=50, gamma=1.0)
    
    best_pos, best_score, history_data = fa.run()
    
    print("-" * 30)
    print(f"KẾT QUẢ CUỐI CÙNG:")
    print(f"Vị trí tối ưu (x, y): {best_pos}")
    print(f"Giá trị hàm mục tiêu: {best_score:.10f}") # Càng gần 0 càng tốt
    print("-" * 30)

    # ==========================================
    # VẼ BIỂU ĐỒ 
    # ==========================================
    # Biểu đồ 1: Sự hội tụ của thuật toán
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(history_data)
    plt.title('Quá trình hội tụ (Convergence)')
    plt.xlabel('Thế hệ (Generation)')
    plt.ylabel('Giá trị tốt nhất (Best Value)')
    plt.grid(True)

    # Biểu đồ 2: Vị trí cuối cùng của bầy đom đóm
    plt.subplot(1, 2, 2)
    # Vẽ các điểm đom đóm cuối cùng
    plt.scatter(fa.fireflies[:, 0], fa.fireflies[:, 1], c='red', label='Đom đóm')
    # Vẽ điểm đích (0,0)
    plt.scatter(0, 0, c='green', marker='x', s=100, label='Đích (Global Min)')
    plt.title('Vị trí cuối cùng trong không gian 2D')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()