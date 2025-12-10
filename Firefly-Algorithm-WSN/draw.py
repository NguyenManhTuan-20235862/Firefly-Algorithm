import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Thêm fval vào chữ ký hàm để hiển thị kết quả
def draw(sx, sy, SensorNum, R, fval=None):
    alpha1 = np.linspace(0, 2*np.pi, 100)
    
    # KHẮC PHỤC LỖI INDEXERROR: Gán trực tiếp vì sx, sy là các vector 1 chiều
    xp = sx 
    yp = sy
    
    # Thiết lập kích thước đồ thị
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # VẼ VÙNG PHỦ SÓNG (Hình tròn cho từng nút cảm biến)
    for Pxy in range(SensorNum):
        center_x = sx[Pxy]
        center_y = sy[Pxy]
        
        # Tạo tọa độ hình tròn
        x = center_x + R * np.cos(alpha1)
        y = center_y + R * np.sin(alpha1)
        
        # Tô màu vùng phủ sóng (màu xanh dương, độ trong suốt 30%)
        ax.fill(x, y, 'b', alpha=0.3)
    
    # VẼ NÚT CẢM BIẾN (Dấu chấm)
    ax.plot(xp, yp, 'k.', markersize=5) 
    
    # Đặt giới hạn và tỷ lệ cho đồ thị
    ax.set_aspect('equal')
    ax.set_xlim([0, 100])
    ax.set_ylim([0, 100])
    ax.grid(True)
    
    # TẠO TIÊU ĐỀ (Bao gồm kết quả vùng phủ sóng)
    title_text = f'WSN Coverage Map (R={R}, N={SensorNum})'
    if fval is not None:
        title_text += f'\nFinal Coverage: {fval:.4f}' # Hiển thị 4 chữ số thập phân
        
    plt.title(title_text) 
    plt.xlabel('X Coordinate (Area)') 
    plt.ylabel('Y Coordinate (Area)') 
    
    # Không cần plt.show() ở đây, nó sẽ được gọi từ file chính