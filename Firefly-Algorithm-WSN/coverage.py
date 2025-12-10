import numpy as np

def coverage(Solution, WSNpoint, R):
    CoverNum = 0
    
    # [Tối ưu hóa] Giảm kích thước lưới kiểm tra 
    GRID_SIZE = 50 
    
    for Px in range(GRID_SIZE):
        for Py in range(GRID_SIZE):
            # [Tối ưu hóa] Biến đổi tọa độ để điểm kiểm tra trải rộng toàn bộ khu vực 100x100
            px_real = Px * (100 / GRID_SIZE)
            py_real = Py * (100 / GRID_SIZE)
            
            for S in range(WSNpoint):
                # Tính khoảng cách bình phương
                Distance = (Solution[0, S] - px_real)**2 + (Solution[1, S] - py_real)**2
                
                if Distance < (R * R):
                    CoverNum += 1
                    break
    
    # Tính tỷ lệ trên lưới kiểm tra mới
    z = CoverNum / (GRID_SIZE * GRID_SIZE) 
    return z