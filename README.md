## Thuật Toán Đom Đóm (Firefly Algorithm)

Cài đặt thuật toán Đom Đóm bằng Python cho các bài toán tối ưu hóa, bao gồm triển khai Mạng Cảm Biến Không Dây (WSN) và các bài toán tối ưu hóa tổng quát.

### Tính Năng

- **Nhiều Biến Thể Thuật Toán Đom Đóm**:
  - Thuật toán Đom Đóm Chuẩn (FFA - Standard Firefly Algorithm)
  - Thuật toán Đom Đóm Cải Tiến (IFA - Improved Firefly Algorithm)
  - Thuật toán Đom Đóm Lévy Flight (LFA - Lévy Flight Firefly Algorithm)
  - Thuật toán Đom Đóm Đa Mục Tiêu (MOFA - Multi-Objective Firefly Algorithm)

- **Ứng Dụng**:
  - Tối ưu hóa triển khai Mạng Cảm Biến Không Dây (WSN)
  - Tối ưu hóa hàm số tổng quát
  - Bài toán tối ưu hóa có ràng buộc

### Cấu Trúc Dự Án

#### Cài Đặt Chính cho WSN (`Firefly-Algorithm-WSN/`)

- `FA.py`: Điểm khởi chạy chính cho triển khai WSN
- `init_ffa.py`: Khởi tạo quần thể đom đóm
- `ffa_wsn.py`, `ifa_wsn.py`, `lfa_wsn.py`, `mofa_wsn.py`: Các cài đặt triển khai WSN
- `ffa_move.py`, `ifa_move.py`, `lfa_move.py`, `mofa_move.py`: Các chiến lược cập nhật di chuyển
- `coverage.py`: Tính toán độ phủ của WSN
- `findlimits.py`: Xử lý ràng buộc biên
- `draw.py`: Công cụ trực quan hóa
- `run_ifa.py`, `run_lfa.py`, `run_mofa.py`: Chạy các biến thể thuật toán cụ thể

#### Mã Nguồn (`Firefly-Algorithm-WSN/src/`)

- `firefly_simple.py`: Cài đặt đơn giản thuật toán đom đóm
- `fa_mincon.py`: Thuật toán đom đóm cho tối ưu hóa có ràng buộc
- `alpha_new.py`: Các chiến lược tham số thích ứng

#### Demo (`demo/`)

- `baitoandemo.py`: Ví dụ minh họa với chú thích chi tiết bằng tiếng Việt

### Yêu Cầu

```python
numpy
matplotlib
```

### Cài Đặt

```bash
pip install numpy matplotlib
```

### Cách Sử Dụng

#### Triển Khai WSN Cơ Bản

```bash
cd Firefly-Algorithm-WSN
python FA.py
```

#### Chạy Các Biến Thể Thuật Toán

```bash
python run_ifa.py  # Thuật toán Đom Đóm Cải Tiến
python run_lfa.py  # Thuật toán Đom Đóm Lévy Flight
python run_mofa.py # Thuật toán Đom Đóm Đa Mục Tiêu
```

#### Ví Dụ Demo

```bash
cd demo
python baitoandemo.py
```

### Tham Số

- `w`: Chiều rộng vùng triển khai (mặc định: 100)
- `d`: Số chiều/số cảm biến (mặc định: 100)
- `r`: Bán kính hoạt động của cảm biến (mặc định: 7)
- `para`: Các tham số thuật toán [n, MaxGeneration, alpha, betamin, gamma]
  - `n`: Kích thước quần thể (mặc định: 25)
  - `MaxGeneration`: Số vòng lặp tối đa (mặc định: 5)
  - `alpha`: Hệ số ngẫu nhiên (mặc định: 0.7)
  - `betamin`: Độ hấp dẫn tối thiểu (mặc định: 0.2)
  - `gamma`: Hệ số hấp thụ ánh sáng (mặc định: 1)

### Giới Thiệu Thuật Toán

Thuật toán Đom Đóm là một thuật toán tối ưu hóa siêu heuristic lấy cảm hứng từ thiên nhiên, dựa trên hành vi phát sáng của đom đóm. Các khái niệm chính:

- **Độ Hấp Dẫn**: Giảm theo khoảng cách và phụ thuộc vào cường độ ánh sáng
- **Di Chuyển**: Các đom đóm di chuyển về phía các giải pháp sáng hơn (tốt hơn)
- **Ngẫu Nhiên Hóa**: Thành phần di chuyển ngẫu nhiên để khám phá không gian tìm kiếm


### Tác Giả

NguyenManhTuan-20235862
NguyenTuanAnh-20235650
