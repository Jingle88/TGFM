import numpy as np

# Load the .npy file
npy_matrix = np.load("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/ld_npy_perwindow/22:15416047-18416047_ld.npy")

# Read the .vcor2 file
with open("/nfs/users/nfs_j/jh59/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/LD_matrix_file_perwindow/22:15416047-18416047.unphased.vcor2", "r") as f:
    vcor_lines = f.readlines()

# Remove comments and empty lines
vcor_data = [line.strip().split() for line in vcor_lines if line.strip() and not line.startswith("#")]
vcor_matrix = np.array(vcor_data, dtype=float)

# Check dimensions
assert npy_matrix.shape == vcor_matrix.shape, "ERROR: Matrix dimensions do not match!"
print("✅ Step 1: Matrix dimensions match.")

# Check element-wise equality with a tolerance for floating-point errors
if np.allclose(npy_matrix, vcor_matrix, atol=1e-6):
    print("✅ Step 2: Matrix values match within tolerance.")
else:
    print("❌ ERROR: Matrix values differ!")

if np.allclose(npy_matrix, npy_matrix.T, atol=1e-6):
    print("✅ Step 3: Matrix is symmetric.")
else:
    print("⚠️ WARNING: Matrix is not perfectly symmetric.")

if np.allclose(np.diag(npy_matrix), np.ones(npy_matrix.shape[0]), atol=1e-6):
    print("✅ Step 4: Diagonal values are correct (≈1s).")
else:
    print("⚠️ WARNING: Diagonal values are incorrect.")

import random

indices = random.sample(range(npy_matrix.shape[0]), min(5, npy_matrix.shape[0]))  # Pick 5 random indices

for i in indices:
    for j in indices:
        if not np.isclose(npy_matrix[i, j], vcor_matrix[i, j], atol=1e-6):
            print(f"❌ ERROR: Entry mismatch at ({i}, {j}) -> npy: {npy_matrix[i, j]}, vcor2: {vcor_matrix[i, j]}")
            break
else:
    print("✅ Step 5: Random subset of values match.")

print("\n📊 Summary statistics comparison:")
print(f"Mean:  npy={npy_matrix.mean():.6f}, vcor2={vcor_matrix.mean():.6f}")
print(f"Std:   npy={npy_matrix.std():.6f}, vcor2={vcor_matrix.std():.6f}")
print(f"Min:   npy={npy_matrix.min():.6f}, vcor2={vcor_matrix.min():.6f}")
print(f"Max:   npy={npy_matrix.max():.6f}, vcor2={vcor_matrix.max():.6f}")

if np.allclose(npy_matrix.mean(), vcor_matrix.mean(), atol=1e-6) and \
   np.allclose(npy_matrix.std(), vcor_matrix.std(), atol=1e-6):
    print("✅ Step 6: Global statistics match.")
else:
    print("⚠️ WARNING: Statistical differences detected!")
