from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy  # type: ignore
import matplotlib.pyplot as plt

# Locate all .vti files
vti_files = sorted(Path("out/").glob("*.vti"))
reader = vtk.vtkXMLImageDataReader()
if not vti_files:
    raise FileNotFoundError("No .vti files found in the current directory.")

# Read first file to get dimensions
reader.SetFileName(str(vti_files[0]))
reader.Update()
image_data = reader.GetOutput()
dims = image_data.GetDimensions()  # (nx, ny, nz)

print(f"Found {len(vti_files)} files")
print(f"Grid dimensions: {dims}")

# Initialize list to store all data
all_data = []

# Read all files
for vti_file in vti_files:
    reader.SetFileName(str(vti_file))
    reader.Update()
    image_data = reader.GetOutput()

    # Extract velocity data
    point_data = image_data.GetPointData()
    if point_data.GetNumberOfArrays() == 0:
        raise ValueError("No point data arrays found in the VTI file.")

    # Array 0 is pressure, Array 1 is velocity vector (u, v, 0)
    velocityArray = point_data.GetArray(1)  # This is the velocity vector with 3 components

    # Convert to NumPy - this gives us shape (n_points, 3) where columns are [u, v, w]
    velocity_np = vtk_to_numpy(velocityArray)

    # Extract u (component 0) and v (component 1)
    u_np = velocity_np[:, 0]  # u component
    v_np = velocity_np[:, 1]  # v component

    # Stack u and v as channels: shape (nx*ny*nz, 2)
    # Then reshape to (ny, nx, 2) for a 2D field
    uv_combined = np.stack([u_np, v_np], axis=-1)  # (nx*ny*nz, 2)
    uv_reshaped = uv_combined.reshape((dims[1], dims[0], 2))  # (ny, nx, 2)

    all_data.append(uv_reshaped)

# Concatenate all data: shape (num_files, ny, nx, 2)
data_matrix = np.array(all_data)
print(f"Combined data shape: {data_matrix.shape}")
print("Shape format: (data_points, field_width, field_height, channels)")

# Save combined data
np.save("trainingData.npy", data_matrix)
print("Preprocessed data saved to trainingData.npy")

# Select a random simulation
random_idx = np.random.randint(0, data_matrix.shape[0])
random_simulation = data_matrix[random_idx]  # Shape: (ny, nx, 2)

u_field = random_simulation[:, :, 0]
v_field = random_simulation[:, :, 1]

# Plot u and v components
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

im1 = axes[0].imshow(u_field, origin="lower", cmap="RdBu_r")
axes[0].set_title(f"U component (simulation {random_idx})")
axes[0].set_xlabel("X index")
axes[0].set_ylabel("Y index")
plt.colorbar(im1, ax=axes[0], label="U velocity")

im2 = axes[1].imshow(v_field, origin="lower", cmap="RdBu_r")
axes[1].set_title(f"V component (simulation {random_idx})")
axes[1].set_xlabel("X index")
axes[1].set_ylabel("Y index")
plt.colorbar(im2, ax=axes[1], label="V velocity")

plt.tight_layout()
# plt.savefig('random_simulation_uv.png', dpi=150, bbox_inches='tight')
# print("Plot saved to random_simulation_uv.png")
plt.show()
plt.close()
