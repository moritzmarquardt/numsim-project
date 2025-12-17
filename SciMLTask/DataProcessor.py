import pyvista as pv
import numpy as np
import os
from sklearn.model_selection import train_test_split
import yaml
import torch


class DataProcessor:

    def __init__(self, dataset_name):
        self.dataset_name = dataset_name

    def extractFromFiles(self, folder_path):
        vti_files = sorted(
            [
                os.path.join(folder_path, f)
                for f in os.listdir(folder_path)
                if f.endswith(".vti")
            ]
        )
        if not vti_files:
            raise FileNotFoundError("No .vti files found in the current directory.")

        print(f"Found {len(vti_files)} files")

        # Read first file to get dimensions
        mesh = pv.read(str(vti_files[0]))
        dims = mesh.dimensions  # (nx, ny, nz)
        print(f"Grid dimensions: {dims}")

        # Initialize list to store all data
        all_data = []

        # Read all files
        for vti_file in vti_files:
            mesh = pv.read(str(vti_file))

            velocityArray = mesh.point_data[
                "velocity"
            ]  # This is the velocity vector with 3 components

            u_np = velocityArray[:, 0]  # u component
            v_np = velocityArray[:, 1]  # v component

            # Stack u and v as channels: shape (nx*ny*nz, 2)
            # Then reshape to (ny, nx, 2) for a 2D field
            uv_combined = np.stack([u_np, v_np], axis=-1)  # (nx*ny*nz, 2)
            uv_reshaped = uv_combined.reshape((dims[1], dims[0], 2))  # (ny, nx, 2)

            all_data.append(uv_reshaped)

        # Concatenate all data: shape (num_files, ny, nx, 2)
        data_matrix = np.array(all_data)

        data_matrix = np.transpose(
            data_matrix, (0, 3, 1, 2)
        )  # Now shape is (num_files, channels=2, ny, nx)
        print(f"Combined data shape: {data_matrix.shape}")

        return np.array(data_matrix)

    """
    Normalize label and input data to [0, 1] range and generate min_max yaml file
    Parameters:
    ----------
    label_data : np.array
        Label data with shape (num_samples, channels=2, ny, nx)
    input_data : np.array
        Input data with shape (num_samples, channels=1, ny, nx)
    Returns:
    -------
    data_scaled : np.array
        Scaled label data with shape (num_samples, channels=2, ny, nx)
    input_scaled : np.array
        Scaled input data with shape (num_samples, channels=1, ny, nx)
    """
    def normalizeData(self, label_data, input_data):
        u_min = label_data[:, 0, :, :].min()
        u_max = label_data[:, 0, :, :].max()
        v_min = label_data[:, 1, :, :].min()
        v_max = label_data[:, 1, :, :].max()
        u_input_min = input_data[:, 0, :, :].min()
        u_input_max = input_data[:, 0, :, :].max()
        print(f"Input U component original range: [{u_input_min}, {u_input_max}]")
        print(f"U component original range: [{u_min}, {u_max}]")
        print(f"V component original range: [{v_min}, {v_max}]")

        # scale labels:
        data_scaled = np.empty_like(label_data)
        data_scaled[:, 0, :, :] = (label_data[:, 0, :, :] - u_min) / (u_max - u_min)
        data_scaled[:, 1, :, :] = (label_data[:, 1, :, :] - v_min) / (v_max - v_min)

        print(
            f"U component scaled to range: [{data_scaled[:,0,:,:].min()}, {data_scaled[:,0,:,:].max()}]"
        )
        print(
            f"V component scaled to range: [{data_scaled[:,1,:,:].min()}, {data_scaled[:,1,:,:].max()}]"
        )

        # scale inputs:
        input_scaled = np.empty_like(input_data)
        input_scaled[:, 0, :, :] = (input_data[:, 0, :, :] - u_input_min) / (
            u_input_max - u_input_min
        )

        self.generateYAML(u_min, u_max, v_min, v_max, u_input_min, u_input_max)

        return data_scaled, input_scaled

    """
    Split data into train, test and validation sets and save as torch datasets
    Parameters:
    ----------
    data_scaled : np.array
        Scaled label data with shape (num_samples, channels=2, ny, nx)
    input_scaled : np.array
        Scaled input data with shape (num_samples, channels=1, ny, nx)
    train_ratio : float
        Ratio of training data
    val_ratio : float
        Ratio of validation data
    """
    def train_test_val_split(self, data_scaled, input_scaled, train_ratio=0.8, test_ratio=0.1):
        indices = np.arange(input_scaled.shape[0])
        np.random.seed(42)
        np.random.shuffle(indices)

        train_count = int(len(indices) * train_ratio)
        test_count = int(len(indices) * test_ratio)
        train_indices = indices[:train_count]
        val_indices = indices[train_count + test_count:]
        test_indices = indices[train_count: train_count + test_count]

        train_inputs = input_scaled[train_indices]
        train_labels = data_scaled[train_indices]

        val_inputs = input_scaled[val_indices]
        val_labels = data_scaled[val_indices]

        test_inputs = input_scaled[test_indices]
        test_labels = data_scaled[test_indices]

        # First split: separate out test set (10%)
        # train_val_input, test_inputs, train_val_labels, test_labels = train_test_split(
        #     input_scaled, data_scaled, test_size=0.1, random_state=42, shuffle=True
        # )

        # # Second split: separate training (80%) and validation (10%) from the remaining 90%
        # train_inputs, val_inputs, train_labels, val_labels = train_test_split(
        #     train_val_input, train_val_labels, test_size=0.1111, random_state=42, shuffle=True
        # )

        # Save datasets as torch datasets
        self.saveAsTorchDataset(test_inputs, test_labels, "test")
        self.saveAsTorchDataset(val_inputs, val_labels, "validation")
        self.saveAsTorchDataset(train_inputs, train_labels, "train")

    def saveAsTorchDataset(self, inputs, labels, name):
        torch_dataset = torch.utils.data.TensorDataset(
            torch.FloatTensor(inputs),
            torch.FloatTensor(labels)
        )
        dir_path = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.join(dir_path, f"{self.dataset_name}_{name}_dataset.pt")
        torch.save(torch_dataset, save_path)

    def generateInputData(
        self, velocities_range=(0.5, 1.5), num_samples=101, grid_shape=(21, 21)
    ):
        input = np.zeros([num_samples, 1, grid_shape[0], grid_shape[1]])
        print(f"Input data shape: {input.shape}")

        velocities = np.linspace(
            velocities_range[0], velocities_range[1], num=num_samples
        )
        for i in range(num_samples):
            n = grid_shape[0]
            input[i, 0,  n - 1, 1:n - 1] = velocities[i]

        return input

    def generateYAML(self, u_min, u_max, v_min, v_max, input_min, input_max):
        min_max = {
            "inputs": {"u": {"max": float(input_max), "min": float(input_min)}},
            "labels": {
                "u": {"max": float(u_max), "min": float(u_min)},
                "v": {"max": float(v_max), "min": float(v_min)},
            },
        }

        # Save to yaml file
        dir_path = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.join(dir_path, f"{self.dataset_name}_min_max.yaml")
        with open(save_path, "w") as f:
            yaml.dump(min_max, f, default_flow_style=False)

        print("min_max.yaml file created successfully!")
        print("\nContents:")
        print(yaml.dump(min_max, default_flow_style=False))
