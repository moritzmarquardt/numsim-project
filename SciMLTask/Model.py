import torch.nn as nn
import torch
import matplotlib.pyplot as plt
import os

# torch.set_default_dtype(torch.float64)


class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        # First layer
        self.conv1 = nn.Conv2d(
            in_channels=1, out_channels=16, kernel_size=7, padding="same", stride=1
        )
        self.relu1 = nn.ReLU()

        # Five hidden layers
        self.conv2 = nn.Conv2d(
            in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1
        )
        self.relu2 = nn.ReLU()

        self.conv3 = nn.Conv2d(
            in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1
        )
        self.relu3 = nn.ReLU()

        self.conv4 = nn.Conv2d(
            in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1
        )
        self.relu4 = nn.ReLU()

        self.conv5 = nn.Conv2d(
            in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1
        )
        self.relu5 = nn.ReLU()

        self.conv6 = nn.Conv2d(
            in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1
        )
        self.relu6 = nn.ReLU()

        # Final layer
        self.conv7 = nn.Conv2d(
            in_channels=16, out_channels=2, kernel_size=7, padding="same", stride=1
        )

    def forward(self, x):
        x = self.relu1(self.conv1(x))
        x = self.relu2(self.conv2(x))
        x = self.relu3(self.conv3(x))
        x = self.relu4(self.conv4(x))
        x = self.relu5(self.conv5(x))
        x = self.relu6(self.conv6(x))
        x = self.conv7(x)
        return x

    """ Load data into the model
    Parameters:
    ----------
    train_labels : Tensor
        Labels for training data with diomension (num_sim in train_set, channels = 2, field_width, field_height)
    train_input : Tensor
        Input data for training with diomension (num_sim, channels = 1, field_width, field_height)
    validation_labels : Tensor
        Labels for validation data with diomension (num_sim in validation_set, channels = 2, field_width, field_height)
    validation_input : Tensor
        Input data for validation
    test_labels : Tensor
        Labels for test data with diomension (num_sim in test_set, channels = 2, field_width, field_height)
    test_input : Tensor
        Input data for test data
    """

    def load_data(
        self,
        dataset_name,
    ):
        self.dataset_name = dataset_name
        dir_path = os.path.dirname(os.path.abspath(__file__)) + "/out_ml"
        print(f"Loading datasets from {dir_path}")
        # Use safe globals to allow loading TensorDataset
        with torch.serialization.safe_globals([torch.utils.data.dataset.TensorDataset]):
            self.train_dataset = torch.load(f"{dir_path}/{self.dataset_name}_train_dataset.pt")
            print(f"shape of loaded train inputs: {self.train_dataset.tensors[0].shape} and labels: {self.train_dataset.tensors[1].shape}")
            self.validation_dataset = torch.load(f"{dir_path}/{self.dataset_name}_validation_dataset.pt")
            self.test_dataset = torch.load(f"{dir_path}/{self.dataset_name}_test_dataset.pt")

        self.train_loader = torch.utils.data.DataLoader(
            self.train_dataset, batch_size=16, shuffle=True
        )

    def initialize_model(self):
        self.criterion = nn.MSELoss()
        self.optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        # ensure all parameters are float64 (double)
        # self.double()
        self.to(self.device)

    def train_model(self, epochs=4000):
        train_losses = []
        val_losses = []
        for epoch in range(epochs):
            self.train()
            train_loss = 0.0

            for batch_inputs, batch_labels in self.train_loader:
                inputs, labels = batch_inputs.to(self.device), batch_labels.to(self.device)

                self.optimizer.zero_grad()
                outputs = self(inputs)
                loss = self.criterion(outputs, labels)
                loss.backward()
                self.optimizer.step()

                train_loss += loss.item() * inputs.size(0)

            train_loss /= len(self.train_dataset)

            # Validation
            self.eval()
            with torch.no_grad():
                val_input_device = self.validation_dataset.tensors[0].to(self.device)
                val_labels_device = self.validation_dataset.tensors[1].to(self.device)
                val_outputs = self(val_input_device)
                val_loss = self.criterion(val_outputs, val_labels_device).item()

            if (epoch + 1) % 100 == 0:
                print(
                    f"Epoch [{epoch+1}/{epochs}], Train Loss: {train_loss:.6f}, Val Loss: {val_loss:.6f}"
                )
            train_losses.append(train_loss)
            val_losses.append(val_loss)

        print("Training completed!")
        self.train_losses = train_losses
        self.val_losses = val_losses

    def save_model(self, filename="model.pt"):
        dir_path = os.path.dirname(os.path.abspath(__file__)) + "/out_ml"
        print(f"Saving model to {dir_path}/{filename}")
        torch.save(self.state_dict(), f"{dir_path}/{filename}")

    def load_model(self, filename="model.pt"):
        dir_path = os.path.dirname(os.path.abspath(__file__)) + "/out_ml"
        print(f"Loading model from {dir_path}/{filename}")
        self.load_state_dict(torch.load(f"{dir_path}/{filename}"))
        self.to(self.device)

    def plot_loss(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.train_losses, label="Train Loss")
        plt.plot(self.val_losses, label="Validation Loss")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Epochs")
        plt.ylabel("Loss")
        plt.title("Training and Validation Loss over Epochs")
        plt.legend()
        plt.show()

    def evaluate_model(self):
        self.eval()
        with torch.no_grad():
            test_input_device = self.test_dataset.tensors[0].to(self.device)
            test_labels_device = self.test_dataset.tensors[1].to(self.device)
            test_outputs = self(test_input_device)
            test_loss = self.criterion(test_outputs, test_labels_device).item()

        print(f"Test Loss: {test_loss:.6f}")
        return test_loss

    """ Predict using the trained model
    Parameters:
    ----------
    input_data : Tensor
        Input data for prediction with diomension (num_sim, channels = 1, field_width, field_height)
    Returns:
    -------
    outputs : Tensor
        Predicted output data with diomension (num_sim, channels = 2, field_width, field_height)
    """

    def predict(self, input_data):
        self.eval()
        with torch.no_grad():
            input_device = input_data.to(self.device)
            outputs = self(input_device)
        return outputs

    def plot_test_results(self, num_samples=3):
        self.eval()
        with torch.no_grad():
            test_inputs = self.test_dataset.tensors[0][:num_samples].to(self.device)
            test_labels = self.test_dataset.tensors[1][:num_samples].to(self.device)
            test_outputs = self(test_inputs)

        test_inputs = test_inputs.cpu().numpy()
        test_labels = test_labels.cpu().numpy()
        test_outputs = test_outputs.cpu().numpy()

        for i in range(num_samples):
            plt.figure(figsize=(15, 5))

            plt.subplot(1, 3, 1)
            plt.title("Input")
            plt.imshow(test_inputs[i, 0, :, :], cmap="viridis")
            plt.colorbar()

            plt.subplot(1, 3, 2)
            plt.title("Ground Truth")
            plt.imshow(test_labels[i, 0, :, :], cmap="viridis")
            plt.colorbar()

            plt.subplot(1, 3, 3)
            plt.title("Prediction")
            plt.imshow(test_outputs[i, 0, :, :], cmap="viridis")
            plt.colorbar()

            plt.show()
