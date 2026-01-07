import torch.nn as nn


class MyModel(nn.Module):
    def __init__(self, in_channels: int, out_channels: int, activation_function: nn.Module):
        super().__init__()
        self.conv1 = nn.Conv2d(in_channels=in_channels, out_channels=16, kernel_size=7, padding="same", stride=1)
        self.af1 = activation_function
        self.conv2 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1)
        self.af2 = activation_function
        self.conv3 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1)
        self.af3 = activation_function
        self.conv4 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1)
        self.af4 = activation_function
        self.conv5 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1)
        self.af5 = activation_function
        self.conv6 = nn.Conv2d(in_channels=16, out_channels=16, kernel_size=7, padding="same", stride=1)
        self.af6 = activation_function
        self.conv7 = nn.Conv2d(in_channels=16, out_channels=out_channels, kernel_size=7, padding="same", stride=1)

    def forward(self, x):
        x = self.af1(self.conv1(x))
        x = self.af2(self.conv2(x))
        x = self.af3(self.conv3(x))
        x = self.af4(self.conv4(x))
        x = self.af5(self.conv5(x))
        x = self.af6(self.conv6(x))
        x = self.conv7(x)
        return x


def init_my_model(in_channels=1, out_channels=2, activation_function=nn.SiLU()):
    return MyModel(in_channels=in_channels, out_channels=out_channels, activation_function=activation_function)
