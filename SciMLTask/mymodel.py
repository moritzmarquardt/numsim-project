from typing import OrderedDict
import torch


def init_my_model(
    in_channels: int = 1, out_channels: int = 2, activation_function=torch.nn.ReLU()
):
    model = torch.nn.Sequential(
        OrderedDict(
            [
                (
                    "conv1",
                    torch.nn.Conv2d(
                        in_channels=in_channels,
                        out_channels=16,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
                ("relu1", activation_function),
                (
                    "conv2",
                    torch.nn.Conv2d(
                        in_channels=16,
                        out_channels=16,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
                ("relu2", activation_function),
                (
                    "conv3",
                    torch.nn.Conv2d(
                        in_channels=16,
                        out_channels=16,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
                ("relu3", activation_function),
                (
                    "conv4",
                    torch.nn.Conv2d(
                        in_channels=16,
                        out_channels=16,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
                ("relu4", activation_function),
                (
                    "conv5",
                    torch.nn.Conv2d(
                        in_channels=16,
                        out_channels=16,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
                ("relu5", activation_function),
                (
                    "conv6",
                    torch.nn.Conv2d(
                        in_channels=16,
                        out_channels=16,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
                ("relu6", activation_function),
                (
                    "conv7",
                    torch.nn.Conv2d(
                        in_channels=16,
                        out_channels=out_channels,
                        kernel_size=7,
                        padding="same",
                        stride=1,
                    ),
                ),
            ]
        )
    )
    return model
