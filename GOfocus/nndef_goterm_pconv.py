import torch
import torch.nn as nn
import torch.nn.functional as F
from math import sqrt, log
import sys

class Maxout1d(nn.Module):

    def __init__(self, in_channels, out_channels, pool_size, kernel_size=1, dilation=1):
        super(Maxout1d, self).__init__()
        self.in_channels, self.out_channels, self.pool_size = in_channels, out_channels, pool_size
        self.lin = nn.Conv1d(in_channels=in_channels, out_channels=out_channels * pool_size, kernel_size=kernel_size, dilation=dilation, padding=dilation*(kernel_size-1)//2)
        self.norm = nn.InstanceNorm1d(out_channels, affine=True)

    def forward(self, inputs):
        x = self.lin(inputs)

        N, C, W = x.size()

        x = x.view(N, C//self.pool_size, self.pool_size, W)
        x = x.max(dim=2)[0]
        x = self.norm(x)

        return x


class CSE(nn.Module):
    def __init__(self, width, reduction=16):
        super(CSE, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            nn.Linear(width, width // reduction, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(width // reduction, width, bias=False),
            nn.Sigmoid()
        )

    def forward(self, x):
        N, C, _ = x.size()
        y = self.avg_pool(x).view(N, C)
        y = self.fc(y).view(N, C, 1)

        return x * y.expand_as(x)


# ResNet Module
class ResNet_Block(nn.Module):
    def __init__(self,width,fsize,dilv):
        super(ResNet_Block, self).__init__()
        self.dropout = nn.Dropout(p=0.1)
        self.molayer = Maxout1d(in_channels=width, out_channels=width, pool_size=2, kernel_size=fsize, dilation=dilv)
        self.cse = CSE(width, 16)

    def forward(self, x):

        residual = x
        x = self.dropout(x)
        x = self.molayer(x)
        x = self.cse(x)
        out = x + residual
        
        return out


# ResNet Module
class ResNet(nn.Module):
    def __init__(self,width):
        super(ResNet, self).__init__()

        layers = []

        for rep in range(4):
            for fsize,dilv in [(3,1), (3,2), (3,4), (3,8)]:
                if fsize > 0:
                    layer = ResNet_Block(width, fsize, dilv)
                    layers.append(layer)
                    
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        
        out = self.model(x)

        return out


# Main Net Module
class PoolNet(nn.Module):
    def __init__(self,width,pwidth,nclass):
        super(PoolNet, self).__init__()

        self.width = width
        self.pwidth = pwidth

        self.embed = nn.Embedding(22, width, 21)

        self.resnet = ResNet(width)

        layers = []

        layers.append(nn.AdaptiveMaxPool1d(pwidth))
        layers.append(nn.Flatten())
        layers.append(nn.Dropout(0.1))
        layers.append(nn.Linear(width * pwidth, width * 4))
        layers.append(nn.ReLU(inplace=True))
        layers.append(nn.Dropout(0.1))
        layers.append(nn.Linear(width * 4, nclass))
        #layers.append(nn.LayerNorm(nclass, elementwise_affine=False))

        self.outlayers = nn.Sequential(*layers)

    def forward(self, x):

        x = self.embed(x)

        x = self.resnet(x.transpose(1,2))

        out = self.outlayers(x)

        return out
