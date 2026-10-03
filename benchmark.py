import  numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from pathlib import Path
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
PROJECT_DIR = Path(__file__).resolve().parent.parent
GCSR_DIR = PROJECT_DIR / "GCSR" / "GCSR"
block_file = GCSR_DIR / 'n300/n300_blocks.txt'
net_file = GCSR_DIR / 'n300/n300_nets.txt'
pl_file = GCSR_DIR / 'n300/n300.pl'
#block_file = GCSR_DIR / 'n10/n10_blocks.txt'
#net_file = GCSR_DIR / 'n10/n10_nets.txt'
#pl_file = GCSR_DIR / 'n10/n10.pl'
class Blocks:
    def __init__(self, block_file):
        self.file_name = block_file

    def read_block(self):
        # 用来存储块信息的字典
        blocks = []
        block_order = []  # 保持块的读取顺序
        block_dict = {}
        sb_count = 0
        p_count = 0
        area = 0
        with open( self.file_name, 'r') as f:
            lines = (line.strip() for line in f if line.strip() and not line.startswith(('#', 'Num', 'UCSC')))
            for line in lines:
                parts = line.replace(", ", ",").split()
                block_name = parts[0]                                                       
                if 'hardrectilinear' in line:
                    vertices = [tuple(map(int, v.strip('()').split(','))) for v in parts[3:3 + int(parts[2])]]
                    xs, ys = zip(*vertices)
                    width = max(xs) - min(xs)
                    height = max(ys) - min(ys)  
                    block_dict[block_name] = sb_count 
                    sb_count += 1
                    blocks.append([ width, height])
                    area = area + width * height
                elif 'terminal' in line:
                    width = height = 0
                    block_dict[block_name] = sb_count + p_count 
                    p_count += 1
                    blocks.append([width, height])
                #block_order.append(block_name)
        return np.array(blocks), block_order, sb_count, block_dict, area


class Nets:
    def __init__(self, net_file, blocks, block_dict):
         self.file_name = net_file
         blocks = blocks
         block_dict = block_dict

    def read_net(self):
        nets_list = []
        nets_find = []
        max_net_degree = 0  
        num_nets = 0
        current_net = None
        with open( self.file_name, 'r') as f:
            lines = (line.strip() for line in f if line.strip() and not line.startswith('#'))
            for line in lines:
                if line.startswith('#') or not line.strip():
                    continue
                if 'NumNets' in line:
                    num_nets = int(line.split(':')[-1])
                    continue
                if 'NetDegree' in line:
                    if current_net:
                        nets_list.append(current_net)
                        if any(cell.startswith('sb') for cell in current_net['cells']):
                            nets_find.append(current_net)
                    net_degree = int(line.split(':')[-1])
                    current_net = {'degree': net_degree, 'cells': []}
                    max_net_degree = max(max_net_degree, net_degree)  
                    continue
                block_name = line.split()[0]
                if current_net is not None:
                    current_net['cells'].append(block_name)
            if current_net:
                nets_list.append(current_net)
                if any(cell.startswith('sb') for cell in current_net['cells']):
                    nets_find.append(current_net)
        nets = -np.ones((len(nets_list), max_net_degree), dtype=int)
        for i, net in enumerate(nets_list):
            for j, block_name in enumerate(net['cells']):
                if block_name in  block_dict:
                    block_id =  block_dict[block_name]
                nets[i, j] = block_id
        return nets,num_nets,nets_list,nets_find
    
class P_Position:
    def __init__(self, pl_file,blocks):
        self.file_name = pl_file
        self.blocks = blocks

    def read_pl(self):
        positions = []
        i = 0
        with open(self.file_name, 'r') as file:
            for line in file:
                if line.startswith('sb') :
                    _, x, y = line.split()
                    positions.append([float(x)+self.blocks[i,0]/2, float(y)+self.blocks[i,1]/2])  # 确保x和y是浮点数
                    i += 1
                if line.startswith('p') :
                    _, x, y = line.split()
                    positions.append([float(x), float(y)]) 
        return np.array(positions)  
    
class P_Position:
    def __init__(self, pl_file,blocks):
        self.file_name = pl_file
        self.blocks = blocks

    def read_pl(self):
        positions = []
        i = 0
        with open(self.file_name, 'r') as file:
            for line in file:
                if line.startswith('sb') :
                    _, x, y = line.split()
                    #positions.append([float(x)+self.blocks[i,0]/2+ 100, float(y)+self.blocks[i,1]/2+ 100]) 
                    positions.append([float(x)+self.blocks[i,0]/2, float(y)+self.blocks[i,1]/2]) 
                    i += 1
                if line.startswith('p') :
                    _, x, y = line.split()
                    positions.append([float(x), float(y)]) 
        return np.array(positions)  

class loss_function: 
    def __init__(self,  W, H):
         W = W
         H = H
        
    def calculate_distance_loss(self, nets,  positions):
        #self.valid_indices = batch_idx
        distance_loss =  0.0
        for net in nets:     
            valid_cells = net[net!= -1]  
            #valid_cells = cells if any(cell in self.valid_indices for cell in cells) else [] # 过滤出有效的细胞
            max_distance1 = 0
            max_distance2 = 0
            if len(valid_cells) >=2:
                for i in range(len(valid_cells)):
                    for j in range(i + 1, len(valid_cells)):
                        cell_i = valid_cells[i]
                        cell_j = valid_cells[j]
                        xi, yi = positions[cell_i]
                        xj, yj = positions[cell_j]
                        distance1 = abs(xi - xj) 
                        distance2 = abs(yi - yj)
                        if distance1 > max_distance1:
                            max_distance1 = distance1
                        if distance2 > max_distance2:
                            max_distance2 = distance2
            distance_loss = distance_loss + max_distance1 + max_distance2
        return distance_loss
        
    
    def calculate_overlap_loss0(self, batch_idx, positions, blocks):
        overlap_loss =  0.0
        for i in range(len(batch_idx)):
            for j in range(i + 1, len(batch_idx)):
                idx = batch_idx[i]
                idy = batch_idx[j]
                wi, hi = blocks[idx]
                wj, hj = blocks[idy]
                xi, yi = positions[idx]
                xj, yj = positions[idy]
                lx =  max(- abs(xi - xj) + (wi + wj) / 2,0) 
                ly =  max(- abs(yi - yj) + (hi + hj) / 2,0) 
                overlap_loss = overlap_loss +  lx * ly 
        return overlap_loss
 
 
blocks,  block_order,  sb_count,  block_dict, area = Blocks(block_file).read_block()
block1 =  blocks[ blocks[:, 0] > 0]
nets, num_nets, nets_list ,  nets_find= Nets(net_file, blocks, block_dict).read_net()
#print("blocks[1,0]",blocks[1,0].type())
num_blocks = len( blocks)
full_idx = np.arange( sb_count)
net_idx = np.arange( num_nets)
positions = P_Position( pl_file,block1).read_pl()
np.save('benchmark.npy', positions)
Loss = loss_function(  800, 800)
wirelength = Loss.calculate_distance_loss( nets,  positions)
overlap = Loss.calculate_overlap_loss0( full_idx, positions,  block1)
plt.close(plt.gcf())
fig, ax = plt.subplots(figsize=(10, 8))
for i in range( sb_count):
    width, height = blocks[i]
    x_center, y_center = positions[i]
    left_x = x_center - width / 2
    left_y = y_center - height / 2
    rect = patches.Rectangle((left_x, left_y), width, height, linewidth=1, edgecolor='blue',  facecolor='skyblue', alpha=0.5)
    ax.add_patch(rect)
    ax.text(x_center, y_center, '' , ha='center', va='center', color='blue')
title_text = f'{"benchmark"}: hpwl = {wirelength:.2f} overlap = {overlap:.2f}'
ax.set_title(f'{title_text}', fontsize=18, color='black', fontweight='bold')
ax.set_xlim(0, 800)
ax.set_ylim(0, 800)
ax.set_xlabel('Width', fontsize=18)  # 修改 x 轴标签字体大小
ax.set_ylabel('Height', fontsize=18)  # 修改 y 轴标签字体大小
ax.tick_params(axis='both', which='major', labelsize=18)  # 修改主刻度标签的字体大小
ax.tick_params(axis='both', which='minor', labelsize=18)  # 修改次刻度标签的字体大小
ax.set_xlabel('Width')
ax.set_ylabel('Height')
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
filename = f'{"benchmark"}.eps'
plt.savefig(filename, transparent=False, bbox_inches='tight')
plt.show()
plt.close()
       


    
