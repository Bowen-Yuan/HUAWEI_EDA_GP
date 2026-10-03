import torch
from torch import nn, optim
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time
import math
import os
from pathlib import Path
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
PROJECT_DIR = Path(__file__).resolve().parent.parent
GCSR_DIR = PROJECT_DIR / "GCSR" / "GCSR"

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
        with open(self.file_name, 'r') as f:
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
        self.blocks = blocks
        self.block_dict = block_dict

    def read_net(self):
        nets_list = []
        nets_find = []
        max_net_degree = 0  
        num_nets = 0
        current_net = None
        with open(self.file_name, 'r') as f:
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
                if block_name in self.block_dict:
                    block_id = self.block_dict[block_name]
                nets[i, j] = block_id
        return nets,num_nets,nets_list,nets_find
    
class P_Position(nn.Module):
    def __init__(self, blocks,sb_count,plfile):
        super(P_Position, self).__init__()
        self.num_blocks = len(blocks)
        num_p_blocks = self.num_blocks - sb_count
        self.fixed_positions = torch.empty((num_p_blocks, 2), dtype=torch.float32, requires_grad=False)
        p_index = 0
        with open(plfile, 'r') as file:
            for line in file:
                if line.startswith('p'):
                    _, x, y = line.split()
                    self.fixed_positions[p_index] = torch.tensor([float(x), float(y)], dtype=torch.float32)
                    p_index += 1

class Position(nn.Module):
    def __init__(self, sb_count, W, H):
        super(Position, self).__init__()
        self.positions = nn.Parameter(torch.empty((sb_count, 2), dtype=torch.float32))

        with torch.no_grad():
            nn.init.uniform_(self.positions[:, 0], 0.15 * W, 0.85 * W)
            nn.init.uniform_(self.positions[:, 1], 0.15 * H, 0.85 * H)
            
    def forward(self):
        return self.positions
    
class loss_function: 
    def __init__(self, gamma1, gamma2, W, H):
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.W = W
        self.H = H
        

    def hat_function(self, x, y, r, t):
        abs_x = torch.abs(x)
        abs_y = torch.abs(y)
        outside = (abs_x > r) | (abs_y > t)
        x_dominant = (abs_y / t) <= (abs_x / r)
        result = torch.where(outside, 
                            torch.tensor(0.0, device=x.device, dtype=x.dtype), 
                            torch.where(x_dominant, 
                                        1 - abs_x / r, 
                                        1 - abs_y / t))
        return result
    

    def calculate_distance_loss(self, nets, batch_idx, positions):
        #print("batch_idx",batch_idx)
        batch_cells = torch.tensor(nets[batch_idx], dtype=torch.long, device=positions.device)  # [batch_size, net_degree]
        valid_mask = (batch_cells != -1)  # [batch_size, net_degree]
        valid_cells = torch.where(valid_mask, batch_cells, torch.tensor(0, dtype=torch.long))
        all_positions = positions[valid_cells]  # [batch_size, net_degree, 2]
        masked_positions = torch.where(valid_mask.unsqueeze(-1), all_positions, torch.tensor(float('inf'), dtype=torch.float32))
        min_pos = torch.min(masked_positions, dim=1)[0]  # [batch_size, 2] 
        masked_positions = torch.where(valid_mask.unsqueeze(-1), all_positions, torch.tensor(float('-inf'), dtype=torch.float32))
        max_pos = torch.max(masked_positions, dim=1)[0]  # [batch_size, 2]
        distance_loss = (max_pos - min_pos).sum(dim=1)  # [batch_size] 
        total_distance_loss = distance_loss.sum()
        return total_distance_loss
        
    
    def calculate_boundary_loss(self, batch_idx, positions, blocks):
        batch_positions = positions[batch_idx]
        batch_blocks = torch.tensor(blocks[batch_idx], dtype=torch.float, device=batch_positions.device)
        w_i = batch_blocks[:, 0] / 2
        h_i = batch_blocks[:, 1] / 2
        left_loss = torch.relu(-(batch_positions[:, 0] - w_i))
        right_loss = torch.relu(batch_positions[:, 0] - (self.W - w_i))
        bottom_loss = torch.relu(-(batch_positions[:, 1] - h_i))
        top_loss = torch.relu(batch_positions[:, 1] - (self.H - h_i))
        total_loss = (left_loss + right_loss + bottom_loss + top_loss)
        gamma_weighted_loss = self.gamma1[batch_idx] * total_loss
        boundary_loss = torch.sum(gamma_weighted_loss)
        return boundary_loss
    
    def calculate_boundary_loss1(self, batch_idx, positions, blocks):
        batch_positions = positions[batch_idx]
        batch_blocks = torch.tensor(blocks[batch_idx], dtype=torch.float, device=batch_positions.device)
        average_position = torch.mean(positions, dim=0)
        r = average_position - batch_positions
        distances = torch.linalg.norm(r, dim=1)
        w_half = batch_blocks[:, 0] / 2
        h_half = batch_blocks[:, 1] / 2
        left_loss = torch.relu(-(batch_positions[:, 0] - w_half))
        right_loss = torch.relu(batch_positions[:, 0] - (self.W - w_half))
        bottom_loss = torch.relu(-(batch_positions[:, 1] - h_half))
        top_loss = torch.relu(batch_positions[:, 1] - (self.H - h_half))
        boundary_loss_terms = left_loss + right_loss + bottom_loss + top_loss
        force_magnitude = 5
        total_loss = self.gamma1[batch_idx] * boundary_loss_terms + force_magnitude * distances
        boundary_loss = torch.sum(total_loss)
        return boundary_loss
    
    def calculate_overlap_loss(self, batch_idx, positions, blocks):
        pos = positions[batch_idx]  
        sizes = torch.tensor(blocks[batch_idx], dtype=torch.float, device=pos.device)  
        diff = pos.unsqueeze(1) - pos.unsqueeze(0)  
        size_sum = sizes.unsqueeze(1) + sizes.unsqueeze(0) 
        rij = size_sum[..., 0] / 2  
        tij = size_sum[..., 1] / 2  
        penalties = self.hat_function(diff[..., 0], diff[..., 1], rij, tij)
        mask = torch.triu(torch.ones_like(penalties), diagonal=1).bool()
        penalties = penalties[mask]  
        N = len(batch_idx)
        gamma_indices = (N - 1) * torch.arange(N).unsqueeze(1) - torch.arange(N).unsqueeze(1) * (torch.arange(N).unsqueeze(1) + 1) // 2 + torch.arange(N).unsqueeze(0) - 1
        gamma_indices = gamma_indices[mask]  
        gamma_values = self.gamma2[gamma_indices]  
        overlap_loss = torch.sum(gamma_values * penalties)
        return overlap_loss

    
    def calculate_overlap_loss0(self, batch_idx, positions, blocks):
        pos = positions[batch_idx]  
        sizes = torch.tensor(blocks[batch_idx], dtype=torch.float, device=pos.device)  
        pos_i = pos.unsqueeze(1)    
        pos_j = pos.unsqueeze(0)   
        sizes_i = sizes.unsqueeze(1)  
        sizes_j = sizes.unsqueeze(0)  
        delta_pos = pos_i - pos_j  
        abs_delta_pos = torch.abs(delta_pos)  
        half_size_sum = (sizes_i + sizes_j) / 2 
        overlap_sizes = torch.relu(-abs_delta_pos + half_size_sum)  
        mask = torch.triu(torch.ones(overlap_sizes.shape[:-1]), diagonal=1).bool()
        overlap_area = overlap_sizes[..., 0] * overlap_sizes[..., 1]  
        overlap_loss = overlap_area[mask].sum() 
        return overlap_loss
 
class TrainModel:
    def __init__(self, model, net_file, block_file, pfile, W, H,  max_epoch, batch_size, max_step):
        self.blocks, self.block_order, self.sb_count, self.block_dict,self.area = Blocks(block_file).read_block()
        self.block1 = self.blocks[self.blocks[:, 0] > 0]
        self.nets,self.num_nets,self.nets_list , self.nets_find= Nets(net_file,self.blocks,self.block_dict).read_net()
        self.model = model
        self.pfile = pfile
        self.W = W      
        self.H = H
        self.max_epoch = max_epoch
        self.batch_size = batch_size
        self.loss_history = []
        self.num_blocks = len(self.blocks)
        self.full_idx = np.arange(self.sb_count)
        self.net_idx = np.arange(self.num_nets)
        p_positions = P_Position(self.blocks, self.sb_count, self.pfile)
        self.fixed_positions = p_positions.fixed_positions
        self.max_setp = max_step
        blocks = torch.tensor(self.block1, dtype=torch.float)
        wi = blocks[:, 0].unsqueeze(1) 
        wj = blocks[:, 0].unsqueeze(0)  
        hi = blocks[:, 1].unsqueeze(1)  
        hj = blocks[:, 1].unsqueeze(0) 
        self.rij_matrix = (wi + wj) / 2
        self.tij_matrix = (hi + hj) / 2
        self.w_half = blocks[:, 0] / 2
        self.h_half = blocks[:, 1] / 2
        self.triu_indices = torch.triu_indices(len(self.full_idx), len(self.full_idx), offset=1)
        self.rij = self.rij_matrix[self.triu_indices[0], self.triu_indices[1]]
        self.tij = self.tij_matrix[self.triu_indices[0], self.triu_indices[1]]
        self.y1 = torch.randn(self.sb_count, requires_grad=True) 
        self.y2 = torch.randn(self.sb_count*(self.sb_count-1)//2, requires_grad=True) 
    
    def sample(self, nets_list):
        terminal_to_nets = {}
        for i, net in enumerate(nets_list):
            for terminal in net['cells']:
                if terminal not in terminal_to_nets:
                    terminal_to_nets[terminal] = []
                terminal_to_nets[terminal].append(i)
        n = len(nets_list)
        adjacency_matrix = np.zeros((n, n), dtype=int)
        for nets in terminal_to_nets.values():
            for i in range(len(nets)):
                for j in range(i + 1, len(nets)):
                    adjacency_matrix[nets[i]][nets[j]] = 1
                    adjacency_matrix[nets[j]][nets[i]] = 1
        degrees = adjacency_matrix.sum(axis=0)
        probabilities = np.exp(degrees) / np.exp(degrees).sum()
        return probabilities
    
    def calculate_gamma1(self, positions):
        positions = positions.clone().detach().requires_grad_(True)
        all_positions = torch.cat((positions, self.fixed_positions), dim=0)
        gamma1 = torch.randn(self.sb_count, requires_grad=True) 
        gamma2 = torch.randn(self.sb_count*(self.sb_count-1)//2, requires_grad=True) 
        Loss = loss_function(1, 1, self.W, self.H)
        distance_loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)
        distance_loss.backward()
        dis_grad = positions.grad.clone()
        positions.grad.zero_()
       
        def boundary_loss(positions):
            boundary_loss = (torch.relu(-positions[:, 0] + self.w_half) + 
                            torch.relu(positions[:, 0] - self.W + self.w_half) +
                            torch.relu(-positions[:, 1] + self.h_half) + 
                            torch.relu(positions[:, 1] - self.H + self.h_half))
            return boundary_loss
        jacobian_matrix = torch.autograd.functional.jacobian(boundary_loss, positions)
        divided_result = torch.where(jacobian_matrix != 0,
                             dis_grad.unsqueeze(0) / jacobian_matrix,
                             torch.zeros_like(dis_grad.unsqueeze(0)))
        gamma1 = torch.ceil(divided_result.max(dim=2)[0].max(dim=1)[0])
        positions.grad.zero_()

        def overlap_loss(positions):
            xi = positions[:, 0]
            yi = positions[:, 1]
            xi_xj_matrix = xi.unsqueeze(1) - xi.unsqueeze(0)
            yi_yj_matrix = yi.unsqueeze(1) - yi.unsqueeze(0)
            xi_xj = xi_xj_matrix[self.triu_indices[0], self.triu_indices[1]]
            yi_yj = yi_yj_matrix[self.triu_indices[0], self.triu_indices[1]]     
            overlap = Loss.hat_function(xi_xj,yi_yj,self.rij,self.tij)       
            return overlap
        jacobian_matrix = torch.autograd.functional.jacobian(overlap_loss, positions)
        divided_result = torch.where(jacobian_matrix != 0,
                             dis_grad.unsqueeze(0) / jacobian_matrix,
                             torch.zeros_like(dis_grad.unsqueeze(0)))
        gamma2 = torch.ceil(divided_result.max(dim=2)[0].max(dim=1)[0])
        positions.grad.zero_()

        return gamma1,gamma2
    
    def calculate_gamma(self, positions):
        positions = positions.clone().detach().requires_grad_(True)
        all_positions = torch.cat((positions, self.fixed_positions), dim=0)
        gamma1 = torch.zeros(self.sb_count, requires_grad=True) 
        gamma2 = torch.zeros(self.sb_count*(self.sb_count-1)//2, requires_grad=True) 
        Loss = loss_function(1, 1, self.W, self.H)
        distance_loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)
        distance_loss.backward()
        dis_grad = positions.grad.clone()
        positions.grad.zero_()
        
        boundary_loss = (torch.relu(-positions[:, 0] + self.w_half) + 
                            torch.relu(positions[:, 0] - self.W + self.w_half) +
                            torch.relu(-positions[:, 1] + self.h_half) + 
                            torch.relu(positions[:, 1] - self.H + self.h_half))
        boundary_mask = boundary_loss != 0
        indices = torch.where(boundary_mask)
        def boundary_loss(positions):
            boundary_loss = (torch.relu(-positions[:, 0] + self.w_half) + 
                            torch.relu(positions[:, 0] - self.W + self.w_half) +
                            torch.relu(-positions[:, 1] + self.h_half) + 
                            torch.relu(positions[:, 1] - self.H + self.h_half))
            boundary_mask = boundary_loss != 0
            valid_boundary = torch.masked_select(boundary_loss, boundary_mask)
            return valid_boundary
        if indices[0].numel() > 0:
            jacobian_matrix = torch.autograd.functional.jacobian(boundary_loss, positions)
            divided_result = torch.where(jacobian_matrix != 0,
                                dis_grad.unsqueeze(0) / jacobian_matrix,
                                torch.zeros_like(dis_grad.unsqueeze(0)))
            max_values = torch.ceil(divided_result.max(dim=2)[0].max(dim=1)[0])
            #max_values = divided_result.max(dim=2)[0].max(dim=1)[0]
            new_gamma1 = gamma1.clone()
            new_gamma1[indices[0]] = max_values
            gamma1 = new_gamma1
            positions.grad.zero_()
        
        xi = positions[:, 0]
        yi = positions[:, 1]
        xi_xj_matrix = xi.unsqueeze(1) - xi.unsqueeze(0)
        yi_yj_matrix = yi.unsqueeze(1) - yi.unsqueeze(0)
        xi_xj = xi_xj_matrix[self.triu_indices[0], self.triu_indices[1]]
        yi_yj = yi_yj_matrix[self.triu_indices[0], self.triu_indices[1]] 
        overlap = Loss.hat_function(xi_xj,yi_yj,self.rij,self.tij)
        overlap_mask = overlap != 0
        indices = torch.where(overlap_mask)
        def overlap_loss(positions):       
            xi = positions[:, 0]
            yi = positions[:, 1]
            xi_xj_matrix = xi.unsqueeze(1) - xi.unsqueeze(0)
            yi_yj_matrix = yi.unsqueeze(1) - yi.unsqueeze(0)
            xi_xj = xi_xj_matrix[self.triu_indices[0], self.triu_indices[1]]
            yi_yj = yi_yj_matrix[self.triu_indices[0], self.triu_indices[1]] 
            overlap = Loss.hat_function(xi_xj,yi_yj,self.rij,self.tij)
            overlap_mask = overlap != 0
            valid_overlap = torch.masked_select(overlap, overlap_mask)
            return valid_overlap
        if indices[0].numel() > 0:
            jacobian_matrix = torch.autograd.functional.jacobian(overlap_loss, positions)
            divided_result = torch.where(jacobian_matrix != 0,
                                dis_grad.unsqueeze(0) / jacobian_matrix,
                                torch.zeros_like(dis_grad.unsqueeze(0)))
            max_values = torch.ceil(divided_result.max(dim=2)[0].max(dim=1)[0])
            #max_values = divided_result.max(dim=2)[0].max(dim=1)[0]
            #print("max_values",max_values)
            new_gamma2 = gamma2.clone()
            new_gamma2[indices[0]] = max_values
            gamma2 = new_gamma2
            positions.grad.zero_()
        #print("gamma1,gamma2",gamma1,gamma2)
        return gamma1,gamma2
    
    

    def update_penalty(self, positions): 
        new_positions = positions.clone().detach().requires_grad_(True)
        for i in range(len(self.full_idx)):
            for j in range(len(self.full_idx)):  # 从0开始，包括整个范围
                if i != j: 
                    wi, hi = self.blocks[i]
                    wj, hj = self.blocks[j]
                    xi, yi = new_positions[i]
                    xj, yj = new_positions[j]
                    rij = (wi + wj) / 2
                    tij = (hi + hj) / 2
                    Loss = loss_function(1, 1, self.W, self.H)
                    penalty = Loss.hat_function(xi - xj, yi - yj, rij, tij)
                    if penalty > 0:
                        penalty.backward(retain_graph=True)
                        penalty_grad = new_positions.grad[i].clone()
                        new_positions.grad.zero_() 
                        grad_norm = penalty_grad.norm()
                        penalty_step = torch.tensor([(rij - abs(xi - xj)), (tij - abs(yi - yj))], dtype=penalty_grad.dtype, device=penalty_grad.device)
                        update = (penalty_step * penalty_grad) / grad_norm
                        updated_value = new_positions[i] - update
                        with torch.no_grad():
                            new_positions[i] = updated_value.detach()
                        new_positions.requires_grad_(True) 
                        break
            pos = new_positions[i]
            wi, hi = self.blocks[i]
            boundary = (torch.relu(-pos[0] + wi / 2) + torch.relu(pos[0] - self.W + wi / 2) +
                        torch.relu(-pos[1] + hi / 2) + torch.relu(pos[1] - self.H + hi / 2))
            boundary.backward(retain_graph=True)
            boundary_grad = new_positions.grad[i].clone()
            grad_norm = boundary_grad.norm()
            if grad_norm > 0:
                penalty_step = torch.zeros(2, dtype=boundary_grad.dtype, device=boundary_grad.device)
                if boundary_grad[0] > 0:
                    penalty_step[0] = pos[0]-(self.W - wi / 2)
                else:
                    penalty_step[0] = wi / 2 - pos[0]
                if boundary_grad[1] > 0:
                    penalty_step[1] = pos[1]-(self.H - hi / 2)
                else:
                    penalty_step[1] = hi / 2 - pos[1]
                update = (penalty_step * boundary_grad )/ grad_norm
                updated_value = new_positions[i] -  update
                with torch.no_grad():
                    new_positions[i] = updated_value.detach()
                new_positions.requires_grad_(True)  #
        return new_positions


                
    def draw_blocks(self, blocks, positions, all_positions, name, block_order, gamma1, gamma2):
        plt.close(plt.gcf())
        fig, ax = plt.subplots(figsize=(10, 8))
        for i in range(self.sb_count):
            width, height = blocks[i]
            x_center, y_center = positions[i].detach().numpy()
            left_x = x_center - width / 2
            left_y = y_center - height / 2
            rect = patches.Rectangle((left_x, left_y), width, height, linewidth=1, edgecolor='blue', facecolor='skyblue', alpha=0.5)
            ax.add_patch(rect)
            ax.text(x_center, y_center, '', ha='center', va='center', color='blue')

        loss_fn = loss_function(gamma1, gamma2, self.W, self.H)
        wirelength = loss_fn.calculate_distance_loss(self.nets, self.net_idx, all_positions)
        overlap = loss_fn.calculate_overlap_loss0(self.full_idx, positions, self.blocks)
        title_text = f'{name}: hpwl = {wirelength:.2f} overlap = {overlap:.2f}'
        print('hpwl,overlap', wirelength, overlap)
        ax.set_title(f'{title_text}', fontsize=18, color='black', fontweight='bold')
        ax.set_xlim(0, self.W)
        ax.set_ylim(0, self.H)
        ax.set_xlabel('Width', fontsize=18)  # 修改 x 轴标签字体大小
        ax.set_ylabel('Height', fontsize=18)  # 修改 y 轴标签字体大小
        ax.tick_params(axis='both', which='major', labelsize=18)  # 修改主刻度标签的字体大小
        ax.tick_params(axis='both', which='minor', labelsize=18)  # 修改次刻度标签的字体大小
        plt.gca().set_aspect('equal', adjustable='box')
        plt.draw()
        filename = f'{name}placement.eps'
        plt.savefig(filename, transparent=False, bbox_inches='tight')
        plt.show()
        plt.close()



    def finalize_training(self, start_time, all_positions, positions, name, gamma1, gamma2):
        end_time = time.time()
        total_training_time = end_time - start_time
        print(total_training_time)
        print(f"Final Loss = {self.loss_history[-1]}")
        plt.figure(figsize=(10, 5))
        plt.semilogy(self.loss_history, label='Training Loss')
        #plt.title('Loss Curve during Training')
        plt.title(f'\nTotal Training Time: {total_training_time:.2f} seconds', fontsize=14)
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        plt.legend()
        plt.grid(True)
        filename = f'{name}trainloss.png'
        plt.savefig(filename, transparent=False, bbox_inches='tight')
        plt.show()
        self.draw_blocks(self.blocks, positions, all_positions, name, self.block_order, gamma1, gamma2 )

    def compute_batches(self,positions):
        # 计算当前网格划分下的块批次
        num_cells = self.W // self.grid_size
        batch_idx = {}
        for i, pos in enumerate(positions):
            grid_x = int(pos[0] // self.grid_size)
            grid_y = int(pos[1] // self.grid_size)
            batch_key = grid_y * num_cells + grid_x
            if batch_key not in batch_idx:
                batch_idx[batch_key] = []
            batch_idx[batch_key].append(i)
        return batch_idx
    
    def update_gradient(self, positions, loss, optimizer, epoch):
        loss.backward(retain_graph=True)
        grad = positions.grad.clone()
        eps = 0.2 / ((epoch + 1) ** 2)
        grad_length = torch.norm(grad).item()
        direction = np.random.normal(0, 1, grad.shape)
        direction = torch.tensor(direction, dtype=grad.dtype, device=grad.device)
        new_direction = grad + eps * grad_length * direction
        positions.grad = new_direction
        optimizer.step()
        positions.grad.zero_()

    def train_RBSM(self, optimizer_distance, optimizer_penalty):
        start_time = time.time()
        lr_logs1, lr_logs2 = [], []
        positions = self.model()
        gamma1 = torch.full((self.sb_count,), 1000.0, requires_grad=True)
        gamma2 = torch.full((self.sb_count*(self.sb_count-1)//2,), 1000.0, requires_grad=True)

        scheduler_distance = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer_distance, T_max=self.max_setp)
        scheduler_penalty = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer_penalty, T_max=self.max_setp)
        #scheduler_distance = torch.optim.lr_scheduler.ExponentialLR(optimizer_distance, gamma=0.95)
        #scheduler_penalty = torch.optim.lr_scheduler.ExponentialLR(optimizer_penalty, gamma=0.95)
        #switch_step = int(self.max_setp * 0.9)  # 前 50% 纯余弦，后 50% 余弦+衰减

        # scheduler_distance = torch.optim.lr_scheduler.LambdaLR(
        #     optimizer_distance,
        #     lr_lambda=lambda k: (
        #         0.5*(1 + math.cos(math.pi*k/self.max_setp))   # 余弦部分
        #         if k < switch_step 
        #         else 0.5*(1 + math.cos(math.pi*k/self.max_setp)) / math.sqrt(k+1)  # 后期乘1/sqrt(k)
        #     )
        # )

        # scheduler_penalty = torch.optim.lr_scheduler.LambdaLR(
        #     optimizer_penalty,
        #     lr_lambda=lambda k: (
        #         (1 + math.cos(math.pi*k/self.max_setp)) 
        #         if k < switch_step
        #         else (1 + math.cos(math.pi*k/self.max_setp)) / math.sqrt(k+1)
        #     )
        # )
    

        batch_size = math.ceil(self.num_nets / 25)
        probabilities = self.sample(self.nets_find)
        k = 1
        last_lengths = []
        convergence_threshold = 0.0001  # 0.01%
        self.overlap_threshold = self.area * 0.01
        print(f"overlap_threshold={self.overlap_threshold}")
        max_convergences = 3
        last_overlaps = []
        stop_training = False 
        while True:
            new_gamma1, new_gamma2 = self.calculate_gamma(positions)
            gamma1 = torch.max(gamma1,  100 *new_gamma1)
            gamma2 = torch.max(gamma2,  100 *new_gamma2)        
            #gamma1 =( 0.5*gamma1 + 0.5*100 *new_gamma1  )
            #gamma2 =( 0.5*gamma2 + 0.5*100 *new_gamma2  )            
            Loss = loss_function(gamma1, gamma2, self.W, self.H)
            #print(gamma1,gamma2)
            for step in range(25):
                optimizer_distance.zero_grad()
                if k < 20:
                    batch_idx = np.random.choice(self.num_nets, size=batch_size, replace=False, p=probabilities)
                else:
                    batch_idx = self.net_idx
                #batch_idx = self.net_idx
                #print("batch_idx",batch_idx)
                all_positions = torch.cat((positions, self.fixed_positions), dim=0)
                distance_loss = Loss.calculate_distance_loss(self.nets, batch_idx, all_positions)
                #distance_loss.backward(retain_graph=True)
                optimizer_distance.step()
                self.update_gradient(positions, distance_loss, optimizer_distance, k)
                penalty_loss = Loss.calculate_boundary_loss1(self.full_idx, positions, self.block1) + \
                                Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
                #penalty_loss.backward(retain_graph=True)
                optimizer_penalty.step()
                self.update_gradient(positions, penalty_loss, optimizer_penalty, k)
                loss = distance_loss + penalty_loss
                overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1)    
                if step % 5 == 0:
                    with torch.no_grad():
                        current_lr1 = scheduler_distance.get_last_lr()[0]
                        current_lr2 = scheduler_penalty.get_last_lr()[0]
                        lr_logs1.append(current_lr1)
                        lr_logs2.append(current_lr2)
                        self.loss_history.append(loss.item())
                        #print(f"Epoch {(k-1)*25+step+1}:  Loss = {loss.item():.4f}\t diatance={distance_loss}\t overlap ={overlap}\t  distance_lr={current_lr1}\t penalty_lr={current_lr2 }")
                        #print(f"Epoch {k}:  Loss = {loss.item():.4f}\t distance_lr={current_lr1}\t penalty_lr={current_lr2 }")
                # Check convergence based on the last 5 total losses
                if len(last_lengths) >= max_convergences and all(abs(last_lengths[i] - last_lengths[i+1]) < convergence_threshold * last_lengths[i] for i in range(len(last_lengths)-1)):
                    if overlap < self.overlap_threshold:  # Check overlap condition only in the last step
                        print("Training stopped due to minimal change in wire length and acceptable overlap.")
                        stop_training = True
                        break
                last_lengths.append(distance_loss.item())
                last_overlaps.append(overlap)
                if len(last_lengths) > max_convergences:
                    last_lengths.pop(0)  # keep only the last few measurements
            k += 1
            scheduler_distance.step()
            scheduler_penalty.step()
            if stop_training or k == self.max_setp:
                self.finalize_training(start_time, all_positions, positions,'RBSM',  gamma1, gamma2) 
                break
        start_time = time.time()
        filename = f'RBSMplacement.npy'  # 文件名
        np.save(filename, positions.cpu().detach().numpy())
        np.save('block.py', self.blocks)
        while True:
            optimizer_distance.zero_grad()
            optimizer_penalty.zero_grad()
            all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
            loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)
            loss.backward(retain_graph=True)
            optimizer_distance.step()
            #print(f"loss={loss}")
            positions = self.update_penalty(positions)
            #loss = Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
            #loss.backward(retain_graph=True)
            #optimizer_penalty.step()
            if step % 5 == 0 or step == self.max_epoch*25-1:
                current_lr1 = scheduler_distance.get_last_lr()[0]
                current_lr2 = scheduler_penalty.get_last_lr()[0]
                lr_logs1.append(current_lr1)
                lr_logs2.append(current_lr2)
                current_index = len(lr_logs1) - 1 
                #print(f"current_index={current_index}")
                all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
                full_loss = loss + Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
                self.loss_history.append(full_loss.item())
                #print(f"Epoch {step}:  Loss = {full_loss.item():.4f}\t distance_lr={lr_logs1[current_index]}\t penalty_lr={lr_logs2[current_index]}")
            overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1) 
            scheduler_distance.step()
            scheduler_penalty.step()
            if overlap == 0:
                break
        all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
        #print(f"gamma1, gamma2={gamma1, gamma2}")
        np.save('legalRBSMplacement.npy', positions.cpu().detach().numpy())
        self.finalize_training(start_time, all_positions, positions, 'legalRBSM', gamma1, gamma2)


    def train_GD(self, optimizer):
        start_time = time.time()
        k = 1
        #gamma1 = torch.randn(self.sb_count, requires_grad=True) 
        #gamma2 = torch.randn(self.sb_count*(self.sb_count-1)//2, requires_grad=True) 
        gamma1 = torch.full((self.sb_count,), 10000.0, requires_grad=True)
        gamma2 = torch.full((self.sb_count*(self.sb_count-1)//2,), 10000.0, requires_grad=True)
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.95)
        #scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=self.max_epoch)
        lr_logs = []
        positions = self.model()
        last_lengths = []
        convergence_threshold = 0.0001  # 0.01%
        self.overlap_threshold = self.area * 0.01
        s = 0
        print(f"overlap_threshold={self.overlap_threshold}")
        max_convergences = 3
        last_overlaps = []
        stop_training = False
        Loss = loss_function(gamma1, gamma2, self.W , self.H)
        while True:
            scheduler.step()
            #new_gamma1, new_gamma2 = self.calculate_gamma(positions) 
            #gamma1 = torch.max(gamma1, 100*new_gamma1)
            #gamma2 = torch.max(gamma2, 100*new_gamma2)  
            #Loss = loss_function(gamma1, gamma2, self.W , self.H)
            for step in range(25):
                optimizer.zero_grad()
                all_positions = torch.cat((positions,  self.fixed_positions), dim=0)        
                loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions) + Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)     
                loss.backward(retain_graph=True)
                overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1)    
                #grad = positions.grad.clone()
                #eps = 0.2/((k) ** 3)
                #grad_length = torch.norm(grad).item()
                #direction =  np.random.normal(0, 1, grad.shape)
                #direction = torch.tensor(direction, dtype=grad.dtype, device=grad.device)  # 转换为 PyTorch 张量
                #new_direction = grad + eps*grad_length*direction
                #positions.grad = new_direction
                optimizer.step()
                if step % 5 == 0: 
                    current_lr = scheduler.get_last_lr()[0]
                    lr_logs.append(current_lr)
                    self.loss_history.append(loss.item())
                    print(f"Epoch {(k-1)*25+step}:  Loss = {loss.item():.4f}\t overlap = {overlap} \t lr={current_lr}")
                if len(last_lengths) >= max_convergences and all(abs(last_lengths[i] - last_lengths[i+1]) < convergence_threshold * last_lengths[i] for i in range(len(last_lengths)-1)):
                    if overlap < self.overlap_threshold:  # Check overlap condition only in the last step
                        print("Training stopped due to minimal change in wire length and acceptable overlap.")
                        stop_training = True
                        break
                last_lengths.append(loss.item())
                last_overlaps.append(overlap)
                if len(last_lengths) > max_convergences:
                    last_lengths.pop(0)  # keep only the last few measurements
            k += 1
            scheduler.step()
            if stop_training or k == self.max_setp:
                self.finalize_training(start_time, all_positions, positions,'GD',  gamma1, gamma2) 
                break
        start_time = time.time()
        filename = f'GDplacement.npy'  # 文件名
        np.save(filename, positions.cpu().detach().numpy())
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.95)    
        while True:
            optimizer.zero_grad()
            all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
            loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)
            loss.backward(retain_graph=True)
            #optimizer.step()
            grad = positions.grad.clone()
            positions = positions -current_lr*grad
            positions = self.update_penalty(positions)
            all_positions = torch.cat((positions,  self.fixed_positions), dim=0) 
            if s % 5 == 0:
                current_lr = scheduler.get_last_lr()[0]
                lr_logs.append(current_lr)
                current_index = len(lr_logs) - 1 
                full_loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)+Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
                self.loss_history.append(full_loss.item())
               # print(f"Epoch {s}: Loss = {full_loss.item():10.3e}\t lr={lr_logs[current_index]}")
            scheduler.step()
            overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1) 
            s += 1
            if overlap == 0 or s == 200:
                break
        all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
        print(f"gamma1, gamma2={gamma1, gamma2}")
        self.finalize_training(start_time, all_positions, positions,'legalGD',  gamma1, gamma2)

    def train_Adam1(self, optimizer):
        start_time = time.time()
        epoch1 = 35
        p_positions = P_Position(self.blocks,self.sb_count,self.pfile)
        #gamma1 = torch.randn(self.sb_count, requires_grad=True)
        #gamma2 = torch.randn(self.sb_count * (self.sb_count - 1) // 2, requires_grad=True)
        gamma1 = torch.full((self.sb_count,), 10000.0, requires_grad=True)
        gamma2 = torch.full((self.sb_count*(self.sb_count-1)//2,), 10000.0, requires_grad=True)
        #scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.95)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=self.max_epoch)
        lr_logs = []
        k = 1
        last_lengths = []
        convergence_threshold = 0.0001  # 0.01%
        self.overlap_threshold = self.area * 0.01
        print(f"overlap_threshold={self.overlap_threshold}")
        max_convergences = 3
        last_overlaps = []
        stop_training = False 
        while True:
            positions = self.model()
            #new_gamma1, new_gamma2 = self.calculate_gamma(positions) 
            #gamma1 = torch.max(gamma1, 100*new_gamma1)
            #gamma2 = torch.max(gamma2, 100*new_gamma2)    
            #print(f"gamma1, gamma2={gamma1, gamma2}")
            #gamma1 =( 0.5*gamma1 + 0.5*100 *new_gamma1  )
           # gamma2 =( 0.5*gamma2 + 0.5*100 *new_gamma2  )  
            Loss = loss_function(gamma1,gamma2, self.W , self.H)
            for step in range(25):
                positions = self.model()
                all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
                optimizer.zero_grad()
                loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)+Loss.calculate_boundary_loss1(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)            
                loss.backward(retain_graph=True)
                overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1)   
                optimizer.step()
                    #print(f"loss={loss}")
                if step % 5 == 0 :
                        #with torch.no_grad():
                        current_lr = scheduler.get_last_lr()[0]
                        lr_logs.append(current_lr)
                        current_index = len(lr_logs) - 1 
                        self.loss_history.append(loss.item())
                       # print(f"Epoch {k*25+step+1}:  Loss = {loss.item():.4f}\t overlap = {overlap} \t lr={lr_logs[current_index]}")
                if len(last_lengths) >= max_convergences and all(abs(last_lengths[i] - last_lengths[i+1]) < convergence_threshold * last_lengths[i] for i in range(len(last_lengths)-1)):
                    if overlap < self.overlap_threshold:  # Check overlap condition only in the last step
                        print("Training stopped due to minimal change in wire length and acceptable overlap.")
                        stop_training = True
                        break
                last_lengths.append(loss.item())
                last_overlaps.append(overlap)
                if len(last_lengths) > max_convergences:
                    last_lengths.pop(0)  # keep only the last few measurements
            k += 1
            scheduler.step()
            if stop_training or k == self.max_setp:
                self.finalize_training(start_time, all_positions, positions,'Adam',  gamma1, gamma2) 
                break
        start_time = time.time()
        #self.finalize_training(start_time, all_positions, positions, 'middleAdam', gamma1, gamma2)
        while True:
            optimizer.zero_grad()
            all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
            loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)
            loss.backward(retain_graph=True)
            optimizer.step()
            positions = self.update_penalty(positions)
            if step % 5 == 0:
                current_lr = scheduler.get_last_lr()[0]
                lr_logs.append(current_lr)
                current_index = len(lr_logs) - 1 
                #print(f"current_index={current_index}")
                full_loss = loss + Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
                self.loss_history.append(full_loss.item())
                #print(f"Epoch {step}: Loss = {full_loss.item():10.3e}\t lr={lr_logs[current_index]}")
            scheduler.step()
            overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1) 
            if overlap == 0:
                break
        all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
       # print(f"gamma1, gamma2={gamma1, gamma2}")
        self.finalize_training(start_time, all_positions, positions,'legalAdam',  gamma1, gamma2)

    def train_Adam(self, optimizer_distance, optimizer_penalty):
        start_time = time.time()
        lr_logs1, lr_logs2 = [], []
        positions = self.model()
        gamma1 = torch.full((self.sb_count,), 1000.0, requires_grad=True)
        gamma2 = torch.full((self.sb_count*(self.sb_count-1)//2,), 1000.0, requires_grad=True)
        scheduler_distance = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer_distance, T_max=self.max_setp)
        scheduler_penalty = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer_penalty, T_max=self.max_setp)
        #scheduler_distance = torch.optim.lr_scheduler.ExponentialLR(optimizer_distance, gamma=0.95)
        #scheduler_penalty = torch.optim.lr_scheduler.ExponentialLR(optimizer_penalty, gamma=0.95)
        batch_size = math.ceil(self.num_nets / 25)
        probabilities = self.sample(self.nets_find)
        k = 1
        last_lengths = []
        convergence_threshold = 0.0001  # 0.01%
        self.overlap_threshold = self.area * 0.01
        print(f"overlap_threshold={self.overlap_threshold}")
        max_convergences = 3
        last_overlaps = []
        stop_training = False 
        while True:
            new_gamma1, new_gamma2 = self.calculate_gamma(positions)
            gamma1 = torch.max(gamma1, 100 * new_gamma1)
            gamma2 = torch.max(gamma2, 100 * new_gamma2)
            #gamma1 =( 0.5*gamma1 + 0.5*100 *new_gamma1  )
            #gamma2 =( 0.5*gamma2 + 0.5*100 *new_gamma2  )  
            Loss = loss_function(gamma1, gamma2, self.W, self.H)
            for step in range(25):
                optimizer_distance.zero_grad()
                if k < 20:
                    batch_idx = np.random.choice(self.num_nets, size=batch_size, replace=False, p=probabilities)
                else:
                    batch_idx = self.net_idx
                batch_idx = self.net_idx
                #print("batch_idx",batch_idx)
                all_positions = torch.cat((positions, self.fixed_positions), dim=0)
                distance_loss = Loss.calculate_distance_loss(self.nets, batch_idx, all_positions)
                distance_loss.backward(retain_graph=True)
                optimizer_distance.step()
                #self.update_gradient(positions, distance_loss, optimizer_distance, k)
                penalty_loss = Loss.calculate_boundary_loss1(self.full_idx, positions, self.block1) + \
                                Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
                penalty_loss.backward(retain_graph=True)
                optimizer_penalty.step()
                #self.update_gradient(positions, penalty_loss, optimizer_penalty, k)
                loss = distance_loss + penalty_loss
                overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1)    
                if step % 5 == 0: 
                    with torch.no_grad():
                        current_lr1 = scheduler_distance.get_last_lr()[0]
                        current_lr2 = scheduler_penalty.get_last_lr()[0]
                        lr_logs1.append(current_lr1)
                        lr_logs2.append(current_lr2)
                        self.loss_history.append(loss.item())
                        print(f"Epoch {(k-1)*25+step+1}:  Loss = {loss.item():.4f}\t diatance={distance_loss}\t overlap ={overlap}\t  distance_lr={current_lr1}\t penalty_lr={current_lr2 }")
                        #print(f"Epoch {k}:  Loss = {loss.item():.4f}\t distance_lr={current_lr1}\t penalty_lr={current_lr2 }")
                # Check convergence based on the last 5 total losses
                if len(last_lengths) >= max_convergences and all(abs(last_lengths[i] - last_lengths[i+1]) < convergence_threshold * last_lengths[i] for i in range(len(last_lengths)-1)):
                    if overlap < self.overlap_threshold:  # Check overlap condition only in the last step
                        print("Training stopped due to minimal change in wire length and acceptable overlap.")
                        stop_training = True
                        break
                last_lengths.append(distance_loss.item())
                last_overlaps.append(overlap)
                if len(last_lengths) > max_convergences:
                    last_lengths.pop(0)  # keep only the last few measurements
            k += 1
            scheduler_distance.step()
            scheduler_penalty.step()
            if stop_training or k == self.max_setp:
                self.finalize_training(start_time, all_positions, positions,'Adam',  gamma1, gamma2) 
                break
        start_time = time.time()
        filename = f'Adamplacement.npy'  # 文件名
        np.save(filename, positions.cpu().detach().numpy())
        while True:
            optimizer_distance.zero_grad()
            optimizer_penalty.zero_grad()
            all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
            loss = Loss.calculate_distance_loss(self.nets, self.net_idx, all_positions)
            loss.backward(retain_graph=True)
            optimizer_distance.step()
            #print(f"loss={loss}")
            positions = self.update_penalty(positions)
            #loss = Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
            #loss.backward(retain_graph=True)
            #optimizer_penalty.step()
            if step % 5 == 0 or step == self.max_epoch*25-1:
                current_lr1 = scheduler_distance.get_last_lr()[0]
                current_lr2 = scheduler_penalty.get_last_lr()[0]
                lr_logs1.append(current_lr1)
                lr_logs2.append(current_lr2)
                current_index = len(lr_logs1) - 1 
                #print(f"current_index={current_index}")
                all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
                full_loss = loss + Loss.calculate_boundary_loss(self.full_idx, positions, self.block1) + Loss.calculate_overlap_loss(self.full_idx, positions, self.block1)
                self.loss_history.append(full_loss.item())
                #print(f"Epoch {step}:  Loss = {full_loss.item():.4f}\t distance_lr={lr_logs1[current_index]}\t penalty_lr={lr_logs2[current_index]}")
            overlap =  Loss.calculate_overlap_loss0(self.full_idx, positions, self.block1) 
            scheduler_distance.step()
            scheduler_penalty.step()
            if overlap == 0:
                break
        all_positions = torch.cat((positions,  self.fixed_positions), dim=0)
        #print(f"gamma1, gamma2={gamma1, gamma2}")
        self.finalize_training(start_time, all_positions, positions, 'legalAdam', gamma1, gamma2)    

class TrainingPipeline:
    def __init__(self, block_file, net_file, pl_file, W, H, max_epoch,  batch_size, max_step):
        self.blocks, self.block_order, self.sb_count, self.block_dict, self.area = Blocks(block_file).read_block()
        self.pl_file = pl_file
        self.block_file = block_file
        self.net_file = net_file
        self.W = W
        self.H = H
        self.max_epoch = max_epoch
        self.batch_size = batch_size
        self.model = Position(self.sb_count,self.W,self.H)
        self.train_model_instance = TrainModel( self.model,self.net_file, self.block_file, self.pl_file, self.W, self.H,  max_epoch, batch_size, max_step)

    def train_rbsm(self, lr_distance, lr_penalty):
        optimizer_distance = optim.SGD([self.model.positions], lr=lr_distance)
        optimizer_penalty = optim.SGD([self.model.positions], lr=lr_penalty)
        self.train_model_instance.train_RBSM(optimizer_distance, optimizer_penalty)

    def train_sgd(self, lr):
        #optimizer_sgd = optim.SGD(self.model.parameters(), lr=lr, weight_decay=0.0, momentum=0.8)
        optimizer_sgd = optim.SGD(self.model.parameters(), lr=lr)
        self.train_model_instance.train_GD(optimizer_sgd)
        
    def train_adam2(self):
        optimizer_adam = optim.Adam(
            self.model.parameters(),
            lr=0.2, 
            #max_steps=(0.85, 0.995),
            eps=1e-08, 
            weight_decay=0.00, 
            amsgrad=True
        )
        self.train_model_instance.train_Adam2(optimizer_adam)

    def train_adam(self):
        optimizer_adam1 = optim.Adam(
            self.model.parameters(),
            lr=0.05, 
            #max_steps=(0.85, 0.995),
            eps=1e-08, 
            weight_decay=0.00, 
            amsgrad=True
        )
        optimizer_adam2 = optim.Adam(
            self.model.parameters(),   
            lr=0.1, 
            #max_steps=(0.85, 0.995),
            eps=1e-08, 
            weight_decay=0.00, 
            amsgrad=True
        )
        self.train_model_instance.train_Adam(optimizer_adam1, optimizer_adam2)
                                             
def main():
    # 定义参数
    W, H = 800, 800  # Example values
    max_epoch = 40
    
    #pipeline = TrainingPipeline(GCSR_DIR / 'n10/n10_blocks.txt', GCSR_DIR / 'n10/n10_nets.txt', GCSR_DIR / 'n10/n10.pl', W, H,  max_epoch, batch_size=2, max_step=100)
    #pipeline = TrainingPipeline(GCSR_DIR / 'n30/n30_blocks.txt', GCSR_DIR / 'n30/n30_nets.txt', GCSR_DIR / 'n30/n30.pl', W, H, max_epoch,  batch_size=6, max_step=100)
    #pipeline = TrainingPipeline(GCSR_DIR / 'n50/n50_blocks.txt', GCSR_DIR / 'n50/n50_nets.txt', GCSR_DIR / 'n50/n50.pl', W, H,  max_epoch,  batch_size=10, max_step=60)
    #pipeline = TrainingPipeline(GCSR_DIR / 'n100/n100_blocks.txt', GCSR_DIR / 'n100/n100_nets.txt', GCSR_DIR / 'n100/n100.pl', W, H, max_epoch, batch_size=20,max_step=65)
    pipeline = TrainingPipeline(GCSR_DIR / 'n200/n200_blocks.txt', GCSR_DIR / 'n200/n200_nets.txt', GCSR_DIR / 'n200/n200.pl', W, H, max_epoch, batch_size=40,max_step=65)
    #pipeline = TrainingPipeline(GCSR_DIR / 'n300/n300_blocks.txt', GCSR_DIR / 'n300/n300_nets.txt', GCSR_DIR / 'n300/n300.pl', W, H,  max_epoch, batch_size=60, max_step=65)
    
    # 选择训练方法
    #method = 'GD'
    method = 'RBSM' 
    #method = 'Adam'  
    #method = 'Adam2'    

    if method == 'RBSM':
        pipeline.train_rbsm(lr_distance=0.1, lr_penalty=0.1)
    elif method == 'GD':    
        pipeline.train_sgd(lr=0.1)
    #elif method == 'Adam':
    #pipeline.train_adam()
    elif method == 'Adam':
       pipeline.train_adam()
    #else:  
    #    print("Unknown training method")
 

if __name__ == "__main__":
    main()
                                                                                                                                                                         
