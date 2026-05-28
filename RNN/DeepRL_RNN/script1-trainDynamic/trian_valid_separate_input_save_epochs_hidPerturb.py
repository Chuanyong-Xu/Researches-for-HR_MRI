import os
import random
from datetime import datetime

import numpy as np
import torch
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.optim as optim
import scipy.io
from scipy.io import savemat
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from tqdm import tqdm
from collections import deque
from dalaloader import extract_behavioral_data, RNNInputDataset, split_data_by_ratio, reshape_own
from model import RNN
from utils import set_seed, reshape_and_mask_tdev


def test_epoch_like_original(model,
                             dataset_test,
                             batch_size=1,
                             device='cpu'):
    """
    与原 test.py 里的 test() 保持一致，用当前 epoch 的模型参数提取：
    1) outputs / rnn_pre
    2) hidden_states / hidden_state
    3) R²

    注意：这里刻意不加入训练时的 reshape_and_mask_tdev，
    因为你原来的 test.py 在最终输出 update_data_switch_weight_diff0.mat 时也没有做这一步。
    """

    test_loader = DataLoader(dataset_test, batch_size=batch_size)
    model.to(device)
    model.eval()

    test_loss = 0.0
    criterion = nn.MSELoss()

    all_outputs = []
    all_targets = []
    all_hidden_states = []
    r2 = np.nan

    with torch.no_grad():
        for diff_tensor, switch_tensor, is_start, targets in tqdm(test_loader, desc="Testing"):
            diff_tensor = diff_tensor.float().to(device)
            switch_tensor = switch_tensor.float().to(device)
            is_start = is_start.float().to(device)
            targets = targets.float().to(device)

            outputs, hidden_states, _, is_switch1 = model(diff_tensor, switch_tensor, is_start)

            all_outputs.append(outputs.cpu())
            all_targets.append(targets.cpu())
            all_hidden_states.append(hidden_states.cpu())

            targets = targets.transpose(1, 2)

            has_inf = torch.isinf(targets[0])
            has_nan = torch.isnan(targets[0])
            if has_inf.any():
                valid_vals = targets[0][~has_inf & ~has_nan]
                max_val = valid_vals.max() if len(valid_vals) > 0 else torch.tensor(0.0)
                targets[0][has_inf] = max_val
                print("已将 inf 替换为最大值：", max_val.item())
            if has_nan.any():
                print("包含 nan，位置：", torch.nonzero(has_nan, as_tuple=True)[0])

            loss = criterion(outputs, targets)

            # 严格保持你原 test.py 的写法：r2_score(outputs.squeeze(), targets.squeeze())
            # 虽然 sklearn 标准用法通常是 r2_score(y_true, y_pred)，这里不改你的原逻辑。
            r2 = r2_score(outputs.squeeze().cpu(), targets.squeeze().cpu())
            test_loss += loss.item()

    avg_test_loss = test_loss / len(test_loader)

    print(f"Test Loss = {avg_test_loss:.4f}")
    print(f"Test R² = {r2:.4f}")

    all_outputs = torch.cat(all_outputs, dim=0)
    all_targets = torch.cat(all_targets, dim=0)
    all_hidden_states = torch.cat(all_hidden_states, dim=0)

    return avg_test_loss, all_outputs, all_hidden_states, r2, all_targets


def _iter_mat_rows(results):
    if getattr(results, 'ndim', None) == 0:
        return [results.item() if hasattr(results, 'item') else results]
    if isinstance(results, np.ndarray):
        return results.tolist()
    return [results]


def get_cond_id_from_mat(matlab_path, var_name='infer_data'):
    # condition order matches MATLAB: D1E0,D2E0,D3E0,D1E1,...,D3E2
    data = scipy.io.loadmat(matlab_path, struct_as_record=False, squeeze_me=True)
    rows = _iter_mat_rows(data[var_name])
    cond_all = []
    for row in rows:
        tDev = np.asarray(getattr(row, 'tDev')).reshape(-1)
        Nback = np.asarray(getattr(row, 'Nback')).reshape(-1)
        cond = np.zeros(len(tDev), dtype=np.int64)
        cond[(tDev == 1) & (Nback == 0)] = 1
        cond[(tDev == 2) & (Nback == 0)] = 2
        cond[(tDev == 3) & (Nback == 0)] = 3
        cond[(tDev == 1) & (Nback == 1)] = 4
        cond[(tDev == 2) & (Nback == 1)] = 5
        cond[(tDev == 3) & (Nback == 1)] = 6
        cond[(tDev == 1) & (Nback >= 2)] = 7
        cond[(tDev == 2) & (Nback >= 2)] = 8
        cond[(tDev == 3) & (Nback >= 2)] = 9
        cond_all.append(cond)
    return np.concatenate(cond_all, axis=0)


def read_target_from_mat(matlab_path, target_col='pr_of_switch', var_name='infer_data'):
    data = scipy.io.loadmat(matlab_path, struct_as_record=False, squeeze_me=True)
    rows = _iter_mat_rows(data[var_name])
    y_all = []
    for row in rows:
        y_all.append(np.asarray(getattr(row, target_col)).reshape(-1))
    return np.concatenate(y_all, axis=0).astype(np.float32)


def remove_condition_mean_geometry(hidden_np, cond_id):
    # Keep the original hidden trajectory/residuals, but remove 9-condition mean geometry.
    H = np.asarray(hidden_np, dtype=np.float32)
    if H.ndim == 3:
        H = H.reshape(-1, H.shape[-1])
    if H.ndim != 2:
        raise ValueError(f'hidden state should be 2D after reshape, got shape {H.shape}')
    cond_id = np.asarray(cond_id).reshape(-1)
    if len(cond_id) != H.shape[0]:
        raise ValueError(f'cond_id length {len(cond_id)} != hidden trials {H.shape[0]}')

    valid = cond_id > 0
    H_pert = H.copy()
    grand_mean = H[valid].mean(axis=0, keepdims=True)
    for cc in range(1, 10):
        idx = cond_id == cc
        if np.any(idx):
            cond_mean = H[idx].mean(axis=0, keepdims=True)
            H_pert[idx] = H[idx] - cond_mean + grand_mean
    return H_pert


def recompute_prediction_from_hidden(model, hidden_np, target_np, device='cpu'):
    H = np.asarray(hidden_np, dtype=np.float32)
    if H.ndim == 2:
        H_t = torch.tensor(H, dtype=torch.float32, device=device).unsqueeze(0)
    else:
        H_t = torch.tensor(H, dtype=torch.float32, device=device)
    model.to(device)
    model.eval()
    with torch.no_grad():
        pred = model.fc(H_t)
        if getattr(model, 'is_sigmoid', 1) == 1:
            pred = torch.sigmoid(pred)
    pred_np = pred.squeeze().detach().cpu().numpy().astype(np.float32)
    y = np.asarray(target_np, dtype=np.float32).reshape(-1)
    mse = float(np.mean((pred_np.reshape(-1) - y) ** 2))
    # standard R2: y_true first, y_pred second
    r2_standard = float(r2_score(y, pred_np.reshape(-1)))
    # keep an original-style R2 for comparability with your previous title values
    r2_original_style = float(r2_score(pred_np.reshape(-1), y))
    return pred_np, mse, r2_standard, r2_original_style


def pre_save_in_mat(matlab_path,
                    pre,
                    col_A="tDev",
                    new_col="rnn_pre",
                    is_train_col="is_train",
                    hidden_col="hidden_state",
                    hidden_states=None,
                    var_name='infer_data',
                    save_path="../data/update_data.mat",
                    epoch=None,
                    r2_epoch=None,
                    train_loss_epoch=None,
                    val_loss_epoch=None,
                    test_loss_epoch=None,
                    pre_geo_rm=None,
                    new_col_geo_rm='rnn_pre_geo_rm',
                    r2_geo_rm_epoch=None,
                    r2_geo_rm_epoch_original_style=None,
                    test_loss_geo_rm_epoch=None,
                    delta_loss_geo_rm_epoch=None):
    """
    基本复制你原 test.py 的 pre_save_in_mat()。
    唯一增加：在 epoch 文件最外层保存 epoch / r2_epoch / loss，方便 MATLAB 直接读取标题信息。
    """
    import scipy.io
    from scipy.io import savemat

    data = scipy.io.loadmat(matlab_path, struct_as_record=False, squeeze_me=True)

    if var_name not in data:
        raise ValueError(f"变量 '{var_name}' 不在 .mat 文件中，实际变量有：{list(data.keys())}")

    results = data[var_name]
    if results.ndim == 0:
        results = [results]
    elif isinstance(results, np.ndarray):
        results = results.tolist()

    total_trials = sum(getattr(row, col_A).shape[0] for row in results)

    n60 = int(total_trials * 0.6)
    n20 = int(total_trials * 0.2)
    is_train_array = np.concatenate([
        np.ones(n60),
        np.zeros(n20),
        np.full(total_trials - n60 - n20, 2)
    ]).astype(int)

    arr_idx = 0
    for row in results:
        trials = getattr(row, col_A).shape[0]

        pred_chunk = pre[arr_idx: arr_idx + trials]
        pred_chunk = pred_chunk.reshape((trials, -1))
        setattr(row, new_col, pred_chunk)

        if pre_geo_rm is not None:
            pred_geo_chunk = pre_geo_rm[arr_idx: arr_idx + trials]
            pred_geo_chunk = pred_geo_chunk.reshape((trials, -1))
            setattr(row, new_col_geo_rm, pred_geo_chunk)

        is_train_chunk = is_train_array[arr_idx: arr_idx + trials].reshape((trials, 1))
        setattr(row, is_train_col, is_train_chunk)

        if hidden_states is not None:
            hidden_chunk = hidden_states[arr_idx: arr_idx + trials]
            setattr(row, hidden_col, hidden_chunk)

        arr_idx += trials

    if arr_idx != len(pre):
        print(f"警告：还有未分配的数据，共剩余 {len(pre) - arr_idx} 个元素")

    save_dict = {var_name: np.array(results)}
    if epoch is not None:
        save_dict['epoch'] = np.array([[epoch]])
    if r2_epoch is not None:
        save_dict['r2_epoch'] = np.array([[r2_epoch]])
    if train_loss_epoch is not None:
        save_dict['train_loss_epoch'] = np.array([[train_loss_epoch]])
    if val_loss_epoch is not None:
        save_dict['val_loss_epoch'] = np.array([[val_loss_epoch]])
    if test_loss_epoch is not None:
        save_dict['test_loss_epoch'] = np.array([[test_loss_epoch]])
    if r2_geo_rm_epoch is not None:
        save_dict['r2_geo_rm_epoch'] = np.array([[r2_geo_rm_epoch]])
    if r2_geo_rm_epoch_original_style is not None:
        save_dict['r2_geo_rm_epoch_original_style'] = np.array([[r2_geo_rm_epoch_original_style]])
    if test_loss_geo_rm_epoch is not None:
        save_dict['test_loss_geo_rm_epoch'] = np.array([[test_loss_geo_rm_epoch]])
    if delta_loss_geo_rm_epoch is not None:
        save_dict['delta_loss_geo_rm_epoch'] = np.array([[delta_loss_geo_rm_epoch]])

    savemat(save_path, save_dict)


def save_epoch_mat(model,
                   dataset_save,
                   batch_size,
                   device,
                   matlab_path,
                   save_mat_dir,
                   save_mat_prefix,
                   epoch,
                   train_loss_epoch=np.nan,
                   val_loss_epoch=np.nan):
    test_loss_epoch, epoch_outputs, epoch_hidden_state, r2_epoch, epoch_targets = test_epoch_like_original(
        model=model,
        dataset_test=dataset_save,
        batch_size=batch_size,
        device=device)

    epoch_hidden_state = epoch_hidden_state.squeeze().cpu().detach().numpy()
    epoch_outputs = epoch_outputs.squeeze().cpu().detach().numpy()

    # Perturb hidden states only before readout: remove 9-condition mean geometry,
    # then use the same trained fc layer to recompute prediction.
    cond_id = get_cond_id_from_mat(matlab_path)
    target_np = read_target_from_mat(matlab_path)
    epoch_hidden_geo_rm = remove_condition_mean_geometry(epoch_hidden_state, cond_id)
    epoch_outputs_geo_rm, test_loss_geo_rm_epoch, r2_geo_rm_epoch, r2_geo_rm_epoch_original_style = recompute_prediction_from_hidden(
        model=model,
        hidden_np=epoch_hidden_geo_rm,
        target_np=target_np,
        device=device)
    delta_loss_geo_rm_epoch = test_loss_geo_rm_epoch - test_loss_epoch

    save_path = os.path.join(save_mat_dir, f"{save_mat_prefix}_epoch{epoch}.mat")
    pre_save_in_mat(matlab_path,
                    epoch_outputs,
                    hidden_states=epoch_hidden_state,
                    save_path=save_path,
                    epoch=epoch,
                    r2_epoch=r2_epoch,
                    train_loss_epoch=train_loss_epoch,
                    val_loss_epoch=val_loss_epoch,
                    test_loss_epoch=test_loss_epoch,
                    pre_geo_rm=epoch_outputs_geo_rm,
                    r2_geo_rm_epoch=r2_geo_rm_epoch,
                    r2_geo_rm_epoch_original_style=r2_geo_rm_epoch_original_style,
                    test_loss_geo_rm_epoch=test_loss_geo_rm_epoch,
                    delta_loss_geo_rm_epoch=delta_loss_geo_rm_epoch)
    print(f"Epoch {epoch} data saved to {save_path}, R² = {r2_epoch:.4f}, geo-removed loss = {test_loss_geo_rm_epoch:.4f}")


def train_val_loop(model,
                   dataset_train,
                   dataset_valid,
                   save_model_path,
                   epochs=1000,
                   batch_size=1,
                   lr=1e-3,
                   device='cpu',
                   shuffle=False,
                   early_stop_patience=5,
                   early_stop_delta=1e-4,
                   model_name="best_model.pt",
                   save_epoch_data=False,
                   save_every=20,
                   dataset_save=None,
                   matlab_path="../data/infer_data.mat",
                   save_mat_dir="../data",
                   save_mat_prefix="update_data_switch_weight_diff", ):
    """
    训练验证函数

    :param model_name:
    :param save_model_path:
    :param model:
    :param dataset_train:
    :param dataset_valid:
    :param epochs:
    :param batch_size:
    :param lr:
    :param device:
    :param shuffle:
    :param early_stop_patience:
    :param early_stop_delta:
    :return:
    """

    os.makedirs(save_model_path, exist_ok=True)
    os.makedirs(save_mat_dir, exist_ok=True)
    if dataset_save is None:
        dataset_save = dataset_train

    # 加载器
    train_loader = DataLoader(dataset_train, batch_size=batch_size, shuffle=shuffle)
    val_loader = DataLoader(dataset_valid, batch_size=batch_size)

    # 模型、损失、优化器
    model.to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    # Early stopping
    best_val_loss = float('inf')
    no_improve_epochs = 0

    # 保存 epoch 0, 1-10, 20,30,40...直到设定的最后 epoch
    save_epochs = set(list(range(1, 11)) + list(range(20, epochs + 1, save_every)) + [epochs])

    # epoch 0：未训练前的随机初始化模型
    if save_epoch_data:
        save_epoch_mat(model=model,
                       dataset_save=dataset_save,
                       batch_size=batch_size,
                       device=device,
                       matlab_path=matlab_path,
                       save_mat_dir=save_mat_dir,
                       save_mat_prefix=save_mat_prefix,
                       epoch=0,
                       train_loss_epoch=np.nan,
                       val_loss_epoch=np.nan)

    for epoch in range(1, epochs + 1):
        model.train()
        train_loss = 0.0

        for diff_tensor, switch_tensor, is_start_tensor, targets in tqdm(train_loader,
                                                                         desc=f"[Epoch {epoch}] Training"):
            diff_tensor, switch_tensor, is_start_tensor, targets = diff_tensor.float().to(
                device), switch_tensor.float().to(
                device), is_start_tensor.float().to(device), targets.float().to(device)
            # print(diff_tensor.shape)
            # print(switch_tensor.shape)
            # asdasd

            # 当switch在当前trial为1时，让diff变为[0, 0, 0]
            diff_tensor[0], switch_tensor[0] = reshape_and_mask_tdev(diff=diff_tensor[0],
                                                                     switch=switch_tensor[0])
            # for n in range(diff_tensor.shape[1]):
            #     print(diff_tensor[0][n])
            #     print(switch_tensor[0][n])
            #     print(targets[0][0][n])
            #     print("--------------------")
            # dfgvdfv

            # print(diff_tensor.shape)
            # print(switch_tensor.shape)
            # print(diff_tensor[0])
            # print(switch_tensor[0])
            # asdasd

            optimizer.zero_grad()
            pre, out, hn, is_switch1 = model(diff_tensor, switch_tensor, is_start_tensor)  # 拆出 difficulty 和 switch

            # # 构造 ground truth: 当前 step 是否为 [0, 1]
            # switch_is_01 = (switch_tensor[:, :, 1] == 1).float().unsqueeze(-1)  # shape: (B, T, 1)
            # # loss_s1：模型对switch==[0,1]的预测准确性（二分类）
            # loss_s1 = nn.BCEWithLogitsLoss()(is_switch1, switch_is_01)

            batch_loss = 0.0
            loss_count = 0
            mistake_count = 0  # 记录连续错的次数

            for t in range(pre.size(1)):
                pre_t = pre[:, t, :]
                target_t = targets[:, :, t]

                current_switch = switch_tensor[0][t]  # shape: (2,)
                if torch.all(current_switch == torch.tensor([0.0, 1.0], device=device)):
                    weight = 1.0
                    mistake_count = 0  # 正确，重置
                elif torch.all(current_switch == torch.tensor([1.0, 0.0], device=device)):
                    mistake_count += 1
                    if mistake_count == 1:
                        weight = 3.0
                    elif mistake_count == 2:
                        weight = 6.0
                    else:
                        weight = 10.0
                else:
                    weight = 1.0  # 如果不是严格 [0,1] 或 [1,0]，默认不加权
                    mistake_count = 0  # 也可考虑不重置，取决于你定义是否严格判断

                loss_t = weight * nn.functional.mse_loss(pre_t, target_t)
                batch_loss += loss_t
                loss_count += 1

            batch_loss = (batch_loss / loss_count)

            # targets = targets.transpose(1, 2)
            #
            # has_inf = torch.isinf(targets[0])
            # has_nan = torch.isnan(targets[0])
            # # 打印结果
            # if has_inf.any():
            #     valid_vals = targets[0][~has_inf & ~has_nan]
            #     max_val = valid_vals.max() if len(valid_vals) > 0 else torch.tensor(0.0)
            #     targets[0][has_inf] = max_val
            #     print("已将 inf 替换为最大值：", max_val.item())
            # if has_nan.any():
            #     print("包含 nan，位置：", torch.nonzero(has_nan, as_tuple=True)[0])
            # #
            # # print(pre[0][0])
            # # print(targets[0][0])
            # loss = criterion(pre, targets)
            # # print(loss)
            # # asdasd
            batch_loss.backward()
            optimizer.step()
            train_loss += batch_loss.item()

        avg_train_loss = train_loss / len(train_loader)

        # 验证
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for diff_tensor, switch_tensor, is_start, targets in val_loader:
                diff_tensor, switch_tensor, is_start, targets = diff_tensor.float().to(
                    device), switch_tensor.float().to(device), is_start.float().to(device), targets.float().to(device)
                outputs, _, _, is_switch2 = model(diff_tensor, switch_tensor, is_start)

                # # 构造 ground truth: 当前 step 是否为 [0, 1]
                # switch_is_01 = (switch_tensor[:, :, 1] == 1).float().unsqueeze(-1)  # shape: (B, T, 1)
                # # loss_s1：模型对switch==[0,1]的预测准确性（二分类）
                # loss_s1 = nn.BCEWithLogitsLoss()(is_switch2, switch_is_01)

                batch_loss = 0.0
                loss_count = 0
                mistake_count = 0

                for t in range(outputs.size(1)):
                    pre_t = outputs[:, t, :]
                    target_t = targets[:, :, t]

                    current_switch = switch_tensor[0][t]  # shape: (2,)
                    if torch.all(current_switch == torch.tensor([0.0, 1.0], device=device)):
                        weight = 1.0
                        mistake_count = 0  # 正确，重置
                    elif torch.all(current_switch == torch.tensor([1.0, 0.0], device=device)):
                        mistake_count += 1
                        if mistake_count == 1:
                            weight = 3.0
                        elif mistake_count == 2:
                            weight = 6.0
                        else:
                            weight = 10.0
                    else:
                        weight = 1.0  # 如果不是严格 [0,1] 或 [1,0]，默认不加权

                    loss_t = weight * nn.functional.mse_loss(pre_t, target_t)

                    batch_loss += loss_t
                    loss_count += 1

                val_loss += (batch_loss / loss_count).item()

                # targets = targets.transpose(1, 2)
                #
                # has_inf = torch.isinf(targets[0])
                # has_nan = torch.isnan(targets[0])
                # # 打印结果
                # if has_inf.any():
                #     valid_vals = targets[0][~has_inf & ~has_nan]
                #     max_val = valid_vals.max() if len(valid_vals) > 0 else torch.tensor(0.0)
                #     targets[0][has_inf] = max_val
                #     print("已将 inf 替换为最大值：", max_val.item())
                # if has_nan.any():
                #     print("包含 nan，位置：", torch.nonzero(has_nan, as_tuple=True)[0])

                # batch_loss = criterion(outputs, targets)
                # val_loss += batch_loss.item()

        avg_val_loss = val_loss / len(val_loader)

        print(f"Epoch {epoch}: Train Loss = {avg_train_loss:.4f}, Val Loss = {avg_val_loss:.4f}")

        # 提取并保存指定 epoch 的 rnn_pre / hidden_state / R² / loss
        if save_epoch_data and epoch in save_epochs:
            save_epoch_mat(model=model,
                           dataset_save=dataset_save,
                           batch_size=batch_size,
                           device=device,
                           matlab_path=matlab_path,
                           save_mat_dir=save_mat_dir,
                           save_mat_prefix=save_mat_prefix,
                           epoch=epoch,
                           train_loss_epoch=avg_train_loss,
                           val_loss_epoch=avg_val_loss)

        # Early stopping 判断
        if avg_val_loss + early_stop_delta < best_val_loss:
            best_val_loss = avg_val_loss
            no_improve_epochs = 0

            # 保存当前最优模型参数
            best_model_path = os.path.join(save_model_path, model_name)
            torch.save(model.state_dict(), best_model_path)
            print(f"Best model saved to {best_model_path}")

        else:
            no_improve_epochs += 1
            if no_improve_epochs >= early_stop_patience:
                # 如果 early stopping 的最后一个 epoch 不在保存列表里，也保存一次，避免最后状态丢失
                if save_epoch_data and epoch not in save_epochs:
                    save_epoch_mat(model=model,
                                   dataset_save=dataset_save,
                                   batch_size=batch_size,
                                   device=device,
                                   matlab_path=matlab_path,
                                   save_mat_dir=save_mat_dir,
                                   save_mat_prefix=save_mat_prefix,
                                   epoch=epoch,
                                   train_loss_epoch=avg_train_loss,
                                   val_loss_epoch=avg_val_loss)
                print(f"Early stopping at epoch {epoch} (no improvement in {early_stop_patience} epochs)")
                break


if __name__ == "__main__":
    set_seed(42)
    model = RNN(is_sigmoid=1,
                hidden_size=32, )

    """
     预测切换概率
    """
    tdev, tf, target, is_start = extract_behavioral_data('../data/infer_data.mat',
                                                         var_name='infer_data')
    """
    预测confidence
    """
    # tdev, tf, target, is_start = extract_behavioral_data('../data/infer_data.mat',
    #                                                      row1='tDev',
    #                                                      row2='TF',
    #                                                      row3='mu_switch_estimated',
    #                                                      var_name='infer_data', )

    """
    无valid和test
    """
    tdev_all, tf_all, target_all, is_start_all = reshape_own(tdev, tf, target, is_start)

    # # 加强数据
    # tdev_all = np.tile(tdev_all, (1, 3))
    # tf_all = np.tile(tf_all, (1, 3))
    # target_all = np.tile(target_all, (1, 3))
    # is_start_all = np.tile(is_start_all, (1, 3))

    dataset_all = RNNInputDataset(tdev_all, tf_all, target_all, is_start_all)

    """
    有valid和test
    """
    # (tdev_train, tf_train, pr_train, is_start_train), (tdev_valid, tf_valid, pr_valid, is_start_valid), (
    # tdev_test, tf_test, pr_test, is_start_test) = split_data_by_ratio(tdev=tdev,
    #                                                                   tf=tf,
    #                                                                   pr_of_switch=pr_of_switch,
    #                                                                   is_start=is_start,
    #                                                                   train_ratio=0.6,
    #                                                                   test_ratio=0.2,
    #                                                                   valid_ratio=0.2, )
    # dataset_train = RNNInputDataset(tdev_train, tf_train, pr_train, is_start_train)
    # dataset_valid = RNNInputDataset(tdev_valid, tf_valid, pr_valid, is_start_valid)

    train_val_loop(model=model,
                   dataset_train=dataset_all,
                   dataset_valid=dataset_all,
                   save_model_path="../model",
                   early_stop_patience=10,
                   model_name="best_model_switch_weight_diff0.pt",
                   epochs=10000,
                   save_epoch_data=True,
                   save_every=20,
                   dataset_save=dataset_all,
                   matlab_path="../data/infer_data.mat",
                   save_mat_dir="../data",
                   save_mat_prefix="update_data_switch_weight_diff", )
