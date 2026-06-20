import os
import random
from datetime import datetime
from sklearn.metrics import r2_score, mean_squared_error

import numpy as np
import torch
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from collections import deque
from dalaloader import extract_behavioral_data_subjectwise, RNNInputDatasetSubjectwise
from model import RNN
from utils import set_seed, reshape_and_mask_tdev

def pre_save_in_mat(matlab_path,
                    pre,
                    col_A="tDev",
                    new_col="rnn_pre",
                    is_train_col="is_train",
                    hidden_col="hidden_state",
                    hidden_states=None,
                    var_name='infer_data',
                    save_path="../data/update_data.mat",
                    is_train_array=None):   # <<< ADD THIS
    """
    保存预测结果和隐藏状态到 .mat 文件中
    """
    import scipy.io
    import numpy as np
    from scipy.io import savemat

    # 加载数据
    data = scipy.io.loadmat(matlab_path, struct_as_record=False, squeeze_me=True)

    if var_name not in data:
        raise ValueError(f"变量 '{var_name}' 不在 .mat 文件中，实际变量有：{list(data.keys())}")

    results = data[var_name]
    if results.ndim == 0:
        results = [results]
    elif isinstance(results, np.ndarray):
        results = results.tolist()

    total_trials = sum(getattr(row, col_A).shape[0] for row in results)

    # 创建/使用 is_train 标签
    if is_train_array is None:
        n60 = int(total_trials * 0.6)
        n20 = int(total_trials * 0.2)
        is_train_array = np.concatenate([
            np.ones(n60),
            np.zeros(n20),
            np.full(total_trials - n60 - n20, 2)
        ]).astype(int)
    else:
        is_train_array = np.asarray(is_train_array).astype(int).reshape(-1)
        if len(is_train_array) != total_trials:
            raise ValueError(f"is_train_array length ({len(is_train_array)}) != total_trials ({total_trials})")

        # 遍历结构体
        arr_idx = 0
        for row in results:
            trials = getattr(row, col_A).shape[0]

            # 处理预测结果
            pred_chunk = pre[arr_idx: arr_idx + trials]
            pred_chunk = pred_chunk.reshape((trials, -1))  # 保证二维
            setattr(row, new_col, pred_chunk)

            # is_train
            is_train_chunk = is_train_array[arr_idx: arr_idx + trials].reshape((trials, 1))
            setattr(row, is_train_col, is_train_chunk)

            # hidden_state
            if hidden_states is not None:
                hidden_chunk = hidden_states[arr_idx: arr_idx + trials]
                setattr(row, hidden_col, hidden_chunk)

            arr_idx += trials

        if arr_idx != len(pre):
            print(f"警告：还有未分配的数据，共剩余 {len(pre) - arr_idx} 个元素")

        # 保存到新 .mat 文件
        savemat(save_path, {var_name: np.array(results)})



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
                   model_name="best_model.pt", ):
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

    for epoch in range(1, epochs + 1):
        model.train()
        train_loss = 0.0

        for outcome_tensor, switch_tensor, is_start_tensor, targets in tqdm(train_loader,
                                                                         desc=f"[Epoch {epoch}] Training"):
            ###### diff_tensor, switch_tensor, is_start_tensor, targets in tqdm(train_loader,
            ######                                                             desc=f"[Epoch {epoch}] Training"):
            outcome_tensor, switch_tensor, is_start_tensor, targets = outcome_tensor.float().to(
                device), switch_tensor.float().to(device), is_start_tensor.float().to(device), targets.float().to(device)
            # print(diff_tensor.shape)
            # print(switch_tensor.shape)
            # asdasd
            
            #print("outcome:", outcome_tensor.shape)
            #print("switch :", switch_tensor.shape)
            #print("start  :", is_start_tensor.shape) 
            
            # 当switch在当前trial为1时，让diff变为[0, 0, 0]
            outcome_tensor = reshape_and_mask_tdev(outcome_tensor,switch_tensor)
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
            pre, out, hn, is_switch1 = model(outcome_tensor, is_start_tensor)  # 拆出 difficulty 和 switch

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
            for outcome_tensor, switch_tensor, is_start, targets in val_loader:
                outcome_tensor, switch_tensor, is_start, targets = outcome_tensor.float().to(
                    device), switch_tensor.float().to(device), is_start.float().to(device), targets.float().to(device)
                outputs, _, _, is_switch2 = model(outcome_tensor, is_start)

                # # 构造 ground truth: 当前 step 是否为 [0, 1]
                # switch_is_01 = (switch_tensor[:, :, 1] == 1).float().unsqueeze(-1)  # shape: (B, T, 1)
                # # loss_s1：模型对switch==[0,1]的预测准确性（二分类）
                # loss_s1 = nn.BCEWithLogitsLoss()(is_switch2, switch_is_01)

                batch_loss = 0.0
                loss_count = 0

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
                            weight = 6.0 # 6.0
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
                print(f"Early stopping at epoch {epoch} (no improvement in {early_stop_patience} epochs)")
                break


if __name__ == "__main__":
    set_seed(42)

    # ===== 1) 按被试读取（保持被试顺序 = mat 里 infer_data 的顺序）=====
    tdev_list, tf_list, target_list, is_start_list = extract_behavioral_data_subjectwise(
        '../data/infer_data.mat',
        row1='tDev',
        row2='TF',
        row3='pr_of_switch',   # 如果你在训 confidence，就改成 'mu_switch_estimated'
        var_name='infer_data'
    )
    n_sub = len(tdev_list)
    print("N subjects =", n_sub)

    # 用来存每个被试的 held-out 预测与 hidden（LOSO）
    pred_by_sub = [None] * n_sub          # each: (T,)
    hidden_by_sub = [None] * n_sub        # each: (T,H)

    # ===== 2) LOSO 训练与测试 =====
    for test_idx in range(n_sub):
        print("\n" + "="*60)
        print(f"[LOSO] Test subject = {test_idx}/{n_sub-1}")

        # train pool = all but test
        train_pool = [i for i in range(n_sub) if i != test_idx]

        # 最小且必要：从 train_pool 挑 1 个做 valid，避免 valid=train
        valid_idx = train_pool[0]
        train_idx = train_pool[1:] if len(train_pool) > 1 else train_pool

        dataset_train = RNNInputDatasetSubjectwise(
            [tdev_list[i] for i in train_idx],
            [tf_list[i] for i in train_idx],
            [target_list[i] for i in train_idx],
            [is_start_list[i] for i in train_idx]
        )
        dataset_valid = RNNInputDatasetSubjectwise(
            [tdev_list[valid_idx]],
            [tf_list[valid_idx]],
            [target_list[valid_idx]],
            [is_start_list[valid_idx]]
        )
        dataset_test = RNNInputDatasetSubjectwise(
            [tdev_list[test_idx]],
            [tf_list[test_idx]],
            [target_list[test_idx]],
            [is_start_list[test_idx]]
        )

        # 每折重新初始化模型（LOSO 标准）
        model = RNN(is_sigmoid=1, hidden_size=32)

        fold_model_name = f"best_model_LOSO_testsub{test_idx}.pt"

        # 训练
        train_val_loop(
            model=model,
            dataset_train=dataset_train,
            dataset_valid=dataset_valid,
            save_model_path="../model",
            early_stop_patience=10,
            model_name=fold_model_name,
            epochs=10000,
            batch_size=1,
            shuffle=False
        )

        # ===== 3) 用该折最优模型对 held-out 被试做预测，取 hidden state =====
        best_path = os.path.join("../model", fold_model_name)
        model.load_state_dict(torch.load(best_path, map_location="cpu"))
        model.eval()

        loader = DataLoader(dataset_test, batch_size=1, shuffle=False)

        with torch.no_grad():
            for outcome_tensor, switch_tensor, is_start_tensor, targets in loader:
                outcome_tensor = outcome_tensor.float()
                is_start_tensor = is_start_tensor.float()

                pred, hidden_states, _, _ = model(outcome_tensor, is_start_tensor)
                # pred: (1,T,1) -> (T,)
                pred_np = pred.squeeze(0).squeeze(-1).cpu().numpy()
                # hidden_states(out): (1,T,H) -> (T,H)
                hid_np = hidden_states.squeeze(0).cpu().numpy()

                pred_by_sub[test_idx] = pred_np
                hidden_by_sub[test_idx] = hid_np

        print(f"[LOSO] saved: pred {pred_by_sub[test_idx].shape}, hidden {hidden_by_sub[test_idx].shape}")

    # ===== 4) 关键：拼回“与原 mat 顺序一致”的长向量/矩阵，然后按原函数保存到同名 mat =====
    # 你的原始拼接顺序是 for subj in infer_data 顺序逐个 extend；这里同样按 0..n_sub-1 拼接即可
    pred_all = np.concatenate(pred_by_sub, axis=0)                     # (total_trials,)
    hidden_all = np.concatenate(hidden_by_sub, axis=0)                 # (total_trials, H)

    # ===== LOSO overall performance (held-out) =====
    target_all = np.concatenate(target_list, axis=0)   # 与 pred_all 顺序一致（按被试顺序拼接）
    assert pred_all.shape[0] == target_all.shape[0]

    # 1) Overall R² / MSE
    r2_all = r2_score(target_all, pred_all)
    mse_all = mean_squared_error(target_all, pred_all)
    rmse_all = np.sqrt(mse_all)

    print("\n" + "="*60)
    print("[LOSO Overall] (all held-out trials concatenated)")
    print(f"R2   = {r2_all:.4f}")
    print(f"MSE  = {mse_all:.6f}")
    print(f"RMSE = {rmse_all:.6f}")

    # 2) Per-subject R² (optional but recommended)
    r2_by_sub = []
    for si in range(n_sub):
        y_true = np.asarray(target_list[si]).reshape(-1)
        y_pred = np.asarray(pred_by_sub[si]).reshape(-1)
        # 如果某个被试序列太短，R²可能不稳定/报错，这里做个保护
        if len(y_true) < 2 or np.allclose(np.var(y_true), 0):
            r2 = np.nan
        else:
            r2 = r2_score(y_true, y_pred)
        r2_by_sub.append(r2)

    r2_by_sub = np.array(r2_by_sub, dtype=float)
    print("\n[LOSO Per-subject R2]")
    print("mean±std (ignoring NaN) =", np.nanmean(r2_by_sub), "±", np.nanstd(r2_by_sub))
    print("R2_by_sub =", r2_by_sub)

    is_train_all = np.full(pred_all.shape[0], 2, dtype=int)   # LOSO: 每个 trial 都是 held-out

    # 注意：不改变你保存文件的形式/名称（沿用你原来的 save_path）
    pre_save_in_mat(
        matlab_path="../data/infer_data.mat",
        pre=pred_all,
        hidden_states=hidden_all,
        save_path="../data/update_data_switch_weight_diff0.mat",
        is_train_array=is_train_all   # <<< ADD THIS
    )

    print("\n[DONE] Saved to ../data/update_data_switch_weight_diff0.mat (same format/name as before)")